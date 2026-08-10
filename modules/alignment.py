from __future__ import annotations

import os
import fcntl
import shutil
import socket
import tempfile
import uuid
from pathlib import Path
import pandas as pd

from .utils import ensure_dir, run_cmd, is_nonempty_file, step


def _star_index_ready(index_dir: Path) -> bool:
    """Return True only for a nonempty, apparently complete STAR index."""
    required = ["Genome", "SA", "SAindex", "genomeParameters.txt"]
    return (
        index_dir.is_dir()
        and all(
            (index_dir / name).is_file() and (index_dir / name).stat().st_size > 0
            for name in required
        )
    )


def _scratch_root(args) -> Path:
    """Return a writable scratch root for large transient BAM files.

    Large STAR and samtools intermediates must not silently fall back to /tmp.
    Priority:
      1. DSRNASEEKER_SCRATCH, when explicitly supplied by the user/SBATCH.
      2. SLURM_TMPDIR only when DSRNASEEKER_USE_NODE_SCRATCH=1.
      3. <output_dir>/00_scratch on the same project filesystem as the run.

    DSRNASEEKER_MIN_SCRATCH_GB controls the minimum free-space preflight
    (default: 20 GiB).
    """
    explicit = os.environ.get("DSRNASEEKER_SCRATCH")
    use_node = os.environ.get("DSRNASEEKER_USE_NODE_SCRATCH", "0") == "1"

    if explicit:
        root = Path(explicit)
    elif use_node and os.environ.get("SLURM_TMPDIR"):
        root = Path(os.environ["SLURM_TMPDIR"])
    else:
        root = Path(args.output_dir) / "00_scratch"

    root.mkdir(parents=True, exist_ok=True)
    if not os.access(root, os.W_OK | os.X_OK):
        raise PermissionError(f"Scratch directory is not writable: {root}")

    min_free_gb = float(os.environ.get("DSRNASEEKER_MIN_SCRATCH_GB", "20"))
    free_gb = shutil.disk_usage(root).free / (1024 ** 3)
    if free_gb < min_free_gb:
        raise OSError(
            28,
            f"Insufficient scratch space: {free_gb:.1f} GiB available at {root}; "
            f"{min_free_gb:.1f} GiB required. Set DSRNASEEKER_SCRATCH to a "
            "large project-backed directory.",
            str(root),
        )
    return root


def _sample_scratch(args, stage: str, sid: str) -> Path:
    """Create a collision-proof per-sample scratch directory.

    Earlier versions used only sample_id + PID. On shared filesystems, PIDs can
    collide across compute nodes or concurrent jobs, and one process could remove
    another process's active scratch directory. That produced intermittent
    samtools errors such as "failed writing ... No such file or directory".

    tempfile.mkdtemp creates a unique directory atomically. We never delete a
    pre-existing path here.
    """
    safe_sid = "".join(c if c.isalnum() or c in "._-" else "_" for c in sid)
    job = os.environ.get("SLURM_JOB_ID", "nojid")
    host = socket.gethostname().split(".")[0]
    root = _scratch_root(args)
    prefix = f"dsRNASeeker_{stage}_{safe_sid}_{job}_{host}_{os.getpid()}_"
    return Path(tempfile.mkdtemp(prefix=prefix, dir=str(root)))


def _assert_dir_exists(path: Path, context: str) -> None:
    if not path.exists():
        raise FileNotFoundError(
            f"Scratch directory disappeared during {context}: {path}. "
            "Use a non-purged scratch location, e.g. export DSRNASEEKER_SCRATCH=/rs01/projects/jadezhoulab/tud03125/dsRNASeeker_scratch/$SLURM_JOB_ID, "
            "or avoid running multiple workflows with the same shared scratch root."
        )


def _move_replace(src: Path, dest: Path) -> None:
    """Move src to dest, replacing incomplete old dest if present.

    shutil.move works across filesystems; Path.replace/os.replace do not.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists():
        dest.unlink()
    shutil.move(str(src), str(dest))


def build_star_index(args, index_dir: str | Path) -> Path:
    """Build or reuse a STAR index without shared ``./_STARtmp`` collisions.

    STAR's default genome-generation temporary directory is relative to the
    process working directory. Concurrent workflows launched from the same
    repository can therefore contend for ``./_STARtmp`` even when their final
    genome directories differ. This implementation gives every build a unique
    project-backed work directory, explicitly supplies ``--outTmpDir``, and
    promotes the completed index only after sentinel-file validation.

    A per-index advisory lock also prevents two jobs from constructing the same
    shared ``--star-index`` concurrently.
    """
    index_dir = Path(index_dir)
    index_dir.parent.mkdir(parents=True, exist_ok=True)

    if _star_index_ready(index_dir) and not getattr(args, "force", False):
        step(f"Step 1a/6 STAR index: reusing existing index at {index_dir}")
        return index_dir

    sjdb = (
        int(args.sjdb_overhang)
        if args.sjdb_overhang is not None
        else max(int(args.read_length or 101) - 1, 1)
    )
    log = Path(args.output_dir) / "pipeline_info" / "logs" / "STAR_genomeGenerate.log"
    lock_path = index_dir.parent / f".{index_dir.name}.dsRNASeeker_build.lock"

    # flock() protects the open lock-file inode and releases automatically if a
    # job exits. The tiny lock file may remain safely on disk.
    with lock_path.open("a+", encoding="utf-8") as lock_handle:
        step(f"Step 1a/6 STAR index: waiting for build lock {lock_path}")
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)

        # Another process may have completed the index while this job waited.
        if _star_index_ready(index_dir) and not getattr(args, "force", False):
            step(f"Step 1a/6 STAR index: reusing completed index at {index_dir}")
            return index_dir

        build_root = _sample_scratch(args, "STARindex", index_dir.name)
        staging_index = build_root / "genome_index"
        star_tmp = build_root / "STARtmp"
        staging_index.mkdir(parents=True, exist_ok=True)
        # Do not pre-create star_tmp: STAR expects to create --outTmpDir itself.

        cmd = [
            args.star_exe,
            "--runThreadN", str(args.threads),
            "--runMode", "genomeGenerate",
            "--genomeDir", str(staging_index),
            "--genomeFastaFiles", str(args.fasta),
            "--sjdbGTFfile", str(args.gtf),
            "--sjdbOverhang", str(sjdb),
            "--outTmpDir", str(star_tmp),
        ]

        step(
            "Step 1a/6 STAR index: building "
            f"{index_dir} with unique scratch {build_root}"
        )
        try:
            run_cmd(cmd, log_path=log, quiet=getattr(args, "quiet", True))

            if not _star_index_ready(staging_index):
                raise RuntimeError(
                    "STAR genomeGenerate exited without creating a complete "
                    f"index in {staging_index}. See log: {log}"
                )

            # Remove an incomplete or forced old index only after a complete
            # replacement has been generated successfully.
            if index_dir.exists():
                if index_dir.is_dir():
                    shutil.rmtree(index_dir)
                else:
                    index_dir.unlink()

            shutil.move(str(staging_index), str(index_dir))

            if not _star_index_ready(index_dir):
                raise RuntimeError(
                    f"STAR index promotion failed validation: {index_dir}"
                )
        finally:
            shutil.rmtree(build_root, ignore_errors=True)

    return index_dir


def run_star_alignment(args, samples: pd.DataFrame, index_dir: str | Path, star_outdir: str | Path) -> dict[str, Path]:
    """Run STAR and return coordinate-sorted BAMs.

    Robust mode used here:
      STAR writes BAM Unsorted to per-sample scratch, not project storage.
      samtools sort then creates the coordinate-sorted BAM in scratch with -T.
      Only validated final files are moved back to star_outdir.

    This avoids STAR's internal BAM SortedByCoordinate step, which can fail on
    large samples/project filesystems with: "number of bytes expected from the
    BAM bin does not agree with the actual size on disk".
    """
    star_outdir = ensure_dir(star_outdir)
    out: dict[str, Path] = {}
    for row in samples.itertuples(index=False):
        sid = str(row.sample_id)
        fq1 = str(row.fastq_1)
        fq2 = str(getattr(row, "fastq_2", "") or "")

        final_bam = star_outdir / f"{sid}_Aligned.sortedByCoord.out.bam"
        if is_nonempty_file(final_bam) and not getattr(args, "force", False):
            step(f"Step 1b/6 STAR alignment: reusing {sid}")
            out[sid] = final_bam
            continue

        # Remove stale failed output, e.g. 0-byte BAM from interrupted STAR sort.
        if final_bam.exists() and not is_nonempty_file(final_bam):
            final_bam.unlink()

        log = Path(args.output_dir) / "pipeline_info" / "logs" / f"STAR_align.{sid}.log"
        sample_tmp = _sample_scratch(args, "STAR", sid)
        sample_prefix_tmp = sample_tmp / f"{sid}_"
        unsorted_bam_tmp = sample_tmp / f"{sid}_Aligned.out.bam"
        sorted_bam_tmp = sample_tmp / f"{sid}_Aligned.sortedByCoord.out.bam"

        cmd = [
            args.star_exe,
            "--runThreadN", str(args.threads),
            "--genomeDir", str(index_dir),
            "--readFilesIn", fq1,
        ]
        if fq2:
            cmd.append(fq2)
        if fq1.endswith((".gz", ".gzip")) or fq2.endswith((".gz", ".gzip")):
            cmd += ["--readFilesCommand", "zcat"]
        cmd += [
            "--outFileNamePrefix", str(sample_prefix_tmp),
            "--outSAMtype", "BAM", "Unsorted",
            "--outSAMstrandField", "intronMotif",
            "--quantMode", "GeneCounts",
        ]

        step(f"Step 1b/6 STAR alignment: running {sid} using scratch {sample_tmp}")
        run_cmd(cmd, log_path=log, quiet=getattr(args, "quiet", True))

        if not unsorted_bam_tmp.exists():
            raise FileNotFoundError(f"STAR did not produce expected unsorted BAM: {unsorted_bam_tmp}. See log: {log}")

        # STAR/HTSlib can emit BGZF write errors yet still reach its normal
        # "finished successfully" message. Validate the unsorted BAM before
        # attempting a costly sort so a full or failed filesystem is reported
        # immediately and cannot propagate a truncated BAM downstream.
        run_cmd(
            [args.samtools_exe, "quickcheck", "-v", str(unsorted_bam_tmp)],
            log_path=log,
            quiet=getattr(args, "quiet", True),
        )

        # Sort outside STAR, with an explicit temporary-file prefix in scratch.
        _assert_dir_exists(sample_tmp, f"STAR/samtools-sort for {sid}")
        sort_tmp_prefix = sample_tmp / f"{sid}.coord_sort_tmp"
        run_cmd([
            args.samtools_exe, "sort",
            "-@", str(args.threads),
            "-m", "1G",
            "-T", str(sort_tmp_prefix),
            "-o", str(sorted_bam_tmp),
            str(unsorted_bam_tmp),
        ], log_path=log, quiet=getattr(args, "quiet", True))
        run_cmd([args.samtools_exe, "quickcheck", "-v", str(sorted_bam_tmp)], log_path=log, quiet=getattr(args, "quiet", True))

        # Move useful STAR text outputs back under the normal prefix.
        for suffix in [
            "Log.out",
            "Log.progress.out",
            "Log.final.out",
            "ReadsPerGene.out.tab",
            "SJ.out.tab",
        ]:
            src = sample_tmp / f"{sid}_{suffix}"
            if src.exists():
                _move_replace(src, star_outdir / f"{sid}_{suffix}")

        _assert_dir_exists(sample_tmp, f"moving STAR outputs for {sid}")
        _move_replace(sorted_bam_tmp, final_bam)
        run_cmd([args.samtools_exe, "quickcheck", "-v", str(final_bam)], log_path=log, quiet=getattr(args, "quiet", True))

        shutil.rmtree(sample_tmp, ignore_errors=True)
        out[sid] = final_bam
    return out


def markdup_bams(args, aligned_bams: dict[str, Path], markdup_outdir: str | Path) -> dict[str, Path]:
    """Mark duplicates with all large intermediates in scratch.

    This avoids writing name-sorted/fixmate/coord-sorted intermediate BAMs
    directly to the project filesystem. Only the final validated BAM and index
    are moved back to markdup_outdir.
    """
    markdup_outdir = ensure_dir(markdup_outdir)
    out: dict[str, Path] = {}
    for sid, bam in aligned_bams.items():
        md_bam = markdup_outdir / f"{sid}.markdup.sorted.bam"
        md_bai = Path(str(md_bam) + ".bai")
        if is_nonempty_file(md_bam) and is_nonempty_file(md_bai) and not getattr(args, "force", False):
            step(f"Step 1c/6 mark duplicates: reusing {sid}")
            out[sid] = md_bam
            continue

        # Remove stale failed final outputs only. Old intermediates from earlier
        # pipeline versions can be ignored or manually archived.
        for stale in [md_bam, md_bai]:
            if stale.exists() and not is_nonempty_file(stale):
                stale.unlink()

        log = Path(args.output_dir) / "pipeline_info" / "logs" / f"samtools_markdup.{sid}.log"
        sample_tmp = _sample_scratch(args, "markdup", sid)
        sorted_bam = sample_tmp / f"{sid}.name_sorted.bam"
        fixmate_bam = sample_tmp / f"{sid}.fixmate.bam"
        coord_bam = sample_tmp / f"{sid}.coord_sorted.bam"
        md_bam_tmp = sample_tmp / f"{sid}.markdup.sorted.bam"
        md_bai_tmp = Path(str(md_bam_tmp) + ".bai")

        step(f"Step 1c/6 mark duplicates: running {sid} using scratch {sample_tmp}")
        _assert_dir_exists(sample_tmp, f"markdup name-sort for {sid}")
        run_cmd([
            args.samtools_exe, "sort",
            "-@", str(args.threads),
            "-m", "1G",
            "-n",
            "-T", str(sample_tmp / f"{sid}.name_sort_tmp"),
            "-o", str(sorted_bam),
            str(bam),
        ], log_path=log, quiet=getattr(args, "quiet", True))
        _assert_dir_exists(sample_tmp, f"markdup fixmate for {sid}")
        run_cmd([args.samtools_exe, "fixmate", "-@", str(args.threads), "-m", str(sorted_bam), str(fixmate_bam)], log_path=log, quiet=getattr(args, "quiet", True))
        _assert_dir_exists(sample_tmp, f"markdup coordinate-sort for {sid}")
        run_cmd([
            args.samtools_exe, "sort",
            "-@", str(args.threads),
            "-m", "1G",
            "-T", str(sample_tmp / f"{sid}.coord_sort_tmp"),
            "-o", str(coord_bam),
            str(fixmate_bam),
        ], log_path=log, quiet=getattr(args, "quiet", True))
        _assert_dir_exists(sample_tmp, f"samtools markdup for {sid}")
        run_cmd([args.samtools_exe, "markdup", "-@", str(args.threads), str(coord_bam), str(md_bam_tmp)], log_path=log, quiet=getattr(args, "quiet", True))
        run_cmd([args.samtools_exe, "index", "-@", str(args.threads), str(md_bam_tmp)], log_path=log, quiet=getattr(args, "quiet", True))
        run_cmd([args.samtools_exe, "quickcheck", "-v", str(md_bam_tmp)], log_path=log, quiet=getattr(args, "quiet", True))

        _move_replace(md_bam_tmp, md_bam)
        _move_replace(md_bai_tmp, md_bai)
        run_cmd([args.samtools_exe, "quickcheck", "-v", str(md_bam)], log_path=log, quiet=getattr(args, "quiet", True))

        shutil.rmtree(sample_tmp, ignore_errors=True)
        out[sid] = md_bam
    return out


def write_bam_samplesheet(samples: pd.DataFrame, bams: dict[str, Path], out_path: str | Path) -> Path:
    rows = []
    for row in samples.itertuples(index=False):
        sid = str(row.sample_id)
        if sid not in bams:
            raise ValueError(f"No BAM produced/found for sample {sid}")
        rows.append({"sample_id": sid, "condition": str(row.condition), "bam_path": str(bams[sid])})
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out_path, sep="\t", index=False)
    return out_path
