from __future__ import annotations

from pathlib import Path
import os
import shutil
import socket
import tempfile
import pandas as pd
import shlex

from .utils import ensure_dir, run_cmd, is_nonempty_file, step


def _scratch_root(args) -> Path:
    """Use explicit or project-backed scratch; never silently default to /tmp."""
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
            f"{min_free_gb:.1f} GiB required.",
            str(root),
        )
    return root


def _sample_scratch(args, stage: str, sid: str) -> Path:
    safe_sid = "".join(c if c.isalnum() or c in "._-" else "_" for c in sid)
    job = os.environ.get("SLURM_JOB_ID", "nojid")
    host = socket.gethostname().split(".")[0]
    prefix = f"dsRNASeeker_{stage}_{safe_sid}_{job}_{host}_{os.getpid()}_"
    return Path(tempfile.mkdtemp(prefix=prefix, dir=str(_scratch_root(args))))


def _move_replace(src: Path, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists():
        dest.unlink()
    shutil.move(str(src), str(dest))


def _all_filtered_exist(df: pd.DataFrame, outdir: Path) -> bool:
    return all(is_nonempty_file(outdir / f"{sid}_filtered_editing_events.txt") for sid in df["sample_id"].astype(str))


def run_reditools2(args, bam_samplesheet: str | Path, outdir: str | Path) -> Path:
    """Run a generic two-stage REDItools2 module.

    Stage 1: REDItools2/editing caller from BAM+FASTA -> raw per-sample table.
    Stage 2: generic A-to-I QC filter -> <sample>_filtered_editing_events.txt.

    No dataset-specific sample IDs or FCCC paths are hard-coded; all samples come
    from the BAM samplesheet created by the workflow.
    """
    outdir = ensure_dir(outdir)
    raw_dir = ensure_dir(outdir / "raw")
    filt_dir = ensure_dir(outdir)
    df = pd.read_csv(bam_samplesheet, sep="\t")

    if _all_filtered_exist(df, filt_dir) and not getattr(args, "force", False):
        step(f"Step 4/6 RNA editing: reusing REDItools2 filtered outputs in {filt_dir}")
        return filt_dir

    summaries = []
    for row in df.itertuples(index=False):
        sid = str(row.sample_id)
        cond = str(row.condition)
        bam = str(row.bam_path)
        raw = raw_dir / f"{sid}_reditools2_raw.txt"
        filt = filt_dir / f"{sid}_filtered_editing_events.txt"
        sample_summary = outdir / "qc" / f"{sid}_reditools2_summary.tsv"
        log = Path(args.output_dir) / "pipeline_info" / "logs" / f"REDItools2.{sid}.log"

        # If the final filtered file is missing, regenerate the raw file even if
        # a previous nonempty raw file exists. REDItools2 can leave a nonempty but
        # truncated raw table after OSError / filesystem failures.
        need_raw = (
            getattr(args, "force", False)
            or not is_nonempty_file(raw)
            or not is_nonempty_file(filt)
        )

        if need_raw:
            if raw.exists():
                raw.unlink()
            sample_tmp = _sample_scratch(args, "REDItools2", sid)
            raw_tmp = sample_tmp / f"{sid}_reditools2_raw.txt"
            # REDItools2 installations differ. This default matches the REDItools2
            # src/cineca/reditools.py convention used in your legacy scripts:
            #   reditools.py -f input.bam -r reference.fa -o output.txt -s <strand>
            cmd = [
                args.reditools_exe,
                "-f", bam,
                "-r", str(args.fasta),
                "-o", str(raw_tmp),
                "-s", str(args.reditools_strand),
            ]
            if args.reditools_extra:
                cmd += shlex.split(args.reditools_extra)
            step(f"Step 4/6 RNA editing: running REDItools2 for {sid} using scratch {sample_tmp}")
            run_cmd(cmd, log_path=log, quiet=getattr(args, "quiet", True))
            if not is_nonempty_file(raw_tmp):
                raise RuntimeError(f"REDItools2 did not produce a nonempty raw file for {sid}: {raw_tmp}")
            _move_replace(raw_tmp, raw)
            shutil.rmtree(sample_tmp, ignore_errors=True)
        else:
            step(f"Step 4/6 RNA editing: reusing REDItools2 raw output for {sid}")

        if is_nonempty_file(filt) and not getattr(args, "force", False):
            step(f"Step 4/6 RNA editing: reusing REDItools2 filtered output for {sid}")
        else:
            rscript = Path(args.reditools_post_rscript) if args.reditools_post_rscript else Path(__file__).resolve().parents[1] / "r" / "reditools_filter_a2i.R"
            cmd = [
                args.rscript_exe, str(rscript),
                "--raw", str(raw),
                "--out", str(filt),
                "--sample", sid,
                "--condition", cond,
                "--strandedness", str(args.strandedness),
                "--min-meanq", str(args.reditools_min_meanq),
                "--min-coverage", str(args.reditools_min_coverage),
                "--min-frequency", str(args.reditools_min_frequency),
                "--summary-out", str(sample_summary),
            ]
            step(f"Step 4/6 RNA editing: filtering A-to-I REDItools2 events for {sid}")
            run_cmd(cmd, log_path=Path(args.output_dir) / "pipeline_info" / "logs" / f"REDItools2_filter.{sid}.log", quiet=getattr(args, "quiet", True))
        if sample_summary.exists():
            summaries.append(pd.read_csv(sample_summary, sep="\t"))

    if summaries:
        pd.concat(summaries, ignore_index=True).to_csv(outdir / "REDItools2_global_editing_index_summary.tsv", sep="\t", index=False)
    return filt_dir
