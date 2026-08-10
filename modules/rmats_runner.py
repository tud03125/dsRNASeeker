from __future__ import annotations

# DSRNASEEKER_RMATS_AUTO_MEL_FIX_V2

from pathlib import Path
import gzip
import json

import pandas as pd

from .utils import ensure_dir, run_cmd, is_nonempty_file, step


def _libtype_from_strandedness(s: str) -> str:
    s = (s or "auto").lower()
    if s in {"reverse", "fr-firststrand", "firststrand"}:
        return "fr-firststrand"
    if s in {"forward", "fr-secondstrand", "secondstrand"}:
        return "fr-secondstrand"
    return "fr-unstranded"


def _open_gtf_text(path: Path):
    """Open plain-text or gzip-compressed GTF without loading it into memory."""
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")


def _scan_gtf_max_exon_length(gtf: str | Path) -> dict:
    """Return the maximum annotated exon length and its source record.

    rMATS --mel means *maximum exon length*. GTF start/end coordinates are
    interpreted as inclusive, so exon length is end - start + 1.
    """
    path = Path(gtf)
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"GTF is missing or empty: {path}")

    max_length = 0
    max_record: dict[str, object] | None = None
    exon_count = 0
    malformed_exons = 0

    with _open_gtf_text(path) as handle:
        for line_number, raw in enumerate(handle, start=1):
            if not raw or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n\r").split("\t")
            if len(fields) < 9 or fields[2].lower() != "exon":
                continue
            try:
                start = int(fields[3])
                end = int(fields[4])
            except ValueError:
                malformed_exons += 1
                continue
            if start < 1 or end < start:
                malformed_exons += 1
                continue
            exon_count += 1
            length = end - start + 1
            if length > max_length:
                max_length = length
                max_record = {
                    "seqname": fields[0],
                    "start": start,
                    "end": end,
                    "length": length,
                    "strand": fields[6],
                    "attributes": fields[8],
                    "line_number": line_number,
                }

    if exon_count == 0 or max_record is None:
        raise ValueError(
            f"Could not calculate rMATS --mel: no valid exon records were found in {path}"
        )

    resolved = max(500, max_length)
    return {
        "gtf": str(path.resolve()),
        "gtf_size_bytes": path.stat().st_size,
        "gtf_mtime_ns": path.stat().st_mtime_ns,
        "valid_exon_records": exon_count,
        "malformed_exon_records_skipped": malformed_exons,
        "maximum_annotated_exon_length": max_length,
        "resolved_mel": resolved,
        "maximum_exon_record": max_record,
    }


def _resolve_positive_int(value: object, option: str) -> int:
    try:
        parsed = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{option} must be a positive integer; received {value!r}") from exc
    if parsed < 1:
        raise ValueError(f"{option} must be >= 1; received {parsed}")
    return parsed


def _resolve_rmats_novel_ss_limits(args, rmats_outdir: Path) -> tuple[int, int]:
    mil = _resolve_positive_int(
        getattr(args, "rmats_min_intron_length", 1),
        "--rmats-min-intron-length/--mil",
    )

    requested_mel = getattr(args, "rmats_max_exon_length", "auto")
    auto = requested_mel is None or str(requested_mel).strip().lower() == "auto"
    if auto:
        scan = _scan_gtf_max_exon_length(args.gtf)
        mel = int(scan["resolved_mel"])
        source = "auto_gtf_maximum_annotated_exon_length"
    else:
        mel = _resolve_positive_int(
            requested_mel,
            "--rmats-max-exon-length/--mel",
        )
        scan = None
        source = "explicit_user_value"

    provenance = {
        "novel_splice_sites_enabled": bool(getattr(args, "rmats_novel_ss", False)),
        "mil": mil,
        "mil_source": "pipeline_default_or_explicit_user_value",
        "mel": mel,
        "mel_source": source,
        "note": "rMATS --mel is maximum exon length, not maximum intron length",
        "gtf_scan": scan,
    }
    provenance_path = Path(rmats_outdir) / "rmats_novel_ss_limits.json"
    provenance_path.write_text(json.dumps(provenance, indent=2) + "\n")

    if auto:
        step(
            "Step 3/6 splicing: auto rMATS novelSS limits "
            f"--mil {mil} --mel {mel} from maximum annotated exon length in {args.gtf}"
        )
    else:
        step(
            "Step 3/6 splicing: explicit rMATS novelSS limits "
            f"--mil {mil} --mel {mel}"
        )
    return mil, mel


def run_rmats_case_control(args, bam_samplesheet: str | Path, rmats_outdir: str | Path) -> Path:
    rmats_outdir = ensure_dir(rmats_outdir)
    tmpdir = ensure_dir(Path(rmats_outdir) / "tmp")
    df = pd.read_csv(bam_samplesheet, sep="\t")
    case_bams = df.loc[df["condition"].astype(str) == args.case_label, "bam_path"].astype(str).tolist()
    ctrl_bams = df.loc[df["condition"].astype(str) == args.control_label, "bam_path"].astype(str).tolist()
    if not case_bams or not ctrl_bams:
        raise ValueError("rMATS requires at least one case and one control BAM")
    b1 = Path(rmats_outdir) / "b1_case.txt"
    b2 = Path(rmats_outdir) / "b2_control.txt"
    b1.write_text(",".join(case_bams) + "\n")
    b2.write_text(",".join(ctrl_bams) + "\n")
    out_file = Path(rmats_outdir) / f"RI.MATS.{args.rmats_track}.txt"
    if is_nonempty_file(out_file) and not getattr(args, "force", False):
        step(f"Step 3/6 splicing: reusing rMATS output {out_file}")
        return Path(rmats_outdir)
    libtype = args.rmats_libtype or _libtype_from_strandedness(args.strandedness)
    mil, mel = _resolve_rmats_novel_ss_limits(args, Path(rmats_outdir))
    cmd = [
        args.rmats_exe,
        "--b1", str(b1),
        "--b2", str(b2),
        "--gtf", str(args.gtf),
        "-t", "paired" if args.paired else "single",
        "--libType", libtype,
        "--readLength", str(args.read_length),
        "--variable-read-length",
        "--allow-clipping",
        "--cstat", str(args.rmats_cstat),
        "--mil", str(mil),
        "--mel", str(mel),
        "--task", "both",
        "--nthread", str(args.threads),
        "--tstat", str(args.threads),
        "--od", str(rmats_outdir),
        "--tmp", str(tmpdir),
    ]
    if args.rmats_novel_ss:
        cmd.append("--novelSS")
    step("Step 3/6 splicing: running rMATS with b1=case and b2=control")
    run_cmd(
        cmd,
        log_path=Path(args.output_dir) / "pipeline_info" / "logs" / "rMATS.log",
        quiet=getattr(args, "quiet", True),
    )
    return Path(rmats_outdir)
