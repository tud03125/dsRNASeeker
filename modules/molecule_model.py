from __future__ import annotations

"""Conservative molecule-origin annotation for dsRNASeeker candidates.

This module does *not* claim to determine RNA molecular topology from short-read
RNA-seq.  It annotates whether an inverted TE pair is more consistent with:

* a putative intramolecular foldback model;
* a putative intermolecular sense-antisense model; or
* an unresolved inverted-TE-pair model.

The classification is deliberately rule based and auditable.  Same-transcript
support comes from a supplied GTF.  Stranded coverage and retained-intron
support are reused from the existing dsRNASeeker summary.  Optional Z-RNA
columns are carried into the output for interpretation but are not used to call
intramolecular versus intermolecular origin.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable
import argparse
import gzip
import re

import numpy as np
import pandas as pd


_EMPTY = {"", "nan", "none", "na", "n/a", "."}


@dataclass(frozen=True)
class Arm:
    te_id: str
    chrom: str
    start: int
    end: int
    repeat_strand: str

    @property
    def length(self) -> int:
        return max(1, self.end - self.start)


@dataclass(frozen=True)
class TranscriptSpan:
    transcript_id: str
    gene_id: str
    gene_name: str
    chrom: str
    start: int
    end: int
    strand: str


def _clean_text(value) -> str:
    if value is None or pd.isna(value):
        return ""
    text = str(value).strip()
    return "" if text.lower() in _EMPTY else text


def _as_bool(value) -> bool:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if value is None or pd.isna(value):
        return False
    return str(value).strip().lower() in {"1", "true", "t", "yes", "y"}


def _first_existing(row: pd.Series, names: Iterable[str], default=None):
    for name in names:
        if name in row.index:
            value = row.get(name)
            if value is not None and not pd.isna(value):
                return value
    return default


def _open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix.lower() == ".gz" else path.open("rt")


def _parse_gtf_attributes(text: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for token in text.strip().strip(";").split(";"):
        token = token.strip()
        if not token:
            continue
        if " " in token:
            key, value = token.split(None, 1)
            attrs[key] = value.strip().strip('"')
        elif "=" in token:
            key, value = token.split("=", 1)
            attrs[key.strip()] = value.strip().strip('"')
    return attrs


def _parse_te_id(te_id: str) -> Arm:
    """Parse RepeatMasker-style IDs such as B1_chr1_100_200_+.

    The split is made from the right so repeat names may contain underscores.
    Chromosome names are expected to contain a ``chr`` token, as in the
    standard mm10/hg38 outputs used by dsRNASeeker.
    """
    text = _clean_text(te_id)
    m = re.search(r"_(\d+)_(\d+)_([+-])$", text)
    if not m:
        raise ValueError(
            f"Cannot parse TE coordinates from {te_id!r}; expected a suffix "
            "like _chr1_100_200_+"
        )
    start, end, strand = int(m.group(1)), int(m.group(2)), m.group(3)
    prefix = text[: m.start()]
    chr_pos = prefix.rfind("_chr")
    if chr_pos < 0:
        raise ValueError(
            f"Cannot identify chromosome in {te_id!r}; expected a '_chr...' token"
        )
    chrom = prefix[chr_pos + 1 :]
    if end < start:
        start, end = end, start
    return Arm(text, chrom, start, end, strand)


def _read_summary(output_dir: Path, subset: str, summary_in: str | None) -> tuple[Path, pd.DataFrame]:
    if summary_in:
        path = Path(summary_in)
    else:
        path = output_dir / subset / "summary" / "TEpair_dsRNA_master.summary.with_RI.csv"
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing/empty dsRNASeeker summary: {path}")
    df = pd.read_csv(path, low_memory=False)
    if "pair_id" not in df.columns:
        raise ValueError(f"Summary lacks pair_id: {path}")
    return path, df


def _arm_from_row(row: pd.Series, side: str) -> Arm:
    col = f"{side}_TE_id"
    value = _clean_text(row.get(col, ""))
    if not value:
        pair_id = _clean_text(row.get("pair_id", ""))
        parts = pair_id.split("__", 1)
        if len(parts) != 2:
            raise ValueError(f"Cannot split pair_id into A/B TE IDs: {pair_id!r}")
        value = parts[0 if side == "A" else 1]
    return _parse_te_id(value)


def _collect_query_chromosomes(df: pd.DataFrame) -> set[str]:
    chroms: set[str] = set()
    for _, row in df.iterrows():
        for side in ("A", "B"):
            chroms.add(_arm_from_row(row, side).chrom)
    return chroms


def _read_transcript_spans(gtf: Path, keep_chroms: set[str]) -> dict[str, list[TranscriptSpan]]:
    """Read transcript spans, deriving them from exon rows when necessary."""
    if not gtf.is_file() or gtf.stat().st_size == 0:
        raise FileNotFoundError(f"Missing/empty GTF: {gtf}")

    agg: dict[tuple[str, str], dict[str, object]] = {}
    with _open_text(gtf) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom, _, feature, start1, end1, _, strand, _, attr_text = fields[:9]
            if chrom not in keep_chroms:
                continue
            if feature.lower() not in {"transcript", "mrna", "exon", "utr", "cds"}:
                continue
            attrs = _parse_gtf_attributes(attr_text)
            tx_id = _clean_text(
                attrs.get("transcript_id")
                or attrs.get("transcript")
                or attrs.get("Parent")
            )
            if not tx_id:
                continue
            gene_id = _clean_text(attrs.get("gene_id") or attrs.get("gene") or attrs.get("geneID"))
            gene_name = _clean_text(
                attrs.get("gene_name")
                or attrs.get("gene_symbol")
                or attrs.get("Name")
            )
            start0 = max(0, int(start1) - 1)
            end0 = int(end1)
            key = (chrom, tx_id)
            rec = agg.setdefault(
                key,
                {
                    "transcript_id": tx_id,
                    "gene_id": gene_id,
                    "gene_name": gene_name,
                    "chrom": chrom,
                    "start": start0,
                    "end": end0,
                    "strand": strand,
                },
            )
            rec["start"] = min(int(rec["start"]), start0)
            rec["end"] = max(int(rec["end"]), end0)
            if not rec["gene_id"] and gene_id:
                rec["gene_id"] = gene_id
            if not rec["gene_name"] and gene_name:
                rec["gene_name"] = gene_name
            if rec["strand"] not in {"+", "-"} and strand in {"+", "-"}:
                rec["strand"] = strand

    by_chrom: dict[str, list[TranscriptSpan]] = {chrom: [] for chrom in keep_chroms}
    for rec in agg.values():
        span = TranscriptSpan(
            transcript_id=str(rec["transcript_id"]),
            gene_id=str(rec["gene_id"]),
            gene_name=str(rec["gene_name"]),
            chrom=str(rec["chrom"]),
            start=int(rec["start"]),
            end=int(rec["end"]),
            strand=str(rec["strand"]),
        )
        by_chrom.setdefault(span.chrom, []).append(span)
    for chrom in by_chrom:
        by_chrom[chrom].sort(key=lambda x: (x.start, x.end, x.transcript_id))
    return by_chrom


def _overlap_len(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    return max(0, min(a_end, b_end) - max(a_start, b_start))


def _map_arm(
    arm: Arm,
    transcripts_by_chrom: dict[str, list[TranscriptSpan]],
    overlap_fraction: float,
    slop: int,
) -> list[TranscriptSpan]:
    hits: list[TranscriptSpan] = []
    query_start = max(0, arm.start - slop)
    query_end = arm.end + slop
    for tx in transcripts_by_chrom.get(arm.chrom, []):
        if tx.start >= query_end:
            break
        if tx.end <= query_start:
            continue
        overlap = _overlap_len(arm.start, arm.end, tx.start - slop, tx.end + slop)
        if overlap / arm.length >= overlap_fraction:
            hits.append(tx)
    return hits


def _join_limited(values: Iterable[str], limit: int) -> str:
    vals = sorted({_clean_text(v) for v in values if _clean_text(v)})
    if len(vals) <= limit:
        return ";".join(vals)
    return ";".join(vals[:limit]) + f";...(+{len(vals) - limit})"


def _annotation_category(row: pd.Series, side: str) -> str:
    return _clean_text(
        _first_existing(
            row,
            [f"{side}_annotation_category", f"{side}_annotation", f"{side}_annotation_from_TE"],
            "",
        )
    ).lower()


def _model_row(
    row: pd.Series,
    case_label: str,
    transcripts_by_chrom: dict[str, list[TranscriptSpan]],
    overlap_fraction: float,
    slop: int,
    max_ids: int,
) -> dict[str, object]:
    arm_a = _arm_from_row(row, "A")
    arm_b = _arm_from_row(row, "B")
    a_hits = _map_arm(arm_a, transcripts_by_chrom, overlap_fraction, slop)
    b_hits = _map_arm(arm_b, transcripts_by_chrom, overlap_fraction, slop)

    a_by_id = {x.transcript_id: x for x in a_hits}
    b_by_id = {x.transcript_id: x for x in b_hits}
    shared_ids = sorted(set(a_by_id).intersection(b_by_id))

    a_gene_ids = {x.gene_id for x in a_hits if x.gene_id}
    b_gene_ids = {x.gene_id for x in b_hits if x.gene_id}
    shared_gene_ids = sorted(a_gene_ids.intersection(b_gene_ids))
    a_gene_names = {x.gene_name for x in a_hits if x.gene_name}
    b_gene_names = {x.gene_name for x in b_hits if x.gene_name}
    shared_gene_names = sorted(a_gene_names.intersection(b_gene_names))

    summary_a_symbol = _clean_text(row.get("A_SYMBOL", ""))
    summary_b_symbol = _clean_text(row.get("B_SYMBOL", ""))
    same_summary_symbol = bool(
        summary_a_symbol
        and summary_b_symbol
        and summary_a_symbol.upper() == summary_b_symbol.upper()
    )
    same_gene_support = bool(shared_gene_ids or shared_gene_names or same_summary_symbol)
    distinct_summary_genes = bool(
        summary_a_symbol
        and summary_b_symbol
        and summary_a_symbol.upper() != summary_b_symbol.upper()
    )

    a_strands = {x.strand for x in a_hits if x.strand in {"+", "-"}}
    b_strands = {x.strand for x in b_hits if x.strand in {"+", "-"}}
    opposite_tx_strands = any(a != b for a in a_strands for b in b_strands)

    genomic_orientation = _clean_text(row.get("genomic_orientation", "")).lower()
    transcript_orientation = _clean_text(row.get("transcript_orientation", "")).lower()
    inverted_compatible = genomic_orientation == "inverted" or transcript_orientation == "inverted"
    direct_control = genomic_orientation == "direct" and transcript_orientation == "direct"

    case_both_strands = _as_bool(row.get(f"{case_label}_both_strands", False))
    case_arm_opposite = _as_bool(row.get(f"{case_label}_arm_opposite", False))
    case_arms_both_cov = _as_bool(row.get(f"{case_label}_arms_both_cov", False))

    ri_exact_both = _as_bool(row.get("RI_exact_overlap_both_arms", False))
    ri_overlap_both = _as_bool(row.get("RI_overlap_both_arms", False))
    ri_gene_match = _as_bool(row.get("RI_gene_match_any", False))
    ri_support = ri_exact_both or ri_overlap_both or ri_gene_match

    a_cat = _annotation_category(row, "A")
    b_cat = _annotation_category(row, "B")
    transcript_body_categories = {"intron", "utr3", "utr5", "exon", "cds"}
    both_transcript_body = a_cat in transcript_body_categories and b_cat in transcript_body_categories

    intra_flags = {
        "shared_transcript": bool(shared_ids),
        "same_gene": same_gene_support,
        "case_both_arms_covered": case_arms_both_cov,
        "RI_or_gene_match_support": ri_support,
        "both_arms_transcript_body_annotation": both_transcript_body,
    }
    inter_flags = {
        "case_both_strands_expressed": case_both_strands,
        "case_arm_dominant_strands_opposite": case_arm_opposite,
        "case_both_arms_covered": case_arms_both_cov,
        "opposite_annotated_transcript_strands": opposite_tx_strands,
        "distinct_summary_gene_symbols": distinct_summary_genes,
    }
    intra_count = sum(intra_flags.values())
    inter_count = sum(inter_flags.values())

    conflict = bool(shared_ids and case_both_strands and case_arm_opposite and case_arms_both_cov)

    if direct_control or not inverted_compatible:
        model = "direct_orientation_control"
        confidence = "not_applicable"
        basis = "pair does not satisfy inverted-orientation compatibility"
    elif shared_ids:
        model = "putative_intramolecular_foldback"
        if case_arms_both_cov and (ri_exact_both or ri_overlap_both):
            confidence = "high"
        elif case_arms_both_cov or ri_support or same_gene_support:
            confidence = "moderate"
        else:
            confidence = "low"
        basis = "both TE arms map to at least one shared annotated transcript"
        if conflict:
            basis += "; opposite-strand coverage also present, so topology remains experimentally unproven"
    elif case_arm_opposite and case_arms_both_cov and case_both_strands and opposite_tx_strands:
        model = "putative_intermolecular_sense_antisense"
        confidence = "high"
        basis = "no shared transcript; both arms covered; both strands expressed; arm-dominant and annotated transcript strands are opposite"
    elif (
        case_arm_opposite
        and case_arms_both_cov
        and case_both_strands
        and bool(a_hits)
        and bool(b_hits)
        and distinct_summary_genes
    ):
        model = "putative_intermolecular_sense_antisense"
        confidence = "moderate"
        basis = "no shared transcript; both arms map to transcripts, both arms are covered, and opposite-strand expression is present"
    else:
        model = "ambiguous_inverted_TE_pair"
        confidence = "unresolved"
        basis = "short-read and annotation evidence do not distinguish same-molecule from two-molecule origin"

    return {
        "pair_id": row.get("pair_id", ""),
        "molecule_model": model,
        "molecule_model_confidence": confidence,
        "molecule_model_basis": basis,
        "molecule_model_conflict_flag": conflict,
        "A_chrom": arm_a.chrom,
        "A_start": arm_a.start,
        "A_end": arm_a.end,
        "A_repeat_strand": arm_a.repeat_strand,
        "B_chrom": arm_b.chrom,
        "B_start": arm_b.start,
        "B_end": arm_b.end,
        "B_repeat_strand": arm_b.repeat_strand,
        "shared_transcript_support": bool(shared_ids),
        "shared_transcript_count": len(shared_ids),
        "shared_transcript_ids": _join_limited(shared_ids, max_ids),
        "A_transcript_count": len(a_by_id),
        "A_transcript_ids": _join_limited(a_by_id, max_ids),
        "B_transcript_count": len(b_by_id),
        "B_transcript_ids": _join_limited(b_by_id, max_ids),
        "shared_gtf_gene_ids": _join_limited(shared_gene_ids, max_ids),
        "shared_gtf_gene_names": _join_limited(shared_gene_names, max_ids),
        "same_gene_support": same_gene_support,
        "opposite_annotated_transcript_strands": opposite_tx_strands,
        f"{case_label}_both_strands": case_both_strands,
        f"{case_label}_arm_opposite": case_arm_opposite,
        f"{case_label}_arms_both_cov": case_arms_both_cov,
        "RI_exact_overlap_both_arms": ri_exact_both,
        "RI_overlap_both_arms": ri_overlap_both,
        "RI_gene_match_any": ri_gene_match,
        "intramolecular_support_count_0_to_5": intra_count,
        "intramolecular_support_features": ";".join(k for k, v in intra_flags.items() if v),
        "intermolecular_support_count_0_to_5": inter_count,
        "intermolecular_support_features": ";".join(k for k, v in inter_flags.items() if v),
        "A_annotation_category_model": a_cat,
        "B_annotation_category_model": b_cat,
        "model_is_prediction_not_ground_truth": True,
    }


def _merge_zrna(result: pd.DataFrame, zrna_summary: Path | None) -> pd.DataFrame:
    if zrna_summary is None:
        return result
    if not zrna_summary.is_file() or zrna_summary.stat().st_size == 0:
        raise FileNotFoundError(f"Missing/empty Z-RNA summary: {zrna_summary}")
    z = pd.read_csv(zrna_summary, low_memory=False)
    id_col = "pair_id" if "pair_id" in z.columns else "TE_pair" if "TE_pair" in z.columns else None
    if id_col is None:
        raise ValueError(f"Z-RNA summary lacks pair_id/TE_pair: {zrna_summary}")
    wanted = [
        id_col,
        "A_RNA_support_score",
        "A_RNA_support_class",
        "ZRNA_sequence_propensity_score",
        "ZRNA_propensity_score",
        "ZRNA_propensity_class",
        "ZRNA_priority_flag",
        "ZRNA_vs_A_interpretation",
    ]
    z = z[[c for c in wanted if c in z.columns]].copy().rename(columns={id_col: "pair_id"})
    return result.merge(z, on="pair_id", how="left", validate="one_to_one")


def run_molecule_model(args) -> None:
    output_dir = Path(args.output_dir)
    summary_path, summary = _read_summary(output_dir, args.analyze_subset, args.summary_in)
    gtf = Path(args.gtf)

    chroms = _collect_query_chromosomes(summary)
    tx_by_chrom = _read_transcript_spans(gtf, chroms)
    n_transcripts = sum(len(v) for v in tx_by_chrom.values())
    if n_transcripts == 0:
        raise ValueError(
            f"No transcript spans from {gtf} matched candidate chromosomes {sorted(chroms)}"
        )

    records = [
        _model_row(
            row,
            case_label=args.case_label,
            transcripts_by_chrom=tx_by_chrom,
            overlap_fraction=float(args.transcript_overlap_fraction),
            slop=int(args.transcript_containment_slop),
            max_ids=int(args.max_transcript_ids),
        )
        for _, row in summary.iterrows()
    ]
    model = pd.DataFrame(records)

    carry = [
        c
        for c in [
            "pair_id",
            "priority_rank",
            "priority_tier",
            "priority_gate_pass",
            "rank_score",
            "structure_priority_score",
            "A_SYMBOL",
            "B_SYMBOL",
            "A_TE_id",
            "B_TE_id",
            "genomic_orientation",
            "transcript_orientation",
            "A_annotation_category",
            "B_annotation_category",
        ]
        if c in summary.columns
    ]
    result = summary[carry].merge(model, on="pair_id", how="left", validate="one_to_one")

    zrna_path: Path | None = None
    if args.zrna_summary:
        zrna_path = Path(args.zrna_summary)
    elif args.include_default_zrna:
        candidate = (
            output_dir
            / args.analyze_subset
            / "summary"
            / f"TEpair_dsRNA_ZRNA_summary.{args.case_label}.csv"
        )
        if candidate.is_file() and candidate.stat().st_size > 0:
            zrna_path = candidate
    result = _merge_zrna(result, zrna_path)

    outdir = output_dir / args.analyze_subset / "summary"
    outdir.mkdir(parents=True, exist_ok=True)
    out_csv = outdir / f"TEpair_dsRNA_molecule_model.{args.case_label}.csv"
    out_counts = outdir / f"TEpair_dsRNA_molecule_model_counts.{args.case_label}.csv"
    out_meta = outdir / f"TEpair_dsRNA_molecule_model_metadata.{args.case_label}.txt"

    model_order = {
        "putative_intramolecular_foldback": 0,
        "putative_intermolecular_sense_antisense": 1,
        "ambiguous_inverted_TE_pair": 2,
        "direct_orientation_control": 3,
    }
    conf_order = {"high": 0, "moderate": 1, "low": 2, "unresolved": 3, "not_applicable": 4}
    result["_model_order"] = result["molecule_model"].map(model_order).fillna(9)
    result["_conf_order"] = result["molecule_model_confidence"].map(conf_order).fillna(9)
    result = result.sort_values(
        ["_model_order", "_conf_order", "intramolecular_support_count_0_to_5", "intermolecular_support_count_0_to_5"],
        ascending=[True, True, False, False],
        kind="mergesort",
    ).drop(columns=["_model_order", "_conf_order"])
    result.to_csv(out_csv, index=False)

    counts = (
        result.groupby(["molecule_model", "molecule_model_confidence"], dropna=False)
        .size()
        .rename("candidate_count")
        .reset_index()
        .sort_values(["candidate_count", "molecule_model"], ascending=[False, True])
    )
    counts.to_csv(out_counts, index=False)

    out_meta.write_text(
        "dsRNASeeker molecule-origin confidence annotation\n"
        f"summary_in={summary_path}\n"
        f"gtf={gtf}\n"
        f"case_label={args.case_label}\n"
        f"control_label={args.control_label}\n"
        f"analyze_subset={args.analyze_subset}\n"
        f"transcript_overlap_fraction={args.transcript_overlap_fraction}\n"
        f"transcript_containment_slop={args.transcript_containment_slop}\n"
        f"candidate_count={len(result)}\n"
        f"transcript_span_count={n_transcripts}\n"
        f"zrna_summary={zrna_path if zrna_path else 'not_used'}\n"
        "interpretation=Predicted molecule-origin compatibility, not direct proof of RNA topology.\n"
        "zrna_role=Z-RNA/A-form columns are carried for interpretation only and do not determine molecule origin.\n"
    )

    print(f"[OK] candidates={len(result)}")
    print(counts.to_string(index=False))
    print(f"[OK] {out_csv}")
    print(f"[OK] {out_counts}")
    print(f"[OK] {out_meta}")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Annotate dsRNASeeker candidates with conservative intramolecular/intermolecular compatibility."
    )
    p.add_argument("--output-dir", required=True)
    p.add_argument("--case-label", required=True)
    p.add_argument("--control-label", required=True)
    p.add_argument("--gtf", required=True)
    p.add_argument("--analyze-subset", default="inverted", choices=["inverted", "hairpin", "allpairs"])
    p.add_argument("--summary-in", default=None, help="Optional explicit TEpair_dsRNA_master.summary.with_RI.csv")
    p.add_argument("--zrna-summary", default=None, help="Optional explicit TEpair_dsRNA_ZRNA_summary.<case>.csv")
    p.add_argument(
        "--include-default-zrna",
        dest="include_default_zrna",
        action="store_true",
        default=True,
        help="Automatically merge the standard Z-RNA summary when present (default).",
    )
    p.add_argument(
        "--no-include-default-zrna",
        dest="include_default_zrna",
        action="store_false",
        help="Do not search for or merge the standard Z-RNA summary.",
    )
    p.add_argument(
        "--transcript-overlap-fraction",
        type=float,
        default=0.80,
        help="Minimum fraction of each TE arm that must lie within a transcript span (default 0.80).",
    )
    p.add_argument(
        "--transcript-containment-slop",
        type=int,
        default=0,
        help="Optional bp padding around transcript spans (default 0).",
    )
    p.add_argument(
        "--max-transcript-ids",
        type=int,
        default=20,
        help="Maximum transcript/gene identifiers retained per output cell (default 20).",
    )
    return p


def main() -> None:
    args = build_parser().parse_args()
    if not 0 < args.transcript_overlap_fraction <= 1:
        raise ValueError("--transcript-overlap-fraction must be in (0, 1]")
    if args.transcript_containment_slop < 0:
        raise ValueError("--transcript-containment-slop must be >= 0")
    if args.max_transcript_ids < 1:
        raise ValueError("--max-transcript-ids must be >= 1")
    run_molecule_model(args)


if __name__ == "__main__":
    main()
