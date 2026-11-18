#!/usr/bin/env python3
"""Prepare a gene-symbol count matrix and sample metadata from an nf-core/rnaseq
Salmon merged gene count matrix.

The preferred input is salmon.merged.gene_counts_length_scaled.tsv. These are
bias-corrected gene-level count-like values produced by tximport. Values are
rounded before DESeq2, as required by DESeqDataSetFromMatrix.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import re
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--matrix", required=True)
    p.add_argument("--gtf", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--case-label", required=True)
    p.add_argument("--control-label", required=True)
    p.add_argument("--case-samples", default="",
                   help="Comma-separated identifiers. Each must match exactly one matrix column by exact name or unique substring.")
    p.add_argument("--control-samples", default="")
    p.add_argument("--case-regex", default="")
    p.add_argument("--control-regex", default="")
    return p.parse_args()


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open("rt")


def parse_attrs(text: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for key, value in re.findall(r'(\S+)\s+"([^"]*)"', text):
        out[key] = value
    return out


def build_gene_map(gtf_path: Path) -> Tuple[Dict[str, str], Dict[str, str]]:
    exact: Dict[str, str] = {}
    stripped: Dict[str, str] = {}
    with open_text(gtf_path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] not in {"gene", "exon", "transcript"}:
                continue
            attrs = parse_attrs(fields[8])
            gid = attrs.get("gene_id", "").strip()
            gname = attrs.get("gene_name", "").strip() or attrs.get("gene", "").strip()
            if not gid:
                continue
            if not gname:
                # Some knownGene-style GTFs use a symbol-like gene_id already.
                gname = gid
            exact.setdefault(gid, gname)
            stripped.setdefault(re.sub(r"\.\d+$", "", gid), gname)
    if not exact:
        raise RuntimeError(f"No gene_id mappings were parsed from {gtf_path}")
    return exact, stripped


def split_tokens(text: str) -> List[str]:
    return [x.strip() for x in text.split(",") if x.strip()]


def resolve_tokens(tokens: Sequence[str], columns: Sequence[str], label: str) -> List[str]:
    selected: List[str] = []
    for token in tokens:
        if token in columns:
            matches = [token]
        else:
            matches = [c for c in columns if token in c]
        if len(matches) != 1:
            raise RuntimeError(
                f"{label} token {token!r} matched {len(matches)} columns: {matches[:20]}"
            )
        if matches[0] not in selected:
            selected.append(matches[0])
    return selected


def resolve_regex(pattern: str, columns: Sequence[str], label: str) -> List[str]:
    if not pattern:
        return []
    rx = re.compile(pattern)
    matches = [c for c in columns if rx.search(c)]
    if not matches:
        raise RuntimeError(f"{label} regex {pattern!r} matched no columns")
    return matches


def main() -> None:
    args = parse_args()
    matrix_path = Path(args.matrix)
    gtf_path = Path(args.gtf)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if not matrix_path.is_file():
        raise FileNotFoundError(matrix_path)
    if not gtf_path.is_file():
        raise FileNotFoundError(gtf_path)

    header = pd.read_csv(matrix_path, sep="\t", nrows=0).columns.tolist()
    if len(header) < 2:
        raise RuntimeError(f"Count matrix has fewer than two columns: {matrix_path}")
    gene_col = header[0]
    sample_columns = header[1:]

    case_cols = resolve_tokens(split_tokens(args.case_samples), sample_columns, "case")
    control_cols = resolve_tokens(split_tokens(args.control_samples), sample_columns, "control")
    case_cols += [x for x in resolve_regex(args.case_regex, sample_columns, "case") if x not in case_cols]
    control_cols += [x for x in resolve_regex(args.control_regex, sample_columns, "control") if x not in control_cols]

    overlap = sorted(set(case_cols) & set(control_cols))
    if overlap:
        raise RuntimeError(f"Columns selected in both groups: {overlap}")
    if not case_cols or not control_cols:
        raise RuntimeError(f"Need at least one case and one control column; case={case_cols}, control={control_cols}")

    print(f"[INFO] Matrix: {matrix_path}")
    print(f"[INFO] Gene identifier column: {gene_col}")
    print(f"[INFO] Case columns ({len(case_cols)}): {case_cols}")
    print(f"[INFO] Control columns ({len(control_cols)}): {control_cols}")

    usecols = [gene_col] + control_cols + case_cols
    df = pd.read_csv(matrix_path, sep="\t", usecols=usecols)
    df[gene_col] = df[gene_col].astype(str)
    for col in control_cols + case_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0.0)
        if (df[col] < 0).any():
            raise RuntimeError(f"Negative values found in {col}")

    exact, stripped = build_gene_map(gtf_path)

    def map_symbol(gid: str) -> str:
        gid = str(gid).strip()
        if gid in exact:
            return exact[gid]
        base = re.sub(r"\.\d+$", "", gid)
        if base in stripped:
            return stripped[base]
        # Preserve symbol-like identifiers; otherwise retain the ID so the
        # downstream MSigDB overlap check can fail explicitly if mapping is poor.
        return gid

    df.insert(0, "gene_symbol", df[gene_col].map(map_symbol))
    df = df.drop(columns=[gene_col])
    df = df[df["gene_symbol"].notna() & (df["gene_symbol"].astype(str).str.strip() != "")]

    numeric_cols = control_cols + case_cols
    grouped = df.groupby("gene_symbol", as_index=False)[numeric_cols].sum()
    grouped[numeric_cols] = np.rint(grouped[numeric_cols]).astype(np.int64)

    counts_out = outdir / "gene_raw_counts.tsv"
    samples_out = outdir / "gene_samples.tsv"
    grouped.to_csv(counts_out, sep="\t", index=False)

    sample_rows = []
    for col in control_cols:
        sample_rows.append({"sample_id": col, "condition": args.control_label})
    for col in case_cols:
        sample_rows.append({"sample_id": col, "condition": args.case_label})
    pd.DataFrame(sample_rows).to_csv(samples_out, sep="\t", index=False)

    print(f"[OK] Genes written: {len(grouped)}")
    print(f"[OK] {counts_out}")
    print(f"[OK] {samples_out}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        raise
