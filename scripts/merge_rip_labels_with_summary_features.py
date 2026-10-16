#!/usr/bin/env python3
"""Merge ternary RIP labels with dsRNASeeker ADPS features for model training."""
from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd

FEATURES = [
    "orientation_adps",
    "annotation_adps",
    "case_expression_adps",
    "energy_adps",
    "interface_adps",
    "case_editing_adps",
    "RI_adps",
    "control_RI_fraction",
]


def read_table(path_text: str) -> pd.DataFrame:
    path = Path(path_text)
    return pd.read_csv(
        path,
        sep="\t" if path.suffix.lower() in {".tsv", ".txt"} else ",",
        low_memory=False,
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", required=True)
    parser.add_argument("--labels", required=True)
    parser.add_argument("--dataset-id", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    summary = read_table(args.summary)
    labels = read_table(args.labels)
    required_summary = {"pair_id"}.union(FEATURES)
    missing = sorted(required_summary.difference(summary.columns))
    if missing:
        raise ValueError(f"Summary lacks required features: {missing}")
    if not {"pair_id", "label"}.issubset(labels.columns):
        raise ValueError("Labels require pair_id and label")

    # Authoritative source policy: labels contribute only labels/provenance,
    # numeric evidence is always taken from the explicitly selected summary.
    if labels["pair_id"].isna().any() or labels["pair_id"].duplicated().any():
        raise ValueError("Label table contains missing or duplicate pair_id")
    if summary["pair_id"].isna().any() or summary["pair_id"].duplicated().any():
        raise ValueError("Summary contains missing or duplicate pair_id")
    original = pd.to_numeric(labels["label"], errors="coerce")
    if ((labels["label"].notna() & original.isna()).any()
            or not original.dropna().isin((0, 1)).all()):
        raise ValueError("Labels must be 0, 1 or missing; fractional/non-numeric labels forbidden")
    # Drop stale feature copies (whether scored earlier or not).
    labels = labels.drop(columns=[c for c in [*FEATURES, "adaptive_priority_score"]
                                  if c in labels.columns]).copy()

    keep_summary = [
        "pair_id", "A_SYMBOL", "B_SYMBOL", "A_TE_id", "B_TE_id",
        *FEATURES,
        "adaptive_priority_score",
    ]
    keep_summary = [c for c in keep_summary if c in summary.columns]
    merged = labels.merge(
        summary[keep_summary],
        on="pair_id",
        how="left",
        validate="one_to_one",
        suffixes=("", "_summary"),
    )
    merged["dataset_id"] = args.dataset_id
    if "group_id" not in merged.columns:
        a = merged.get("A_SYMBOL", pd.Series("", index=merged.index)).fillna("").astype(str).str.upper()
        b = merged.get("B_SYMBOL", pd.Series("", index=merged.index)).fillna("").astype(str).str.upper()
        merged["group_id"] = [
            "|".join(sorted({x for x in (x1, x2) if x and x != "NAN"})) or pair
            for x1, x2, pair in zip(a, b, merged["pair_id"].astype(str))
        ]

    if merged[FEATURES].isna().any().any():
        raise ValueError("Missing one or more source summary ADPS features after pair-ID merge")
    if "adaptive_priority_score" not in merged.columns or merged["adaptive_priority_score"].isna().any():
        raise ValueError("Missing canonical adaptive_priority_score in source summary")
    missing_feature_rows = merged[FEATURES].isna().all(axis=1)
    if missing_feature_rows.any():
        examples = merged.loc[missing_feature_rows, "pair_id"].head().tolist()
        raise ValueError(f"Labels failed to match summary features; examples={examples}")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(out, sep="\t", index=False)
    print(f"[OK] rows={len(merged)} positives={(pd.to_numeric(merged.label, errors='coerce') == 1).sum()} "
          f"negatives={(pd.to_numeric(merged.label, errors='coerce') == 0).sum()}")
    print(f"[OK] {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
