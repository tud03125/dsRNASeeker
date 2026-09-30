#!/usr/bin/env python3
from __future__ import annotations

import argparse
import bisect
from pathlib import Path
import pandas as pd


def _numeric_int(series: pd.Series, label: str) -> pd.Series:
    x = pd.to_numeric(series, errors="raise")
    if x.isna().any():
        raise SystemExit(f"{label} contains missing numeric values")
    return x.astype("int64")


def read_candidate_arms(path: str | Path) -> pd.DataFrame:
    raw = pd.read_csv(
        path,
        sep="\t",
        header=None,
        comment="#",
        dtype=str,
        keep_default_na=False,
    )
    if raw.shape[1] < 5:
        raise SystemExit(
            "candidate-arm BED must have >=5 columns: chrom,start,end,pair_id,arm"
        )

    # Build a NEW DataFrame rather than assigning integers into columns that
    # pandas 3 may have created with StringDtype.
    out = pd.DataFrame({
        "chrom": raw.iloc[:, 0].astype(str).to_numpy(),
        "start": _numeric_int(raw.iloc[:, 1], "candidate start").to_numpy(),
        "end": _numeric_int(raw.iloc[:, 2], "candidate end").to_numpy(),
        "pair_id": raw.iloc[:, 3].astype(str).to_numpy(),
        "arm": raw.iloc[:, 4].astype(str).to_numpy(),
    })

    if (out["end"] < out["start"]).any():
        raise SystemExit("candidate-arm BED contains end < start")
    return out


def read_positive_bed(path: str | Path, id_col: int) -> pd.DataFrame:
    raw = pd.read_csv(
        path,
        sep="\t",
        header=None,
        comment="#",
        dtype=str,
        keep_default_na=False,
    )
    if raw.shape[1] < 3:
        raise SystemExit("positive BED must have >=3 columns")
    if id_col < 1:
        raise SystemExit("--positive-id-col is 1-based and must be >=1")

    if id_col <= raw.shape[1]:
        positive_id = raw.iloc[:, id_col - 1].astype(str)
    else:
        positive_id = pd.Series(
            [f"positive_{i+1:06d}" for i in range(len(raw))],
            dtype="object",
        )

    # Same pandas-3-safe construction: never mutate StringDtype columns in place.
    out = pd.DataFrame({
        "chrom": raw.iloc[:, 0].astype(str).to_numpy(),
        "start": _numeric_int(raw.iloc[:, 1], "positive start").to_numpy(),
        "end": _numeric_int(raw.iloc[:, 2], "positive end").to_numpy(),
        "positive_id": positive_id.astype(str).to_numpy(),
    })

    if (out["end"] < out["start"]).any():
        raise SystemExit("positive BED contains end < start")
    return out


def index_intervals(d: pd.DataFrame, idcol: str):
    out = {}
    for chrom, g in d.groupby("chrom", sort=False):
        g = g.sort_values(["start", "end"]).reset_index(drop=True)
        out[chrom] = (
            g["start"].astype("int64").tolist(),
            g["end"].astype("int64").tolist(),
            g[idcol].astype(str).tolist(),
        )
    return out


def overlaps(idx, chrom: str, start: int, end: int) -> set[str]:
    if chrom not in idx:
        return set()
    starts, ends, ids = idx[chrom]
    stop = bisect.bisect_left(starts, end)
    return {ids[i] for i in range(stop) if ends[i] > start}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--positive-bed", required=True)
    ap.add_argument("--positive-id-col", type=int, default=4)
    ap.add_argument(
        "--positive-id-list",
        help="Optional one-column list of assay-positive IDs. "
             "BED may contain all author clusters.",
    )
    ap.add_argument("--candidate-arms-bed", required=True)
    ap.add_argument("--dataset-id", required=True)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()

    pos = read_positive_bed(a.positive_bed, a.positive_id_col)

    if a.positive_id_list:
        wanted_df = pd.read_csv(
            a.positive_id_list,
            header=None,
            dtype=str,
            keep_default_na=False,
        )
        wanted = {x for x in wanted_df.iloc[:, 0].astype(str) if x != ""}

        pos_ids = set(pos["positive_id"].astype(str))
        found = wanted & pos_ids
        missing = sorted(wanted - pos_ids)

        pos = pos[pos["positive_id"].isin(wanted)].copy()

        print(
            f"[INFO] positive ID filter: {len(wanted)} requested; "
            f"{len(found)} found in BED; {len(missing)} missing"
        )

        if missing:
            od = Path(a.outdir)
            od.mkdir(parents=True, exist_ok=True)
            (od / "positive_ids_missing_from_bed.txt").write_text(
                "\n".join(missing) + "\n"
            )

        if not found:
            raise SystemExit("none of the requested assay-positive IDs were found in BED")

    arms = read_candidate_arms(a.candidate_arms_bed)

    bad = arms.groupby("pair_id")["chrom"].nunique()
    if (bad > 1).any():
        raise SystemExit("candidate pair crosses chromosomes")

    span = (
        arms.groupby(["pair_id", "chrom"], as_index=False)
        .agg(start=("start", "min"), end=("end", "max"))
    )

    armidx = index_intervals(arms, "pair_id")
    spanidx = index_intervals(span, "pair_id")

    rows = []
    # One record per unique author positive ID, even when a cluster has
    # multiple BED blocks.
    for pid, g in pos.groupby("positive_id", sort=False):
        arm_pairs: set[str] = set()
        span_pairs: set[str] = set()

        for r in g.itertuples(index=False):
            arm_pairs |= overlaps(armidx, r.chrom, int(r.start), int(r.end))
            span_pairs |= overlaps(spanidx, r.chrom, int(r.start), int(r.end))

        rows.append({
            "positive_id": str(pid),
            "n_blocks": int(len(g)),
            "chromosomes": ";".join(sorted(set(g["chrom"].astype(str)))),
            "captured_by_candidate_arm": bool(arm_pairs),
            "captured_by_pair_span": bool(span_pairs),
            "n_arm_overlapping_pairs": int(len(arm_pairs)),
            "n_span_overlapping_pairs": int(len(span_pairs)),
            "arm_pair_ids": ";".join(sorted(arm_pairs)),
            "span_pair_ids": ";".join(sorted(span_pairs)),
        })

    out = pd.DataFrame(rows)
    od = Path(a.outdir)
    od.mkdir(parents=True, exist_ok=True)

    bypos = od / "candidate_generation_capture_by_positive_cluster.tsv"
    summary_path = od / "candidate_generation_capture_summary.tsv"
    missed_path = od / "missed_positive_clusters.arm.tsv"

    out.to_csv(bypos, sep="\t", index=False)

    n = len(out)
    arm = int(out["captured_by_candidate_arm"].sum()) if n else 0
    span_n = int(out["captured_by_pair_span"].sum()) if n else 0

    summary = pd.DataFrame([{
        "dataset_id": a.dataset_id,
        "n_positive_clusters": n,
        "n_captured_by_candidate_arm": arm,
        "candidate_arm_capture_fraction": arm / n if n else float("nan"),
        "n_captured_by_pair_span": span_n,
        "pair_span_capture_fraction": span_n / n if n else float("nan"),
        "n_candidate_pairs": int(arms["pair_id"].nunique()),
        "positive_unit": (
            "unique author cluster ID; split BED blocks grouped before capture"
        ),
    }])

    summary.to_csv(summary_path, sep="\t", index=False)
    out.loc[~out["captured_by_candidate_arm"]].to_csv(
        missed_path, sep="\t", index=False
    )

    print(summary.to_string(index=False))
    print(f"[OK] {summary_path}")
    print(f"[OK] {bypos}")
    print(f"[OK] {missed_path}")


if __name__ == "__main__":
    main()
