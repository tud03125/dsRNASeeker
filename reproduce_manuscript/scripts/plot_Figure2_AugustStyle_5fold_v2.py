#!/usr/bin/env python3
"""Rebuild manuscript Figure 2 using the August five-fold presentation.

Scientific design
-----------------
* Three principal benchmark panels:
    GSE308488/9 FLAG-ZBP1, 1 kb
    GSE184962/4 Z22, 5 kb expanded quantitative benchmark
    GSE244103 J2, 1 kb
* All methods are evaluated on the SAME five StratifiedGroupKFold outer folds.
* ADPS and external-method scores are fixed (not trained) and are merely evaluated
  within each held-out fold.
* The supervised curve uses current fixed-L2 same-study OOF predictions.
* Main dsRNAscan comparator is the prespecified binary two-arm projection.
* Mean PR/ROC curves are obtained by interpolating each fold curve to a common grid
  and averaging across folds. Legend values are mean ± sample SD of fold-specific
  AP or ROC-AUC, reproducing the August manuscript presentation style.

This script intentionally does NOT choose a window based on observed performance.
The Z22 5-kb row is the prespecified adequately-powered quantitative analysis
(>=10 positives and >=10 confident negatives); the 1-kb Z22 result remains a
separate sensitivity analysis.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import (
    average_precision_score,
    precision_recall_curve,
    roc_auc_score,
    roc_curve,
)
from sklearn.model_selection import StratifiedGroupKFold


DATASETS = [
    {
        "tag": "GSE308489__w1000_standardized",
        "dataset_id": "GSE308489",
        "variant": "w1000_standardized",
        "title": "GSE308488/9 HSV-1–FLAG-ZBP1",
        "subtitle": "1-kb common-radius",
    },
    {
        "tag": "GSE184962_4__w5000_expanded",
        "dataset_id": "GSE184962_4",
        "variant": "w5000_expanded",
        "title": "GSE184962/4 IFNβ–Z22",
        "subtitle": "5-kb expanded quantitative benchmark",
    },
    {
        "tag": "GSE244103__w1000_standardized",
        "dataset_id": "GSE244103",
        "variant": "w1000_standardized",
        "title": "GSE244103 DHX9–J2",
        "subtitle": "1-kb common-radius",
    },
]

METHODS = [
    ("ADPS", "score_integrated_adps", "dsRNASeeker integrated ADPS", "-"),
    ("Supervised", "supervised_score", "dsRNASeeker supervised re-ranking", "--"),
    ("EER", "eer_delta_arm_count", "EER projected score", "-."),
    ("dsRNAscan", "dsrnascan_two_arm_match", "dsRNAscan two-arm match", ":"),
]

NATIVE_METHODS = [
    ("dsRNAscan stability", "dsrnascan_stability_max", "dsRNAscan stability", (0, (5, 2))),
    ("dsRNAscan probing", "dsrnascan_probing_max", "dsRNAscan probing", (0, (1, 1))),
]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--root",
        default="reproduce_manuscript/input_root",
        help="Portable reproduction input root (default: reproduce_manuscript/input_root).",
    )
    p.add_argument("--supervised-dir", default=None,
                   help="Override supervised benchmark directory.")
    p.add_argument("--audit-dir", default=None,
                   help="Override supervised input/audit directory.")
    p.add_argument("--benchmark-table", default=None,
                   help="Override all_method_candidate_scores.tsv.gz path.")
    p.add_argument("--outdir", default=None)
    p.add_argument("--grid-points", type=int, default=201)
    p.add_argument("--include-native-dsrnascan", action="store_true")
    p.add_argument("--dpi", type=int, default=400)
    return p.parse_args()


def pick_dir(root: Path, preferred: str, fallback: str) -> Path:
    p = root / preferred
    if p.exists():
        return p
    q = root / fallback
    if q.exists():
        return q
    raise FileNotFoundError(f"Neither {p} nor {q} exists")


def read_tsv_or_csv(path: Path) -> pd.DataFrame:
    if str(path).endswith(".tsv") or str(path).endswith(".tsv.gz"):
        return pd.read_csv(path, sep="\t", low_memory=False)
    return pd.read_csv(path, low_memory=False)


def make_fold_assignments(audit: pd.DataFrame, *, n_splits: int, seed: int) -> pd.DataFrame:
    d = audit.loc[audit["label"].notna(), ["pair_id", "label", "group_id"]].copy()
    d["pair_id"] = d["pair_id"].astype(str)
    d["label"] = d["label"].astype(int)
    d["group_id"] = d["group_id"].astype(str)
    d["outer_fold"] = -1

    cv = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    for fold, (_, test_idx) in enumerate(
        cv.split(d, d["label"].to_numpy(), groups=d["group_id"].to_numpy()), start=1
    ):
        d.iloc[test_idx, d.columns.get_loc("outer_fold")] = fold

    if (d["outer_fold"] < 1).any():
        raise RuntimeError("Failed to assign all labeled candidates to an outer fold")
    return d


def merge_current_scores(
    fold_df: pd.DataFrame,
    pred: pd.DataFrame,
    ext: pd.DataFrame,
    *,
    tag: str,
    dataset_id: str,
    variant: str,
) -> pd.DataFrame:
    # Fixed L2 is the prespecified same-study supervised model used for Figure 2.
    p = pred.loc[(pred["dataset"].astype(str) == tag) &
                 (pred["model"].astype(str) == "fixed_l2_raw7")].copy()
    if p.empty:
        raise ValueError(f"No fixed_l2_raw7 OOF predictions for {tag}")
    p["pair_id"] = p["pair_id"].astype(str)

    required_pred = ["pair_id", "score_integrated_adps", "supervised_score"]
    missing = [c for c in required_pred if c not in p.columns]
    if missing:
        raise ValueError(f"{tag} supervised predictions missing columns: {missing}")

    e = ext.loc[(ext["dataset_id"].astype(str) == dataset_id) &
                (ext["analysis_variant"].astype(str) == variant)].copy()
    if e.empty:
        raise ValueError(f"No external-method candidate scores for {dataset_id}/{variant}")
    e["pair_id"] = e["pair_id"].astype(str)

    ext_cols = [
        "pair_id", "eer_delta_arm_count", "dsrnascan_two_arm_match",
        "dsrnascan_stability_max", "dsrnascan_probing_max",
    ]
    ext_cols = [c for c in ext_cols if c in e.columns]

    m = fold_df.merge(p[required_pred], on="pair_id", how="left", validate="one_to_one")
    m = m.merge(e[ext_cols], on="pair_id", how="left", validate="one_to_one")

    # Every labeled candidate must have the ADPS + supervised scores.
    for c in ["score_integrated_adps", "supervised_score"]:
        if m[c].isna().any():
            bad = m.loc[m[c].isna(), "pair_id"].head(10).tolist()
            raise ValueError(f"{tag}: missing {c} for labeled candidates, e.g. {bad}")

    return m


def fold_curves_and_metrics(d: pd.DataFrame, score_col: str, recall_grid, fpr_grid):
    pr_curves = []
    roc_curves = []
    rows = []

    for fold in sorted(d["outer_fold"].unique()):
        x = d.loc[d["outer_fold"] == fold, ["label", score_col]].dropna()
        y = x["label"].astype(int).to_numpy()
        s = pd.to_numeric(x[score_col], errors="coerce").to_numpy()

        if len(np.unique(y)) < 2:
            raise ValueError(f"Fold {fold} for {score_col} does not contain both label classes")

        ap = average_precision_score(y, s)
        auc = roc_auc_score(y, s)
        rows.append({
            "outer_fold": int(fold),
            "n": int(len(y)),
            "n_positive": int(y.sum()),
            "n_negative": int((1-y).sum()),
            "average_precision": float(ap),
            "roc_auc": float(auc),
        })

        precision, recall, _ = precision_recall_curve(y, s)
        # precision_recall_curve returns decreasing recall. Reverse and collapse
        # duplicate recall coordinates by retaining the maximum precision.
        tmp = pd.DataFrame({"recall": recall[::-1], "precision": precision[::-1]})
        tmp = tmp.groupby("recall", as_index=False)["precision"].max().sort_values("recall")
        interp_p = np.interp(recall_grid, tmp["recall"], tmp["precision"])
        pr_curves.append(interp_p)

        fpr, tpr, _ = roc_curve(y, s)
        interp_t = np.interp(fpr_grid, fpr, tpr)
        interp_t[0] = 0.0
        interp_t[-1] = 1.0
        roc_curves.append(interp_t)

    met = pd.DataFrame(rows)
    return (
        np.mean(np.vstack(pr_curves), axis=0),
        np.mean(np.vstack(roc_curves), axis=0),
        met,
    )


def metric_summary(met: pd.DataFrame) -> dict:
    return {
        "ap_mean": float(met["average_precision"].mean()),
        "ap_sd": float(met["average_precision"].std(ddof=1)),
        "auc_mean": float(met["roc_auc"].mean()),
        "auc_sd": float(met["roc_auc"].std(ddof=1)),
    }


def main():
    args = parse_args()
    root = Path(args.root)
    supervised_dir = Path(args.supervised_dir) if args.supervised_dir else pick_dir(
        root, "supervised_benchmark_expanded", "supervised_benchmark"
    )
    audit_dir = Path(args.audit_dir) if args.audit_dir else pick_dir(
        root, "supervised_input_expanded", "supervised_input"
    )
    benchmark_table = Path(args.benchmark_table) if args.benchmark_table else (
        root / "benchmark" / "all_method_candidate_scores.tsv.gz"
    )
    outdir = Path(args.outdir) if args.outdir else (
        Path("reproduce_manuscript/reproduced/Figure2")
    )
    outdir.mkdir(parents=True, exist_ok=True)

    pred_path = supervised_dir / "same_study_oof_predictions.csv"
    meta_path = supervised_dir / "run_metadata.json"

    for p in [pred_path, meta_path, benchmark_table]:
        if not p.exists():
            raise FileNotFoundError(p)

    pred = pd.read_csv(pred_path, low_memory=False)
    ext = pd.read_csv(benchmark_table, sep="\t", low_memory=False)
    meta = json.load(open(meta_path))

    n_splits = int(meta.get("outer_folds", 5))
    seed = int(meta.get("seed", 20260925))
    if n_splits != 5:
        raise ValueError(f"August-style Figure 2 requires 5 outer folds; metadata says {n_splits}")

    recall_grid = np.linspace(0.0, 1.0, args.grid_points)
    fpr_grid = np.linspace(0.0, 1.0, args.grid_points)

    methods = list(METHODS)
    if args.include_native_dsrnascan:
        methods += NATIVE_METHODS

    fig, axes = plt.subplots(3, 2, figsize=(12.6, 14.2), constrained_layout=True)
    fold_rows = []
    summary_rows = []
    provenance_rows = []

    for row_i, spec in enumerate(DATASETS):
        tag = spec["tag"]
        audit_path = audit_dir / f"{tag}.audit.tsv"
        if not audit_path.exists():
            raise FileNotFoundError(audit_path)
        audit = pd.read_csv(audit_path, sep="\t", low_memory=False)

        fold_df = make_fold_assignments(audit, n_splits=n_splits, seed=seed)
        d = merge_current_scores(
            fold_df, pred, ext,
            tag=tag, dataset_id=spec["dataset_id"], variant=spec["variant"]
        )

        n_pos = int(d["label"].sum())
        n_neg = int((1-d["label"]).sum())
        prevalence = n_pos / len(d)

        ax_pr, ax_roc = axes[row_i, 0], axes[row_i, 1]

        for method_key, score_col, display, ls in methods:
            if score_col not in d.columns or d[score_col].notna().sum() == 0:
                continue
            mean_pr, mean_roc, met = fold_curves_and_metrics(d, score_col, recall_grid, fpr_grid)
            sm = metric_summary(met)

            ax_pr.plot(
                recall_grid, mean_pr, linestyle=ls, linewidth=2,
                label=f"{display} (AP={sm['ap_mean']:.3f} ± {sm['ap_sd']:.3f})"
            )
            ax_roc.plot(
                fpr_grid, mean_roc, linestyle=ls, linewidth=2,
                label=f"{display} (AUC={sm['auc_mean']:.3f} ± {sm['auc_sd']:.3f})"
            )

            for _, rr in met.iterrows():
                fold_rows.append({
                    "dataset": tag,
                    "display_dataset": spec["title"],
                    "analysis_variant": spec["variant"],
                    "method": method_key,
                    "score_column": score_col,
                    **rr.to_dict(),
                })
            summary_rows.append({
                "dataset": tag,
                "display_dataset": spec["title"],
                "analysis_variant": spec["variant"],
                "method": method_key,
                "score_column": score_col,
                "n_labeled": len(d),
                "n_positive": n_pos,
                "n_negative": n_neg,
                "positive_prevalence": prevalence,
                **sm,
            })

        ax_pr.axhline(prevalence, linestyle=":", linewidth=1.2, alpha=0.8,
                      label=f"Positive prevalence={prevalence:.3f}")
        ax_roc.plot([0, 1], [0, 1], linestyle=":", linewidth=1.2, alpha=0.8,
                    label="Chance")

        title = f"{spec['title']} (+{n_pos}/−{n_neg})\n{spec['subtitle']}"
        ax_pr.set_title(title + " — precision–recall", fontsize=11)
        ax_roc.set_title(title + " — ROC", fontsize=11)
        ax_pr.set_xlabel("Recall")
        ax_pr.set_ylabel("Precision")
        ax_roc.set_xlabel("False-positive rate")
        ax_roc.set_ylabel("True-positive rate")
        ax_pr.set_xlim(0, 1); ax_pr.set_ylim(0, 1.03)
        ax_roc.set_xlim(0, 1); ax_roc.set_ylim(0, 1.03)
        ax_pr.grid(alpha=0.16)
        ax_roc.grid(alpha=0.16)
        ax_pr.legend(fontsize=7.4, loc="best", frameon=True)
        ax_roc.legend(fontsize=7.4, loc="best", frameon=True)

        provenance_rows.append({
            "dataset": tag,
            "audit_table": str(audit_path),
            "supervised_predictions": str(pred_path),
            "external_scores": str(benchmark_table),
            "outer_folds": n_splits,
            "fold_seed": seed,
            "fold_algorithm": "StratifiedGroupKFold(shuffle=True)",
            "n_labeled": len(d),
            "n_positive": n_pos,
            "n_negative": n_neg,
        })

    fig.suptitle(
        "Benchmark performance of dsRNASeeker and external methods across three selected studies\n"
        "Mean five-fold grouped PR/ROC curves; legends show fold-specific metric mean ± SD",
        fontsize=14,
    )

    png = outdir / "Figure2_AugustStyle_updated_5fold.png"
    pdf = outdir / "Figure2_AugustStyle_updated_5fold.pdf"
    svg = outdir / "Figure2_AugustStyle_updated_5fold.svg"
    fig.savefig(png, dpi=args.dpi, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(svg, bbox_inches="tight")
    plt.close(fig)

    pd.DataFrame(fold_rows).to_csv(outdir / "Figure2_fold_metrics.csv", index=False)
    pd.DataFrame(summary_rows).to_csv(outdir / "Figure2_metric_summary_meanSD.csv", index=False)
    pd.DataFrame(provenance_rows).to_csv(outdir / "Figure2_input_provenance.csv", index=False)

    # Save exactly which labeled pair went to which fold for auditability.
    fold_assignments = []
    for spec in DATASETS:
        audit = pd.read_csv(audit_dir / f"{spec['tag']}.audit.tsv", sep="\t", low_memory=False)
        f = make_fold_assignments(audit, n_splits=n_splits, seed=seed)
        f.insert(0, "dataset", spec["tag"])
        fold_assignments.append(f)
    pd.concat(fold_assignments, ignore_index=True).to_csv(
        outdir / "Figure2_outer_fold_assignments.tsv", sep="\t", index=False
    )

    print(f"[OK] {png}")
    print(f"[OK] {pdf}")
    print(f"[OK] {svg}")
    print(f"[OK] {outdir / 'Figure2_fold_metrics.csv'}")
    print(f"[OK] {outdir / 'Figure2_metric_summary_meanSD.csv'}")
    print(f"[OK] fold seed={seed}; n_splits={n_splits}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
