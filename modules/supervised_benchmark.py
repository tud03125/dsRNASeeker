from __future__ import annotations

"""Manifest-driven supervised evaluation for dsRNASeeker.

This module generalizes the manuscript Stage 05 analysis. No accessions or
filesystem paths are embedded in Python: users supply a TSV/CSV manifest.
"""

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

from modules.supervised import (
    RAW7_FEATURES,
    ROUTE6_FEATURES,
    ROUTE_DEFINITIONS,
    _class_balanced_weights,
    _family_class_balanced_weights,
    _fit_checked,
    _fit_selected_with_fallback,
    _model,
    _safe_metric,
    _select_configuration,
    _choose_cv_splits,
    ensure_route_features,
)

COMPONENT_ALIASES = {
    "orientation_adps": ("orientation_adps", "orientation_adps_z01", "component_orientation"),
    "annotation_adps": ("annotation_adps", "annotation_adps_z01", "component_annotation"),
    "case_expression_adps": ("case_expression_adps", "case_expression_adps_z01", "component_expression"),
    "energy_adps": ("energy_adps", "energy_adps_z01", "component_energy"),
    "interface_adps": ("interface_adps", "interface_adps_z01", "component_interface"),
    "case_editing_adps": ("case_editing_adps", "case_editing_adps_z01", "component_editing"),
    "RI_adps": ("RI_adps", "RI_adps_z01", "component_RI"),
}
BASELINE_ALIASES = {
    "score_integrated_adps": (
        "score_integrated_adps", "adaptive_priority_score", "case_priority_score",
        "rank_score", "evidence_score_raw",
    ),
    "score_equal_weight_all": ("score_equal_weight_all",),
    "score_duplex_structure_only": ("score_duplex_structure_only",),
    "score_route_any_max": ("score_route_any_max",),
}
REQUIRED_MANIFEST = {
    "dataset_id", "analysis_variant", "study_family", "audit_table",
    "same_study_eligible", "loso_representative", "loso_target",
}


@dataclass(frozen=True)
class VariantSpec:
    dataset_id: str
    analysis_variant: str
    study_family: str
    audit_table: Path
    same_study_eligible: bool
    loso_representative: bool
    loso_target: bool


def _read_table(path: Path) -> pd.DataFrame:
    if path.suffix.lower() in {".tsv", ".txt"}:
        return pd.read_csv(path, sep="\t", low_memory=False)
    return pd.read_csv(path, low_memory=False)


def _parse_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    if text in {"1", "true", "t", "yes", "y"}:
        return True
    if text in {"0", "false", "f", "no", "n", "", "nan", "none"}:
        return False
    raise ValueError(f"Cannot parse boolean value: {value!r}")


def _first_existing(columns: Iterable[str], aliases: Sequence[str]) -> str | None:
    present = set(columns)
    return next((name for name in aliases if name in present), None)


def read_manifest(path: str | Path, input_root: str | Path | None = None) -> list[VariantSpec]:
    path = Path(path)
    manifest = _read_table(path)
    missing = sorted(REQUIRED_MANIFEST.difference(manifest.columns))
    if missing:
        raise ValueError(f"Manifest {path} missing columns: {missing}")
    root = Path(input_root) if input_root is not None else path.parent
    specs: list[VariantSpec] = []
    for _, row in manifest.iterrows():
        audit = Path(str(row["audit_table"]))
        if not audit.is_absolute():
            audit = root / audit
        specs.append(VariantSpec(
            dataset_id=str(row["dataset_id"]),
            analysis_variant=str(row["analysis_variant"]),
            study_family=str(row["study_family"]),
            audit_table=audit,
            same_study_eligible=_parse_bool(row["same_study_eligible"]),
            loso_representative=_parse_bool(row["loso_representative"]),
            loso_target=_parse_bool(row["loso_target"]),
        ))
    keys = [(x.dataset_id, x.analysis_variant) for x in specs]
    if len(keys) != len(set(keys)):
        raise ValueError("Manifest contains duplicate dataset_id/analysis_variant rows")
    families = {x.study_family for x in specs if x.loso_target}
    for family in families:
        reps = [x for x in specs if x.study_family == family and x.loso_representative]
        if len(reps) != 1:
            raise ValueError(
                f"Study family {family!r} must have exactly one loso_representative; found {len(reps)}"
            )
    return specs


def normalize_audit_table(spec: VariantSpec) -> pd.DataFrame:
    if not spec.audit_table.is_file() or spec.audit_table.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty audit table: {spec.audit_table}")
    raw = _read_table(spec.audit_table)
    for col in ("pair_id", "label", "group_id"):
        if col not in raw.columns:
            raise ValueError(f"{spec.audit_table} missing required column {col!r}")
    out = raw[["pair_id", "label", "group_id"]].copy()
    out["label"] = pd.to_numeric(out["label"], errors="coerce")
    out = out.loc[out["label"].notna()].copy()
    out["label"] = out["label"].astype(int)
    if not set(out["label"].unique()).issubset({0, 1}):
        raise ValueError(f"Unexpected labels in {spec.audit_table}: {sorted(out['label'].unique())}")
    out["group_id"] = out["group_id"].astype(str)
    for canonical, aliases in COMPONENT_ALIASES.items():
        source = _first_existing(raw.columns, aliases)
        out[canonical] = (
            pd.to_numeric(raw.loc[out.index, source], errors="coerce")
            if source else np.nan
        )
    out = ensure_route_features(out)
    source = _first_existing(raw.columns, BASELINE_ALIASES["score_integrated_adps"])
    out["score_integrated_adps"] = (
        pd.to_numeric(raw.loc[out.index, source], errors="coerce")
        if source else out[RAW7_FEATURES].mean(axis=1, skipna=True)
    )
    source = _first_existing(raw.columns, BASELINE_ALIASES["score_equal_weight_all"])
    out["score_equal_weight_all"] = (
        pd.to_numeric(raw.loc[out.index, source], errors="coerce")
        if source else out[RAW7_FEATURES].mean(axis=1, skipna=True)
    )
    source = _first_existing(raw.columns, BASELINE_ALIASES["score_duplex_structure_only"])
    out["score_duplex_structure_only"] = (
        pd.to_numeric(raw.loc[out.index, source], errors="coerce")
        if source else out[["energy_adps", "interface_adps"]].mean(axis=1, skipna=True)
    )
    source = _first_existing(raw.columns, BASELINE_ALIASES["score_route_any_max"])
    out["score_route_any_max"] = (
        pd.to_numeric(raw.loc[out.index, source], errors="coerce")
        if source else out[ROUTE6_FEATURES].max(axis=1, skipna=True)
    )
    out["dataset_id"] = spec.dataset_id
    out["analysis_variant"] = spec.analysis_variant
    out["study_family"] = spec.study_family
    out["prefixed_group_id"] = (
        spec.dataset_id + "::" + spec.analysis_variant + "::" + out["group_id"].astype(str)
    )
    return out.reset_index(drop=True)


def safe_metrics(y: Sequence[int], score: Sequence[float]) -> dict[str, float | int]:
    y = pd.Series(y).astype(int)
    s = pd.to_numeric(pd.Series(score), errors="coerce")
    keep = s.notna()
    y, s = y.loc[keep], s.loc[keep]
    result: dict[str, float | int] = {
        "n": int(len(y)), "n_positive": int(y.eq(1).sum()),
        "n_negative": int(y.eq(0).sum()),
        "prevalence": float(y.mean()) if len(y) else np.nan,
        "average_precision": np.nan, "roc_auc": np.nan,
    }
    if y.nunique() == 2:
        result["average_precision"] = float(average_precision_score(y, s))
        result["roc_auc"] = float(roc_auc_score(y, s))
    return result


def group_bootstrap(
    frame: pd.DataFrame,
    score_col: str,
    *,
    reference_col: str | None,
    n_boot: int,
    seed: int,
) -> dict[str, float | int]:
    columns = ["label", "group_id", score_col]
    if reference_col and reference_col != score_col:
        columns.append(reference_col)
    labeled = frame[columns].copy()
    groups = labeled["group_id"].astype(str).unique()
    if int(n_boot) <= 0 or len(groups) == 0:
        return {
            "ap_ci_low": np.nan, "ap_ci_high": np.nan,
            "auc_ci_low": np.nan, "auc_ci_high": np.nan,
            "delta_ap_ci_low": np.nan, "delta_ap_ci_high": np.nan,
            "delta_auc_ci_low": np.nan, "delta_auc_ci_high": np.nan,
            "bootstrap_valid": 0,
        }
    rng = np.random.default_rng(seed)
    aps: list[float] = []; aucs: list[float] = []
    daps: list[float] = []; daucs: list[float] = []
    for _ in range(int(n_boot)):
        sampled = rng.choice(groups, size=len(groups), replace=True)
        boot = pd.concat([
            labeled.loc[labeled["group_id"].astype(str).eq(g)] for g in sampled
        ], ignore_index=True)
        metric = safe_metrics(boot["label"], boot[score_col])
        if not np.isfinite(metric["average_precision"]):
            continue
        aps.append(float(metric["average_precision"])); aucs.append(float(metric["roc_auc"]))
        if reference_col:
            ref = safe_metrics(boot["label"], boot[reference_col])
            if np.isfinite(ref["average_precision"]):
                daps.append(float(metric["average_precision"]) - float(ref["average_precision"]))
                daucs.append(float(metric["roc_auc"]) - float(ref["roc_auc"]))
    def ci(values: list[float]) -> tuple[float, float]:
        if not values: return np.nan, np.nan
        return float(np.quantile(values, .025)), float(np.quantile(values, .975))
    ap_lo, ap_hi = ci(aps); auc_lo, auc_hi = ci(aucs)
    dap_lo, dap_hi = ci(daps); dauc_lo, dauc_hi = ci(daucs)
    return {
        "ap_ci_low": ap_lo, "ap_ci_high": ap_hi,
        "auc_ci_low": auc_lo, "auc_ci_high": auc_hi,
        "delta_ap_ci_low": dap_lo, "delta_ap_ci_high": dap_hi,
        "delta_auc_ci_low": dauc_lo, "delta_auc_ci_high": dauc_hi,
        "bootstrap_valid": int(len(aps)),
    }


def _selection_args(*, inner_folds: int, seed: int, metric: str, c_grid: str, l1_grid: str):
    return argparse.Namespace(
        supervised_feature_panel="auto_prespecified",
        supervised_tune=True,
        supervised_model="elasticnet",
        supervised_inner_cv_folds=int(inner_folds),
        supervised_random_state=int(seed),
        supervised_selection_metric=str(metric),
        supervised_c_grid=str(c_grid),
        supervised_l1_ratio_grid=str(l1_grid),
        supervised_c=1.0,
        supervised_l1_ratio=0.0,
    )


def _available_features(train: pd.DataFrame, target: pd.DataFrame, columns: Sequence[str]) -> list[str]:
    features: list[str] = []
    for col in columns:
        if col not in train.columns or col not in target.columns:
            continue
        x = pd.to_numeric(train[col], errors="coerce")
        if x.notna().sum() >= 2 and x.nunique(dropna=True) > 1:
            features.append(col)
    return features


def fixed_grouped_oof(data: pd.DataFrame, *, outer_folds: int, seed: int) -> tuple[np.ndarray, pd.DataFrame]:
    y = data["label"].astype(int).reset_index(drop=True)
    groups = data["group_id"].astype(str).reset_index(drop=True)
    splits, method = _choose_cv_splits(data, y, groups, cv_folds=outer_folds, random_state=seed)
    if not splits:
        raise ValueError("Could not construct grouped outer folds for fixed L2 model")
    features = _available_features(data, data, RAW7_FEATURES)
    if not features:
        raise ValueError("No nonconstant raw7 features available")
    oof = np.full(len(data), np.nan)
    rows: list[dict[str, Any]] = []
    for fold, (train_idx, test_idx) in enumerate(splits, start=1):
        train, test = data.iloc[train_idx], data.iloc[test_idx]
        y_train, y_test = train["label"].astype(int), test["label"].astype(int)
        model = _model(seed + fold, model_type="legacy_l2", C=1.0, l1_ratio=None)
        model, converged = _fit_checked(
            model, train[features].apply(pd.to_numeric, errors="coerce"),
            y_train, _class_balanced_weights(y_train),
        )
        if not converged:
            raise RuntimeError(f"Fixed L2 failed to converge in outer fold {fold}")
        score = model.predict_proba(test[features].apply(pd.to_numeric, errors="coerce"))[:, 1]
        oof[test_idx] = score
        rows.append({
            "outer_fold": fold, "status": "used", "cv_method": method,
            "n_train": len(train_idx), "n_test": len(test_idx),
            "n_test_positive": int(y_test.eq(1).sum()),
            "n_test_negative": int(y_test.eq(0).sum()),
            "feature_panel": "raw7", "features_used": ";".join(features),
            "C": 1.0, "l1_ratio": 0.0,
            "average_precision": _safe_metric(average_precision_score, y_test, score),
            "roc_auc": _safe_metric(roc_auc_score, y_test, score),
        })
    if not np.isfinite(oof).all():
        raise RuntimeError(f"Incomplete fixed-L2 OOF coverage: {int((~np.isfinite(oof)).sum())} rows")
    return oof, pd.DataFrame(rows)


def nested_grouped_oof(
    data: pd.DataFrame, *, outer_folds: int, inner_folds: int, seed: int,
    selection_metric: str, c_grid: str, l1_grid: str,
    output_dir: Path, prefix: str,
) -> tuple[np.ndarray, pd.DataFrame]:
    y = data["label"].astype(int).reset_index(drop=True)
    groups = data["group_id"].astype(str).reset_index(drop=True)
    splits, method = _choose_cv_splits(data, y, groups, cv_folds=outer_folds, random_state=seed)
    if not splits:
        raise ValueError("Could not construct grouped outer folds for nested model")
    oof = np.full(len(data), np.nan)
    rows: list[dict[str, Any]] = []
    for fold, (train_idx, test_idx) in enumerate(splits, start=1):
        train = data.iloc[train_idx].reset_index(drop=True)
        test = data.iloc[test_idx].reset_index(drop=True)
        y_train, y_test = train["label"].astype(int), test["label"].astype(int)
        args = _selection_args(
            inner_folds=inner_folds, seed=seed + fold * 1000,
            metric=selection_metric, c_grid=c_grid, l1_grid=l1_grid,
        )
        config, tuning = _select_configuration(
            train, test, y_train, train["group_id"].astype(str), args
        )
        model, config, fit_status = _fit_selected_with_fallback(
            train, test, y_train, config, tuning,
            random_state=seed + fold, sample_weight=_class_balanced_weights(y_train),
        )
        features = list(config["features"])
        score = model.predict_proba(test[features].apply(pd.to_numeric, errors="coerce"))[:, 1]
        oof[test_idx] = score
        rows.append({
            "outer_fold": fold, "status": "used", "fit_status": fit_status,
            "cv_method": method, "n_train": len(train_idx), "n_test": len(test_idx),
            "n_test_positive": int(y_test.eq(1).sum()),
            "n_test_negative": int(y_test.eq(0).sum()),
            "feature_panel": config["feature_panel"],
            "features_used": ";".join(features), "C": config["C"],
            "l1_ratio": config.get("l1_ratio"),
            "selection_status": config["selection_status"],
            "inner_cv_score": config.get("inner_cv_score"),
            "average_precision": _safe_metric(average_precision_score, y_test, score),
            "roc_auc": _safe_metric(roc_auc_score, y_test, score),
        })
        tuning.assign(outer_fold=fold).to_csv(
            output_dir / f"tuning_{prefix}_fold{fold}.csv", index=False
        )
    if not np.isfinite(oof).all():
        raise RuntimeError(f"Incomplete nested OOF coverage: {int((~np.isfinite(oof)).sum())} rows")
    return oof, pd.DataFrame(rows)


def add_metric_row(
    rows: list[dict[str, Any]], frame: pd.DataFrame, score_col: str, *,
    model_name: str, dataset_id: str, analysis_variant: str,
    evaluation: str, n_boot: int, seed: int,
    extra: dict[str, Any] | None = None,
) -> None:
    metric = safe_metrics(frame["label"], frame[score_col])
    reference = safe_metrics(frame["label"], frame["score_integrated_adps"])
    bootstrap = group_bootstrap(
        frame, score_col, reference_col="score_integrated_adps",
        n_boot=n_boot, seed=seed,
    )
    row: dict[str, Any] = {
        "evaluation": evaluation, "dataset": dataset_id,
        "analysis_variant": analysis_variant, "model": model_name, **metric,
        "delta_average_precision_vs_integrated": (
            float(metric["average_precision"] - reference["average_precision"])
            if np.isfinite(metric["average_precision"]) else np.nan
        ),
        "delta_roc_auc_vs_integrated": (
            float(metric["roc_auc"] - reference["roc_auc"])
            if np.isfinite(metric["roc_auc"]) else np.nan
        ),
        **bootstrap,
    }
    if extra: row.update(extra)
    rows.append(row)


def run_same_study(
    specs: list[VariantSpec], data: dict[tuple[str, str], pd.DataFrame],
    output_dir: Path, *, outer_folds: int, inner_folds: int,
    min_positive: int, min_negative: int, selection_metric: str,
    c_grid: str, l1_grid: str, bootstrap: int, seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    metric_rows: list[dict[str, Any]] = []
    prediction_rows: list[pd.DataFrame] = []
    selection_rows: list[pd.DataFrame] = []
    baselines = {
        "integrated_adps": "score_integrated_adps",
        "equal_weight_all": "score_equal_weight_all",
        "duplex_structure_only": "score_duplex_structure_only",
        "route_any_max": "score_route_any_max",
    }
    for spec in specs:
        if not spec.same_study_eligible: continue
        d = data[(spec.dataset_id, spec.analysis_variant)].copy()
        if int(d["label"].eq(1).sum()) < min_positive or int(d["label"].eq(0).sum()) < min_negative:
            continue
        for model_name in ("fixed_l2_raw7", "nested_regularized_prespecified"):
            if model_name == "fixed_l2_raw7":
                oof, folds = fixed_grouped_oof(d, outer_folds=outer_folds, seed=seed)
            else:
                oof, folds = nested_grouped_oof(
                    d, outer_folds=outer_folds, inner_folds=inner_folds,
                    seed=seed, selection_metric=selection_metric,
                    c_grid=c_grid, l1_grid=l1_grid, output_dir=output_dir,
                    prefix=f"{spec.dataset_id}_{spec.analysis_variant}_{model_name}",
                )
            pred = d[[
                "pair_id", "label", "group_id", "score_integrated_adps",
                "score_equal_weight_all", "score_duplex_structure_only", "score_route_any_max",
            ]].copy()
            pred["supervised_score"] = oof
            pred["dataset"] = spec.dataset_id
            pred["analysis_variant"] = spec.analysis_variant
            pred["model"] = model_name
            prediction_rows.append(pred)
            add_metric_row(
                metric_rows, pred.rename(columns={"supervised_score": "_score"}), "_score",
                model_name=model_name, dataset_id=spec.dataset_id,
                analysis_variant=spec.analysis_variant,
                evaluation="same_study_nested_grouped_cv", n_boot=bootstrap,
                seed=seed + len(metric_rows) * 101,
            )
            folds["dataset"] = spec.dataset_id
            folds["analysis_variant"] = spec.analysis_variant
            folds["model"] = model_name
            selection_rows.append(folds)
        for name, col in baselines.items():
            add_metric_row(
                metric_rows, d, col, model_name=name,
                dataset_id=spec.dataset_id, analysis_variant=spec.analysis_variant,
                evaluation="same_study_fixed_baseline", n_boot=bootstrap,
                seed=seed + len(metric_rows) * 101,
            )
    return (
        pd.DataFrame(metric_rows),
        pd.concat(prediction_rows, ignore_index=True) if prediction_rows else pd.DataFrame(),
        pd.concat(selection_rows, ignore_index=True) if selection_rows else pd.DataFrame(),
    )


def run_loso(
    specs: list[VariantSpec], data: dict[tuple[str, str], pd.DataFrame],
    output_dir: Path, *, inner_folds: int, selection_metric: str,
    c_grid: str, l1_grid: str, bootstrap: int, seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    representatives = [s for s in specs if s.loso_representative]
    targets = [s for s in specs if s.loso_target]
    families = sorted({s.study_family for s in targets})
    metric_rows: list[dict[str, Any]] = []
    prediction_rows: list[pd.DataFrame] = []
    selection_rows: list[dict[str, Any]] = []
    coefficient_rows: list[dict[str, Any]] = []
    baselines = {
        "integrated_adps": "score_integrated_adps",
        "equal_weight_all": "score_equal_weight_all",
        "duplex_structure_only": "score_duplex_structure_only",
        "route_any_max": "score_route_any_max",
    }
    for family_index, held_family in enumerate(families):
        parts = [
            data[(s.dataset_id, s.analysis_variant)]
            for s in representatives if s.study_family != held_family
        ]
        if not parts: continue
        train = pd.concat(parts, ignore_index=True)
        if train["label"].nunique() < 2: continue
        args = _selection_args(
            inner_folds=inner_folds, seed=seed + family_index * 10000,
            metric=selection_metric, c_grid=c_grid, l1_grid=l1_grid,
        )
        config, tuning = _select_configuration(
            train, train, train["label"].astype(int), train["prefixed_group_id"],
            args, family=train["study_family"],
        )
        weights = _family_class_balanced_weights(train["label"], train["study_family"])
        model, config, fit_status = _fit_selected_with_fallback(
            train, train, train["label"].astype(int), config, tuning,
            random_state=seed + family_index, sample_weight=weights,
        )
        features = list(config["features"])
        tuning.assign(held_out_family=held_family).to_csv(
            output_dir / f"loso_tuning_holdout_{held_family}.csv", index=False
        )
        selection_rows.append({
            "held_out_family": held_family, "n_train": len(train),
            "n_train_families": int(train["study_family"].nunique()),
            "training_families": ";".join(sorted(train["study_family"].unique())),
            "final_fit_status": fit_status, "feature_panel": config["feature_panel"],
            "features_used": ";".join(features), "C": config["C"],
            "l1_ratio": config.get("l1_ratio"),
            "selection_status": config["selection_status"],
            "inner_cv_score": config.get("inner_cv_score"),
        })
        lr = model.named_steps["logistic"]
        for feature, coef in zip(features, lr.coef_[0]):
            coefficient_rows.append({
                "held_out_family": held_family, "feature": feature,
                "standardized_coefficient": float(coef),
                "feature_panel": config["feature_panel"],
                "C": config["C"], "l1_ratio": config.get("l1_ratio"),
            })
        for spec in [x for x in targets if x.study_family == held_family]:
            target = data[(spec.dataset_id, spec.analysis_variant)].copy()
            if target["label"].nunique() < 2: continue
            score = model.predict_proba(
                target[features].apply(pd.to_numeric, errors="coerce")
            )[:, 1]
            pred = target[[
                "pair_id", "label", "group_id", "score_integrated_adps",
                "score_equal_weight_all", "score_duplex_structure_only", "score_route_any_max",
            ]].copy()
            pred["supervised_score"] = score
            pred["dataset"] = spec.dataset_id
            pred["analysis_variant"] = spec.analysis_variant
            pred["held_out_family"] = held_family
            prediction_rows.append(pred)
            add_metric_row(
                metric_rows, pred.rename(columns={"supervised_score": "_score"}), "_score",
                model_name="loso_regularized_prespecified",
                dataset_id=spec.dataset_id, analysis_variant=spec.analysis_variant,
                evaluation="leave_one_study_family_out", n_boot=bootstrap,
                seed=seed + 500000 + len(metric_rows) * 101,
                extra={
                    "held_out_family": held_family,
                    "selected_panel": config["feature_panel"],
                    "selected_C": config["C"],
                    "selected_l1_ratio": config.get("l1_ratio"),
                },
            )
            for name, col in baselines.items():
                add_metric_row(
                    metric_rows, target, col, model_name=name,
                    dataset_id=spec.dataset_id, analysis_variant=spec.analysis_variant,
                    evaluation="loso_target_fixed_baseline", n_boot=bootstrap,
                    seed=seed + 500000 + len(metric_rows) * 101,
                    extra={"held_out_family": held_family},
                )
    return (
        pd.DataFrame(metric_rows),
        pd.concat(prediction_rows, ignore_index=True) if prediction_rows else pd.DataFrame(),
        pd.DataFrame(selection_rows), pd.DataFrame(coefficient_rows),
    )


def run_supervised_benchmark(args) -> None:
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    specs = read_manifest(args.manifest, args.input_root)
    data = {(s.dataset_id, s.analysis_variant): normalize_audit_table(s) for s in specs}
    same_metrics, same_predictions, same_folds = run_same_study(
        specs, data, output_dir,
        outer_folds=int(args.outer_folds), inner_folds=int(args.inner_folds),
        min_positive=int(args.same_study_min_positive),
        min_negative=int(args.same_study_min_negative),
        selection_metric=str(args.selection_metric), c_grid=str(args.c_grid),
        l1_grid=str(args.l1_ratio_grid), bootstrap=int(args.bootstrap),
        seed=int(args.seed),
    )
    loso_metrics, loso_predictions, loso_selections, coefficients = run_loso(
        specs, data, output_dir, inner_folds=int(args.inner_folds),
        selection_metric=str(args.selection_metric), c_grid=str(args.c_grid),
        l1_grid=str(args.l1_ratio_grid), bootstrap=int(args.bootstrap),
        seed=int(args.seed),
    )
    outputs = {
        "same_study_nested_metrics.csv": same_metrics,
        "same_study_oof_predictions.csv": same_predictions,
        "same_study_fold_selections.csv": same_folds,
        "leave_one_study_family_out_metrics.csv": loso_metrics,
        "leave_one_study_family_out_predictions.csv": loso_predictions,
        "leave_one_study_family_out_model_selections.csv": loso_selections,
        "leave_one_study_family_out_coefficients.csv": coefficients,
    }
    for name, table in outputs.items():
        table.to_csv(output_dir / name, index=False)
    metadata = {
        "manifest": str(args.manifest), "input_root": str(args.input_root),
        "output_dir": str(output_dir), "outer_folds": int(args.outer_folds),
        "inner_folds": int(args.inner_folds),
        "selection_metric": str(args.selection_metric),
        "c_grid": str(args.c_grid), "l1_ratio_grid": str(args.l1_ratio_grid),
        "bootstrap": int(args.bootstrap), "seed": int(args.seed),
        "same_study_min_positive": int(args.same_study_min_positive),
        "same_study_min_negative": int(args.same_study_min_negative),
        "feature_panels": {"raw7": RAW7_FEATURES, "routes6": ROUTE6_FEATURES},
        "route_definitions": ROUTE_DEFINITIONS,
        "solver_policy": {
            "l1_ratio_0": "lbfgs_l2", "l1_ratio_gt_0": "saga_elasticnet",
            "inner_selection_requires_all_valid_folds": True,
            "final_fit_fallback": "best_inner_evaluated_l2_then_prespecified_raw7_l2",
        },
    }
    (output_dir / "run_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    report = [
        "# dsRNASeeker supervised benchmark report", "",
        "- Same-study estimates use true outer grouped cross-validation.",
        "- Feature-panel and regularization selection occur only inside each outer-training fold.",
        "- Leave-one-study-family-out evaluation removes the complete target family before model selection.",
        "- Supervised outputs are ranking scores, not calibrated probabilities.",
        f"- Manifest variants: {len(specs)}",
        f"- Same-study metric rows: {len(same_metrics)}",
        f"- Leave-one-family-out metric rows: {len(loso_metrics)}",
    ]
    (output_dir / "SUPERVISED_BENCHMARK_REPORT.md").write_text("\n".join(report) + "\n")
    print(f"[OK] supervised benchmark completed: {output_dir}")
