from __future__ import annotations

"""Leakage-aware supervised ranking for dsRNASeeker.

This module is a drop-in replacement for modules/supervised.py.  It keeps the
existing apply_supervised_priority() interface, but changes the label and
validation behavior in ways that are important for RIP-seq benchmarking:

* labels may be 1 (positive), 0 (confident negative), or missing (unlabeled);
* unlabeled candidates are never silently converted to negatives;
* a --training-labels table may be either target-specific pair labels or an
  external training matrix containing label plus ADPS feature columns;
* cross-validation is grouped by gene/locus to reduce candidate-level leakage;
* final ranking may therefore be trained on other RIP-seq datasets and applied
  to the current dataset;
* tier display fields and tier-aware sorting are recomputed after prediction.
"""

from pathlib import Path
from typing import Any, Iterable
import json
import re
import warnings

import numpy as np
import pandas as pd

try:
    from sklearn.pipeline import Pipeline
    from sklearn.impute import SimpleImputer
    from sklearn.preprocessing import StandardScaler
    from sklearn.linear_model import LogisticRegression
    from sklearn.exceptions import ConvergenceWarning
    from sklearn.metrics import (
        average_precision_score,
        balanced_accuracy_score,
        roc_auc_score,
    )
    from sklearn.model_selection import GroupKFold, GroupShuffleSplit
    try:
        from sklearn.model_selection import StratifiedGroupKFold
    except Exception:  # older scikit-learn
        StratifiedGroupKFold = None
except Exception as e:  # pragma: no cover
    Pipeline = None
    SimpleImputer = None
    StandardScaler = None
    LogisticRegression = None
    ConvergenceWarning = Warning
    average_precision_score = None
    balanced_accuracy_score = None
    roc_auc_score = None
    GroupKFold = None
    GroupShuffleSplit = None
    StratifiedGroupKFold = None
    _SKLEARN_IMPORT_ERROR = e
else:
    _SKLEARN_IMPORT_ERROR = None


RAW7_FEATURES = [
    "orientation_adps",
    "annotation_adps",
    "case_expression_adps",
    "energy_adps",
    "interface_adps",
    "case_editing_adps",
    "RI_adps",
]

# Backward-compatible superset. control_RI_fraction is optional and is not
# silently added to the manuscript's pre-specified raw7 panel.
SUPERVISED_FEATURES = [*RAW7_FEATURES, "control_RI_fraction"]

ROUTE_DEFINITIONS = {
    "score_route_substrate_context": [
        "orientation_adps", "annotation_adps", "case_expression_adps",
    ],
    "score_route_inverted_structure": [
        "orientation_adps", "energy_adps", "interface_adps",
    ],
    "score_route_condition_structure": [
        "case_expression_adps", "energy_adps", "interface_adps",
    ],
    "score_route_editing_marked": [
        "energy_adps", "interface_adps", "case_editing_adps",
    ],
    "score_route_intron_persistence": [
        "annotation_adps", "case_expression_adps", "RI_adps",
        "energy_adps", "interface_adps",
    ],
    "score_route_context_structure": [
        "annotation_adps", "case_expression_adps", "energy_adps",
        "interface_adps",
    ],
}
ROUTE6_FEATURES = list(ROUTE_DEFINITIONS)

# Dataset-agnostic feature panels. auto_prespecified selects only the two
# panels prospectively evaluated in the final benchmark: raw evidence versus
# pre-specified route composites. Older named panels remain available for
# backward compatibility, but are not searched by auto_prespecified.
FEATURE_PANELS = {
    "raw7": list(RAW7_FEATURES),
    "routes6": list(ROUTE6_FEATURES),
    "orientation_only": ["orientation_adps"],
    "annotation_only": ["annotation_adps"],
    "orientation_annotation": ["orientation_adps", "annotation_adps"],
    "structure_only": ["energy_adps", "interface_adps"],
    "condition_only": [
        "case_expression_adps", "case_editing_adps", "RI_adps",
        "control_RI_fraction",
    ],
    "structure_condition": [
        "case_expression_adps", "energy_adps", "interface_adps",
        "case_editing_adps", "RI_adps", "control_RI_fraction",
    ],
    "no_orientation_annotation": [
        "case_expression_adps", "energy_adps", "interface_adps",
        "case_editing_adps", "RI_adps", "control_RI_fraction",
    ],
    "annotation_structure_condition": [
        "annotation_adps", "case_expression_adps", "energy_adps",
        "interface_adps", "case_editing_adps", "RI_adps",
        "control_RI_fraction",
    ],
    "compact_v1": ["annotation_adps", "case_expression_adps", "energy_adps"],
    "all_current": list(SUPERVISED_FEATURES),
}

PRESPECIFIED_PANEL_ORDER = ["raw7", "routes6"]

SAGA_MAX_ITER = 100000
SAGA_TOL = 1e-3
LBFGS_MAX_ITER = 20000
LBFGS_TOL = 1e-7

TIER_ORDER = {
    "tier1_strict_high": 0,
    "tier2_strict": 1,
    "tier3_relaxed": 2,
    "not_prioritized": 3,
}


def _split_table_spec(path_or_spec: str | Path) -> tuple[Path, str | int | None]:
    """Support '/path/file.xlsx::sheet name' without changing CLI arguments."""
    text = str(path_or_spec)
    if "::" in text:
        path_text, sheet_text = text.rsplit("::", 1)
        sheet: str | int = int(sheet_text) if sheet_text.isdigit() else sheet_text
        return Path(path_text), sheet
    return Path(text), None


def _read_table(path_or_spec: str | Path) -> pd.DataFrame:
    path, sheet = _split_table_spec(path_or_spec)
    suffix = path.suffix.lower()
    if suffix in {".xlsx", ".xls", ".xlsm"}:
        return pd.read_excel(path, sheet_name=0 if sheet is None else sheet)
    if suffix in {".tsv", ".txt", ".bed"}:
        return pd.read_csv(path, sep="\t")
    return pd.read_csv(path, sep=None, engine="python")


def _clean_symbol(x: Any) -> str:
    if pd.isna(x):
        return ""
    value = str(x).strip()
    if value.lower() in {"", "nan", "none", "na", "n/a"}:
        return ""
    return value


def _symbol_key(x: Any) -> str:
    return _clean_symbol(x).upper()


def _ensure_symbol_col(T: pd.DataFrame, requested: str, source: str | Path) -> str:
    if requested in T.columns:
        return requested
    candidates = [
        "Symbol", "symbol", "gene_symbol", "gene_name", "Gene", "gene",
        "external_gene_name", "GeneSymbol",
    ]
    for col in candidates:
        if col in T.columns:
            return col
    raise ValueError(
        f"Could not find gene-symbol column '{requested}' in {source}. "
        f"Available columns: {list(T.columns)}"
    )


def _numeric_label(series: pd.Series) -> pd.Series:
    """Return nullable labels; invalid/missing values remain unlabeled."""
    out = pd.to_numeric(series, errors="coerce")
    bad = out.notna() & ~out.isin([0, 1])
    if bad.any():
        vals = sorted(set(out.loc[bad].tolist()))
        raise ValueError(f"Labels must be 0, 1, or missing; found {vals[:10]}")
    return out.astype("Float64")


def derive_pair_labels_from_truth_table(
    M: pd.DataFrame,
    truth_table: str | Path,
    *,
    symbol_col: str = "Symbol",
    truth_label_mode: str = "positive_logfc_padj",
    truth_label_col: str | None = None,
    padj_col: str = "padj",
    logfc_col: str = "log2FoldChange",
    padj_max: float = 0.05,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Convert gene-level evidence to conservative pair labels.

    Pair labels are deliberately ternary:
      1  at least one arm gene is a positive gene;
      0  every annotated arm gene is a *confident negative* gene;
      NA otherwise (unknown, not a negative).

    For positive_logfc_padj, significant positive genes are positives,
    significant negative genes are confident negatives, and nonsignificant
    genes are unlabeled.  This avoids the previous positive-vs-everything-else
    mislabeling.
    """
    T = _read_table(truth_table).copy()
    symbol_col = _ensure_symbol_col(T, symbol_col, truth_table)
    T["truth_symbol"] = T[symbol_col].map(_clean_symbol)
    T["truth_symbol_key"] = T["truth_symbol"].map(_symbol_key)
    T = T[T["truth_symbol_key"].ne("")].copy()

    mode = str(truth_label_mode)
    T["truth_padj"] = np.nan
    T["truth_log2FoldChange"] = np.nan
    T["truth_gene_label"] = pd.Series(pd.NA, index=T.index, dtype="Float64")

    if mode == "all_table_rows":
        T["truth_gene_label"] = 1.0
        T["truth_positive_rule"] = "all non-empty truth-table symbols"

    elif mode == "explicit_label_col":
        if not truth_label_col:
            raise ValueError("--truth-label-col is required for explicit_label_col")
        if truth_label_col not in T.columns:
            raise ValueError(
                f"Could not find label column '{truth_label_col}' in {truth_table}. "
                f"Available columns: {list(T.columns)}"
            )
        T["truth_gene_label"] = _numeric_label(T[truth_label_col])
        T["truth_positive_rule"] = f"explicit ternary label column: {truth_label_col}"

    elif mode == "padj_only":
        if padj_col not in T.columns:
            raise ValueError(f"Missing adjusted-P column '{padj_col}' in {truth_table}")
        T["truth_padj"] = pd.to_numeric(T[padj_col], errors="coerce")
        T.loc[T["truth_padj"].le(float(padj_max)), "truth_gene_label"] = 1.0
        T["truth_positive_rule"] = f"{padj_col}<={padj_max}; all other genes unlabeled"
        if logfc_col in T.columns:
            T["truth_log2FoldChange"] = pd.to_numeric(T[logfc_col], errors="coerce")

    elif mode == "positive_logfc_padj":
        missing = [c for c in (padj_col, logfc_col) if c not in T.columns]
        if missing:
            raise ValueError(
                f"Missing columns {missing} for positive_logfc_padj in {truth_table}"
            )
        T["truth_padj"] = pd.to_numeric(T[padj_col], errors="coerce")
        T["truth_log2FoldChange"] = pd.to_numeric(T[logfc_col], errors="coerce")
        significant = T["truth_padj"].le(float(padj_max))
        T.loc[significant & T["truth_log2FoldChange"].gt(0), "truth_gene_label"] = 1.0
        T.loc[significant & T["truth_log2FoldChange"].lt(0), "truth_gene_label"] = 0.0
        T["truth_positive_rule"] = (
            f"positive: {padj_col}<={padj_max} and {logfc_col}>0; "
            f"negative: {padj_col}<={padj_max} and {logfc_col}<0; otherwise unlabeled"
        )
    else:
        raise ValueError(f"Unknown truth_label_mode: {mode}")

    # Resolve duplicate symbols conservatively: positive dominates; otherwise
    # negative only if there is an explicit negative and no positive.
    symbol_status: dict[str, float] = {}
    for key, grp in T.groupby("truth_symbol_key", sort=False):
        vals = set(pd.to_numeric(grp["truth_gene_label"], errors="coerce").dropna().astype(int))
        if 1 in vals:
            symbol_status[key] = 1.0
        elif 0 in vals:
            symbol_status[key] = 0.0

    labels: list[dict[str, Any]] = []
    for _, row in M.iterrows():
        a = _clean_symbol(row.get("A_SYMBOL", ""))
        b = _clean_symbol(row.get("B_SYMBOL", ""))
        keys = sorted({_symbol_key(x) for x in (a, b) if _symbol_key(x)})
        states = [symbol_status.get(k, np.nan) for k in keys]
        positive = [k for k, state in zip(keys, states) if state == 1]
        negative = [k for k, state in zip(keys, states) if state == 0]

        if positive:
            label: float | None = 1.0
            source = f"truth_table_positive:{mode}"
        elif keys and len(negative) == len(keys):
            label = 0.0
            source = f"truth_table_confident_negative:{mode}"
        else:
            label = None
            source = f"truth_table_unlabeled:{mode}"

        labels.append({
            "pair_id": row.get("pair_id", ""),
            "label": label,
            "label_source": source,
            "matched_positive_symbols": ";".join(positive),
            "matched_negative_symbols": ";".join(negative),
            "A_SYMBOL": a,
            "B_SYMBOL": b,
        })

    labels_df = pd.DataFrame(labels)
    truth_cols = [
        c for c in [
            "truth_symbol", "truth_symbol_key", "truth_gene_label",
            "truth_positive_rule", "truth_padj", "truth_log2FoldChange",
            symbol_col,
        ] if c in T.columns
    ]
    truth_full = T[truth_cols].copy()
    truth_positive = truth_full.loc[
        pd.to_numeric(truth_full["truth_gene_label"], errors="coerce").eq(1)
    ].copy()
    return labels_df, truth_full, truth_positive


def read_pair_labels(
    labels_path: str | Path,
    M: pd.DataFrame,
) -> tuple[pd.DataFrame, None, None]:
    """Read pair labels or an external labeled feature matrix.

    Required columns are label and, for target-specific labels, pair_id.  When
    two or more model feature columns are present, the file is treated as an
    external training matrix and does not need to share pair IDs with M.
    """
    L = _read_table(labels_path).copy()
    if "label" not in L.columns:
        raise ValueError(f"Training-labels table must contain a 'label' column: {labels_path}")
    L["label"] = _numeric_label(L["label"])
    external_features = [c for c in SUPERVISED_FEATURES if c in L.columns]
    if len(external_features) < 2 and "pair_id" not in L.columns:
        raise ValueError(
            "Target-specific labels need pair_id; external training matrices need "
            "label plus at least two ADPS feature columns."
        )
    if "label_source" not in L.columns:
        L["label_source"] = "pair_level_label_file"
    return L, None, None


def _model(
    random_state: int,
    *,
    model_type: str = "legacy_l2",
    C: float = 1.0,
    l1_ratio: float | None = None,
) -> Pipeline:
    """Build a scaled logistic ranking model with a stable solver policy.

    Pure L2 configurations use LBFGS. Genuine L1/Elastic-Net configurations
    use SAGA. This mirrors the final manuscript benchmark and avoids using
    SAGA unnecessarily when l1_ratio is zero.
    """
    if _SKLEARN_IMPORT_ERROR is not None:
        raise ImportError(
            "scikit-learn is required for supervised ranking"
        ) from _SKLEARN_IMPORT_ERROR
    if float(C) <= 0:
        raise ValueError(f"C must be positive; received {C}")
    model_type = str(model_type).lower()
    if model_type not in {"legacy_l2", "elasticnet"}:
        raise ValueError(
            f"Unknown supervised model '{model_type}'. Choose legacy_l2 or elasticnet."
        )
    ratio = None if l1_ratio is None else float(l1_ratio)
    if ratio is not None and not 0.0 <= ratio <= 1.0:
        raise ValueError(f"l1_ratio must be within [0, 1]; received {ratio}")
    use_l2 = model_type == "legacy_l2" or ratio is None or np.isclose(ratio, 0.0)
    if use_l2:
        logistic = LogisticRegression(
            penalty="l2", solver="lbfgs", C=float(C),
            max_iter=LBFGS_MAX_ITER, tol=LBFGS_TOL,
            random_state=int(random_state),
        )
    else:
        logistic = LogisticRegression(
            penalty="elasticnet", solver="saga", l1_ratio=ratio, C=float(C),
            max_iter=SAGA_MAX_ITER, tol=SAGA_TOL,
            random_state=int(random_state),
        )
    return Pipeline([
        ("imputer", SimpleImputer(strategy="median")),
        ("scaler", StandardScaler()),
        ("logistic", logistic),
    ])

def _class_balanced_weights(y: pd.Series) -> np.ndarray:
    y = pd.Series(y).astype(int).reset_index(drop=True)
    counts = y.value_counts()
    if len(counts) < 2:
        return np.ones(len(y), dtype=float)
    weights = y.map({c: 1.0 / (len(counts) * counts[c]) for c in counts.index}).to_numpy(float)
    return weights / np.mean(weights)


def _family_class_balanced_weights(y: pd.Series, family: pd.Series) -> np.ndarray:
    y = pd.Series(y).astype(int).reset_index(drop=True)
    family = pd.Series(family).fillna("").astype(str).reset_index(drop=True)
    families = sorted(x for x in family.unique() if x and x != "nan")
    if len(families) < 2:
        return _class_balanced_weights(y)
    out = np.zeros(len(y), dtype=float)
    for fam in families:
        fam_mask = family.eq(fam)
        for cls in (0, 1):
            mask = fam_mask & y.eq(cls)
            n = int(mask.sum())
            if n:
                out[mask.to_numpy()] = 1.0 / (len(families) * 2.0 * n)
    if np.any(out <= 0):
        return _class_balanced_weights(y)
    return out / np.mean(out)


def _fit_checked(
    model: Pipeline,
    X: pd.DataFrame,
    y: pd.Series,
    sample_weight: np.ndarray | None = None,
) -> tuple[Pipeline, bool]:
    """Fit and return an explicit convergence status."""
    weights = _class_balanced_weights(y) if sample_weight is None else np.asarray(sample_weight, float)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", ConvergenceWarning)
        model.fit(X, y, logistic__sample_weight=weights)
    converged = not any(issubclass(w.category, ConvergenceWarning) for w in caught)
    return model, converged


def ensure_route_features(df: pd.DataFrame) -> pd.DataFrame:
    """Add the six pre-specified route scores and their maximum when absent."""
    out = df.copy()
    for route, columns in ROUTE_DEFINITIONS.items():
        if route in out.columns:
            out[route] = pd.to_numeric(out[route], errors="coerce")
            continue
        available = [c for c in columns if c in out.columns]
        if not available:
            out[route] = np.nan
            continue
        values = out[available].apply(pd.to_numeric, errors="coerce")
        out[route] = values.mean(axis=1, skipna=True)
    if "score_route_any_max" not in out.columns:
        out["score_route_any_max"] = out[ROUTE6_FEATURES].max(axis=1, skipna=True)
    return out

def _safe_metric(fn, y_true: Iterable[Any], y_score_or_pred: Iterable[Any]) -> float:
    try:
        y = pd.Series(y_true).astype(int)
        if y.nunique() < 2:
            return np.nan
        return float(fn(y, y_score_or_pred))
    except Exception:
        return np.nan


def _derive_group_id(df: pd.DataFrame) -> pd.Series:
    if "group_id" in df.columns:
        given = df["group_id"].astype(str).replace({"": np.nan, "nan": np.nan})
    else:
        given = pd.Series(np.nan, index=df.index, dtype=object)

    a = df.get("A_SYMBOL", pd.Series("", index=df.index)).map(_symbol_key)
    b = df.get("B_SYMBOL", pd.Series("", index=df.index)).map(_symbol_key)
    gene_group = pd.Series(index=df.index, dtype=object)
    for idx in df.index:
        symbols = sorted({x for x in (a.loc[idx], b.loc[idx]) if x})
        gene_group.loc[idx] = "GENE:" + "|".join(symbols) if symbols else ""

    pair = df.get("pair_id", pd.Series("", index=df.index)).astype(str)
    locus_group = pair.str.extract(r"_(chr[^_]+)_([0-9]+)_[0-9]+_[+-]", expand=True)
    locus = pd.Series("", index=df.index, dtype=object)
    if locus_group.shape[1] == 2:
        chrom = locus_group[0].fillna("")
        start = pd.to_numeric(locus_group[1], errors="coerce")
        # 100-kb locus bins reduce leakage among many nearby TE pairs.
        locus = "LOCUS:" + chrom + ":" + (start // 100000).fillna(-1).astype(int).astype(str)
        locus = locus.where(chrom.ne(""), "")

    fallback = pd.Series([f"ROW:{i}" for i in range(len(df))], index=df.index)
    out = given.fillna("")
    out = out.where(out.ne(""), gene_group)
    out = out.where(out.ne(""), locus)
    out = out.where(out.ne(""), fallback)
    return out.astype(str)


def _nonconstant_features(
    training: pd.DataFrame,
    target: pd.DataFrame,
    candidates: Iterable[str] | None = None,
) -> list[str]:
    features: list[str] = []
    for col in list(SUPERVISED_FEATURES if candidates is None else candidates):
        if col not in training.columns or col not in target.columns:
            continue
        values = pd.to_numeric(training[col], errors="coerce")
        if values.notna().sum() >= 2 and values.nunique(dropna=True) > 1:
            features.append(col)
    return features


def _parse_float_grid(value: Any, default: list[float]) -> list[float]:
    if value is None:
        return list(default)
    if isinstance(value, (list, tuple, np.ndarray)):
        raw = list(value)
    else:
        raw = [x.strip() for x in str(value).split(",") if x.strip()]
    parsed: list[float] = []
    for item in raw:
        number = float(item)
        if not np.isfinite(number):
            raise ValueError(f"Non-finite hyperparameter value: {item}")
        parsed.append(number)
    if not parsed:
        raise ValueError("Hyperparameter grid cannot be empty")
    return sorted(set(parsed))


def _metric_from_name(name: str, y_true: pd.Series, probability: np.ndarray) -> float:
    name = str(name)
    if y_true.nunique() < 2:
        return np.nan
    if name == "average_precision":
        return _safe_metric(average_precision_score, y_true, probability)
    if name == "roc_auc":
        return _safe_metric(roc_auc_score, y_true, probability)
    if name == "balanced_accuracy":
        return _safe_metric(
            balanced_accuracy_score, y_true, (probability >= 0.5).astype(int)
        )
    raise ValueError(f"Unknown supervised selection metric: {name}")


def _select_configuration(
    training: pd.DataFrame,
    target: pd.DataFrame,
    y: pd.Series,
    groups: pd.Series,
    args,
    *,
    family: pd.Series | None = None,
) -> tuple[dict[str, Any], pd.DataFrame]:
    """Select a pre-specified panel and regularization inside training data.

    A configuration is eligible only if every valid inner fold converges and
    produces the selected metric. Target labels are never used.
    """
    training = ensure_route_features(training)
    target = ensure_route_features(target)
    requested_panel = str(getattr(args, "supervised_feature_panel", "raw7"))
    model_type = str(getattr(args, "supervised_model", "legacy_l2"))
    tune = bool(getattr(args, "supervised_tune", False))
    if requested_panel == "auto_prespecified":
        tune = True
        panels = list(PRESPECIFIED_PANEL_ORDER)
        model_type = "elasticnet"
    else:
        if requested_panel not in FEATURE_PANELS:
            raise ValueError(
                f"Unknown supervised feature panel '{requested_panel}'. "
                f"Available: {sorted(FEATURE_PANELS)} plus auto_prespecified"
            )
        panels = [requested_panel]

    fixed_c = float(getattr(args, "supervised_c", 1.0))
    fixed_l1 = float(getattr(args, "supervised_l1_ratio", 0.5))
    if tune:
        c_values = _parse_float_grid(
            getattr(args, "supervised_c_grid", "0.03,0.3,3"), [0.03, 0.3, 3.0]
        )
        if any(c <= 0 for c in c_values):
            raise ValueError("All C values must be positive")
        l1_values = (
            _parse_float_grid(
                getattr(args, "supervised_l1_ratio_grid", "0,0.5,1"),
                [0.0, 0.5, 1.0],
            ) if model_type == "elasticnet" else [None]
        )
    else:
        c_values = [fixed_c]
        l1_values = [fixed_l1 if model_type == "elasticnet" else None]

    metric_name = str(getattr(args, "supervised_selection_metric", "average_precision"))
    inner_folds = int(getattr(args, "supervised_inner_cv_folds", 3))
    random_state = int(getattr(args, "supervised_random_state", 1))
    y = pd.Series(y).astype(int).reset_index(drop=True)
    groups = pd.Series(groups).astype(str).reset_index(drop=True)
    family_series = None if family is None else pd.Series(family).astype(str).reset_index(drop=True)
    splits, split_method = _choose_cv_splits(
        training, y, groups, cv_folds=inner_folds, random_state=random_state
    )
    expected_folds = sum(
        int(y.iloc[tr].nunique() == 2 and y.iloc[te].nunique() == 2)
        for tr, te in splits
    )

    rows: list[dict[str, Any]] = []
    for panel_index, panel in enumerate(panels):
        features = _nonconstant_features(training, target, FEATURE_PANELS[panel])
        if not features:
            rows.append({
                "panel": panel, "panel_order": panel_index,
                "model_type": model_type,
                "status": "skipped_no_nonconstant_shared_features",
                "features_used": "", "n_features": 0,
                "inner_cv_method": split_method,
                "inner_cv_folds_expected": expected_folds,
            })
            continue
        X = training[features].apply(pd.to_numeric, errors="coerce").reset_index(drop=True)
        for C in c_values:
            for l1_ratio in l1_values:
                fold_scores: list[float] = []
                used_folds = 0
                nonconverged_folds = 0
                if tune and splits:
                    for fold_number, (train_idx, test_idx) in enumerate(splits, start=1):
                        y_train = y.iloc[train_idx]
                        y_test = y.iloc[test_idx]
                        if y_train.nunique() < 2 or y_test.nunique() < 2:
                            continue
                        weights = (
                            _family_class_balanced_weights(y_train, family_series.iloc[train_idx])
                            if family_series is not None else _class_balanced_weights(y_train)
                        )
                        model = _model(
                            random_state + fold_number, model_type=model_type,
                            C=C, l1_ratio=l1_ratio,
                        )
                        model, converged = _fit_checked(
                            model, X.iloc[train_idx], y_train, weights
                        )
                        if not converged:
                            nonconverged_folds += 1
                            continue
                        probability = model.predict_proba(X.iloc[test_idx])[:, 1]
                        score = _metric_from_name(metric_name, y_test, probability)
                        if np.isfinite(score):
                            fold_scores.append(float(score))
                            used_folds += 1
                mean_score = float(np.mean(fold_scores)) if fold_scores else np.nan
                sd_score = float(np.std(fold_scores, ddof=1)) if len(fold_scores) > 1 else np.nan
                fully_evaluated = (
                    not tune or (
                        expected_folds > 0 and used_folds == expected_folds
                        and nonconverged_folds == 0
                    )
                )
                rows.append({
                    "panel": panel, "panel_order": panel_index,
                    "model_type": model_type, "C": float(C),
                    "l1_ratio": np.nan if l1_ratio is None else float(l1_ratio),
                    "status": "evaluated" if fully_evaluated else (
                        "skipped_nonconverged_inner_fit" if nonconverged_folds
                        else "skipped_incomplete_inner_folds"
                    ),
                    "features_used": ";".join(features), "n_features": len(features),
                    "selection_metric": metric_name,
                    "inner_cv_method": split_method,
                    "inner_cv_folds_expected": expected_folds,
                    "inner_cv_folds_used": used_folds,
                    "inner_cv_nonconverged_folds": nonconverged_folds,
                    "mean_selection_score": mean_score,
                    "sd_selection_score": sd_score,
                    "fold_scores": ";".join(f"{x:.10g}" for x in fold_scores),
                })

    results = pd.DataFrame(rows)
    usable = results.loc[results["status"].eq("evaluated")].copy()
    if usable.empty:
        features = _nonconstant_features(training, target, RAW7_FEATURES)
        if not features:
            raise ValueError("No supervised configuration had usable shared features")
        return ({
            "model_type": "elasticnet", "feature_panel": "raw7",
            "features": features, "C": 0.3, "l1_ratio": 0.0,
            "selection_status": "fallback_stable_l2_no_fully_evaluated_configuration",
            "selection_metric": metric_name, "inner_cv_method": split_method,
            "inner_cv_folds_requested": inner_folds,
            "inner_cv_folds_expected": expected_folds,
            "inner_cv_score": None, "tuning_enabled": tune,
        }, results)

    if tune and usable["mean_selection_score"].notna().any():
        ranked = usable.loc[usable["mean_selection_score"].notna()].copy()
        ranked["_l1_tie"] = pd.to_numeric(ranked["l1_ratio"], errors="coerce").fillna(-1.0)
        ranked = ranked.sort_values(
            ["mean_selection_score", "n_features", "C", "_l1_tie", "panel_order"],
            ascending=[False, True, True, False, True], kind="mergesort",
        )
        chosen = ranked.iloc[0]
        selection_status = "selected_by_inner_grouped_cv"
    else:
        chosen = usable.iloc[0]
        selection_status = "fixed_configuration"
    l1_value = pd.to_numeric(pd.Series([chosen.get("l1_ratio")]), errors="coerce").iloc[0]
    return ({
        "model_type": str(chosen["model_type"]),
        "feature_panel": str(chosen["panel"]),
        "features": str(chosen["features_used"]).split(";") if str(chosen["features_used"]) else [],
        "C": float(chosen["C"]),
        "l1_ratio": None if pd.isna(l1_value) else float(l1_value),
        "selection_status": selection_status,
        "selection_metric": metric_name,
        "inner_cv_method": split_method,
        "inner_cv_folds_requested": inner_folds,
        "inner_cv_folds_expected": expected_folds,
        "inner_cv_score": None if pd.isna(chosen.get("mean_selection_score")) else float(chosen["mean_selection_score"]),
        "tuning_enabled": tune,
    }, results)

def _best_l2_configuration(
    training: pd.DataFrame,
    target: pd.DataFrame,
    tuning: pd.DataFrame,
) -> dict[str, Any]:
    if not tuning.empty:
        z = tuning.loc[
            tuning["status"].eq("evaluated")
            & pd.to_numeric(tuning["l1_ratio"], errors="coerce").fillna(0.0).eq(0.0)
            & tuning["mean_selection_score"].notna()
        ].copy()
        if not z.empty:
            row = z.sort_values(
                ["mean_selection_score", "n_features", "C", "panel_order"],
                ascending=[False, True, True, True], kind="mergesort",
            ).iloc[0]
            return {
                "model_type": "elasticnet", "feature_panel": str(row["panel"]),
                "features": str(row["features_used"]).split(";"),
                "C": float(row["C"]), "l1_ratio": 0.0,
                "selection_status": "inner_selected_l2_convergence_fallback",
            }
    features = _nonconstant_features(training, target, RAW7_FEATURES)
    return {
        "model_type": "elasticnet", "feature_panel": "raw7", "features": features,
        "C": 0.3, "l1_ratio": 0.0,
        "selection_status": "prespecified_raw7_l2_convergence_fallback",
    }


def _fit_selected_with_fallback(
    training: pd.DataFrame,
    target: pd.DataFrame,
    y: pd.Series,
    config: dict[str, Any],
    tuning: pd.DataFrame,
    *,
    random_state: int,
    sample_weight: np.ndarray | None = None,
) -> tuple[Pipeline, dict[str, Any], str]:
    features = list(config["features"])
    model = _model(
        random_state, model_type=config["model_type"],
        C=config["C"], l1_ratio=config.get("l1_ratio"),
    )
    model, converged = _fit_checked(
        model, training[features].apply(pd.to_numeric, errors="coerce"),
        y, sample_weight,
    )
    if converged:
        return model, config, "selected_configuration_converged"
    fallback = _best_l2_configuration(training, target, tuning)
    model = _model(
        random_state, model_type="elasticnet", C=fallback["C"], l1_ratio=0.0
    )
    model, converged = _fit_checked(
        model, training[fallback["features"]].apply(pd.to_numeric, errors="coerce"),
        y, sample_weight,
    )
    if not converged:
        raise RuntimeError("Selected configuration and training-only L2 fallback did not converge")
    return model, fallback, "selected_configuration_nonconverged_used_training_only_l2_fallback"


def _nested_grouped_cv(
    training: pd.DataFrame,
    target: pd.DataFrame,
    y: pd.Series,
    groups: pd.Series,
    args,
) -> tuple[dict[str, Any], pd.DataFrame, pd.DataFrame]:
    """True outer grouped CV with selection repeated inside each outer fold."""
    training = ensure_route_features(training).reset_index(drop=True)
    target = ensure_route_features(target)
    y = pd.Series(y).astype(int).reset_index(drop=True)
    groups = pd.Series(groups).astype(str).reset_index(drop=True)
    outer_folds = int(getattr(args, "cv_folds", 5))
    random_state = int(getattr(args, "supervised_random_state", 1))
    splits, method = _choose_cv_splits(
        training, y, groups, cv_folds=outer_folds, random_state=random_state
    )
    if not splits:
        return ({
            "nested_cv_enabled": False, "nested_cv_method": method,
            "nested_cv_folds_used": 0,
        }, pd.DataFrame(), pd.DataFrame())
    oof = np.full(len(y), np.nan, dtype=float)
    fold_id = np.full(len(y), -1, dtype=int)
    rows: list[dict[str, Any]] = []
    for fold, (train_idx, test_idx) in enumerate(splits, start=1):
        train = training.iloc[train_idx].reset_index(drop=True)
        test = training.iloc[test_idx].reset_index(drop=True)
        y_train = y.iloc[train_idx].reset_index(drop=True)
        y_test = y.iloc[test_idx].reset_index(drop=True)
        if y_train.nunique() < 2 or y_test.nunique() < 2:
            rows.append({"outer_fold": fold, "status": "skipped_single_class"})
            continue
        config, tuning = _select_configuration(
            train, test, y_train,
            groups.iloc[train_idx].reset_index(drop=True), args,
        )
        model, config, fit_status = _fit_selected_with_fallback(
            train, test, y_train, config, tuning,
            random_state=random_state + fold,
            sample_weight=_class_balanced_weights(y_train),
        )
        features = list(config["features"])
        probability = model.predict_proba(
            test[features].apply(pd.to_numeric, errors="coerce")
        )[:, 1]
        oof[test_idx] = probability
        fold_id[test_idx] = fold
        rows.append({
            "outer_fold": fold, "status": "used", "fit_status": fit_status,
            "n_train": len(train_idx), "n_test": len(test_idx),
            "n_test_positive": int(y_test.eq(1).sum()),
            "n_test_negative": int(y_test.eq(0).sum()),
            "feature_panel": config["feature_panel"],
            "features_used": ";".join(features), "C": config["C"],
            "l1_ratio": config.get("l1_ratio"),
            "selection_status": config["selection_status"],
            "inner_cv_score": config.get("inner_cv_score"),
            "average_precision": _safe_metric(average_precision_score, y_test, probability),
            "roc_auc": _safe_metric(roc_auc_score, y_test, probability),
        })
    valid = np.isfinite(oof)
    pred = pd.DataFrame({
        "row_index": np.arange(len(y)), "group_id": groups.to_numpy(),
        "label": y.to_numpy(), "outer_fold": fold_id, "oof_score": oof,
    })
    fold_df = pd.DataFrame(rows)
    if not valid.all():
        return ({
            "nested_cv_enabled": False, "nested_cv_method": method,
            "nested_cv_reason": f"incomplete_oof_coverage:{int((~valid).sum())}",
            "nested_cv_folds_used": int(fold_df.get("status", pd.Series(dtype=str)).eq("used").sum()),
        }, pred, fold_df)
    return ({
        "nested_cv_enabled": True, "nested_cv_method": method,
        "nested_cv_folds_requested": outer_folds,
        "nested_cv_folds_used": int(fold_df.get("status", pd.Series(dtype=str)).eq("used").sum()),
        "nested_cv_oof_n": int(valid.sum()),
        "nested_cv_oof_average_precision": _safe_metric(average_precision_score, y, oof),
        "nested_cv_oof_roc_auc": _safe_metric(roc_auc_score, y, oof),
        "nested_cv_interpretation": "outer-fold estimate of the complete panel/hyperparameter selection procedure",
    }, pred, fold_df)


def _choose_cv_splits(
    X: pd.DataFrame,
    y: pd.Series,
    groups: pd.Series,
    *,
    cv_folds: int,
    random_state: int,
) -> tuple[list[tuple[np.ndarray, np.ndarray]], str]:
    if cv_folds is None or int(cv_folds) <= 1:
        return [], "disabled"
    y = pd.Series(y).astype(int).reset_index(drop=True)
    groups = pd.Series(groups).astype(str).reset_index(drop=True)
    pos_groups = groups.loc[y.eq(1)].nunique()
    neg_groups = groups.loc[y.eq(0)].nunique()
    n_splits = min(int(cv_folds), int(pos_groups), int(neg_groups), int(groups.nunique()))
    if n_splits < 2:
        return [], "too_few_positive_or_negative_groups"
    if StratifiedGroupKFold is not None:
        for offset in range(100):
            splitter = StratifiedGroupKFold(
                n_splits=n_splits, shuffle=True, random_state=random_state + offset
            )
            splits = list(splitter.split(X, y, groups))
            if all(
                y.iloc[tr].nunique() == 2 and y.iloc[te].nunique() == 2
                for tr, te in splits
            ):
                method = "StratifiedGroupKFold" if offset == 0 else f"StratifiedGroupKFold_seed_offset_{offset}"
                return splits, method
        splitter = StratifiedGroupKFold(
            n_splits=n_splits, shuffle=True, random_state=random_state
        )
        return list(splitter.split(X, y, groups)), "StratifiedGroupKFold_single_class_possible"
    splitter = GroupKFold(n_splits=n_splits)
    return list(splitter.split(X, y, groups)), "GroupKFold_fallback"

def _grouped_cv(
    X: pd.DataFrame,
    y: pd.Series,
    groups: pd.Series,
    *,
    cv_folds: int,
    random_state: int,
    model_config: dict[str, Any],
) -> tuple[dict[str, Any], pd.DataFrame, pd.DataFrame]:
    splits, method = _choose_cv_splits(
        X, y, groups, cv_folds=cv_folds, random_state=random_state
    )
    if not splits:
        return ({
            "cv_enabled": False,
            "cv_method": method,
            "cv_folds_used": 0,
        }, pd.DataFrame(), pd.DataFrame())

    oof = np.full(len(y), np.nan, dtype=float)
    fold_rows: list[dict[str, Any]] = []
    fold_id = np.full(len(y), -1, dtype=int)

    for number, (train_idx, test_idx) in enumerate(splits, start=1):
        y_train = y.iloc[train_idx]
        y_test = y.iloc[test_idx]
        if y_train.nunique() < 2 or y_test.nunique() < 2:
            fold_rows.append({
                "fold": number,
                "status": "skipped_single_class",
                "n_train": len(train_idx),
                "n_test": len(test_idx),
            })
            continue
        model = _model(
            random_state + number,
            model_type=model_config["model_type"],
            C=model_config["C"],
            l1_ratio=model_config.get("l1_ratio"),
        )
        model, converged = _fit_checked(model, X.iloc[train_idx], y_train)
        if not converged:
            fold_rows.append({
                "fold": number,
                "status": "skipped_nonconverged",
                "n_train": len(train_idx),
                "n_test": len(test_idx),
            })
            continue
        prob = model.predict_proba(X.iloc[test_idx])[:, 1]
        pred = (prob >= 0.5).astype(int)
        oof[test_idx] = prob
        fold_id[test_idx] = number
        fold_rows.append({
            "fold": number,
            "status": "used",
            "n_train": len(train_idx),
            "n_test": len(test_idx),
            "n_test_positive": int(y_test.eq(1).sum()),
            "n_test_negative": int(y_test.eq(0).sum()),
            "roc_auc": _safe_metric(roc_auc_score, y_test, prob),
            "average_precision": _safe_metric(average_precision_score, y_test, prob),
            "balanced_accuracy": _safe_metric(balanced_accuracy_score, y_test, pred),
        })

    valid = np.isfinite(oof)
    oof_df = pd.DataFrame({
        "row_index": np.arange(len(y)),
        "group_id": groups.to_numpy(),
        "label": y.to_numpy(),
        "cv_fold": fold_id,
        "oof_probability": oof,
    })
    fold_df = pd.DataFrame(fold_rows)

    if valid.sum() == 0:
        return ({
            "cv_enabled": False,
            "cv_method": method,
            "cv_reason": "all_grouped_folds_were_single_class",
            "cv_folds_used": 0,
        }, oof_df, fold_df)

    y_valid = y.iloc[np.where(valid)[0]]
    p_valid = oof[valid]
    report = {
        "cv_enabled": True,
        "cv_method": method,
        "cv_folds_requested": int(cv_folds),
        "cv_folds_used": int((fold_df.get("status", pd.Series(dtype=str)) == "used").sum()),
        "cv_oof_n": int(valid.sum()),
        "cv_oof_roc_auc": _safe_metric(roc_auc_score, y_valid, p_valid),
        "cv_oof_average_precision": _safe_metric(average_precision_score, y_valid, p_valid),
        "cv_oof_balanced_accuracy": _safe_metric(
            balanced_accuracy_score, y_valid, (p_valid >= 0.5).astype(int)
        ),
    }
    return report, oof_df, fold_df



def _leave_one_dataset_out_cv(
    X: pd.DataFrame,
    y: pd.Series,
    dataset_ids: pd.Series,
    *,
    random_state: int,
    model_config: dict[str, Any],
) -> tuple[dict[str, Any], pd.DataFrame, pd.DataFrame]:
    """Evaluate transfer across studies rather than across rows from one study."""
    dataset_ids = dataset_ids.fillna("").astype(str).reset_index(drop=True)
    valid_dataset = dataset_ids.ne("") & dataset_ids.ne("nan")
    unique = sorted(dataset_ids.loc[valid_dataset].unique().tolist())
    if len(unique) < 2:
        return ({
            "loso_dataset_cv_enabled": False,
            "loso_dataset_cv_reason": "fewer_than_two_dataset_ids",
            "loso_dataset_cv_folds_used": 0,
        }, pd.DataFrame(), pd.DataFrame())

    oof = np.full(len(y), np.nan, dtype=float)
    fold_rows: list[dict[str, Any]] = []
    for number, held_dataset in enumerate(unique, start=1):
        test_mask = dataset_ids.eq(held_dataset).to_numpy()
        train_mask = valid_dataset.to_numpy() & ~test_mask
        train_idx = np.where(train_mask)[0]
        test_idx = np.where(test_mask)[0]
        y_train = y.iloc[train_idx]
        y_test = y.iloc[test_idx]
        row: dict[str, Any] = {
            "fold": number,
            "held_out_dataset": held_dataset,
            "n_train": int(len(train_idx)),
            "n_test": int(len(test_idx)),
            "n_test_positive": int(y_test.eq(1).sum()),
            "n_test_negative": int(y_test.eq(0).sum()),
        }
        if len(train_idx) == 0 or len(test_idx) == 0:
            row["status"] = "skipped_empty"
            fold_rows.append(row)
            continue
        if y_train.nunique() < 2 or y_test.nunique() < 2:
            row["status"] = "skipped_single_class"
            fold_rows.append(row)
            continue
        model = _model(
            random_state + number,
            model_type=model_config["model_type"],
            C=model_config["C"],
            l1_ratio=model_config.get("l1_ratio"),
        )
        model, converged = _fit_checked(model, X.iloc[train_idx], y_train)
        if not converged:
            row["status"] = "skipped_nonconverged"
            fold_rows.append(row)
            continue
        prob = model.predict_proba(X.iloc[test_idx])[:, 1]
        pred = (prob >= 0.5).astype(int)
        oof[test_idx] = prob
        row.update({
            "status": "used",
            "roc_auc": _safe_metric(roc_auc_score, y_test, prob),
            "average_precision": _safe_metric(average_precision_score, y_test, prob),
            "balanced_accuracy": _safe_metric(balanced_accuracy_score, y_test, pred),
        })
        fold_rows.append(row)

    fold_df = pd.DataFrame(fold_rows)
    valid = np.isfinite(oof)
    pred_df = pd.DataFrame({
        "row_index": np.arange(len(y)),
        "dataset_id": dataset_ids.to_numpy(),
        "label": y.to_numpy(),
        "loso_probability": oof,
    })
    used = int((fold_df.get("status", pd.Series(dtype=str)) == "used").sum())
    if valid.sum() == 0:
        return ({
            "loso_dataset_cv_enabled": False,
            "loso_dataset_cv_reason": "no_two_class_held_out_dataset",
            "loso_dataset_cv_folds_used": used,
        }, pred_df, fold_df)

    y_valid = y.iloc[np.where(valid)[0]]
    p_valid = oof[valid]
    return ({
        "loso_dataset_cv_enabled": True,
        "loso_dataset_cv_folds_used": used,
        "loso_dataset_cv_oof_n": int(valid.sum()),
        "loso_dataset_cv_oof_roc_auc": _safe_metric(roc_auc_score, y_valid, p_valid),
        "loso_dataset_cv_oof_average_precision": _safe_metric(
            average_precision_score, y_valid, p_valid
        ),
        "loso_dataset_cv_oof_balanced_accuracy": _safe_metric(
            balanced_accuracy_score, y_valid, (p_valid >= 0.5).astype(int)
        ),
    }, pred_df, fold_df)

def _grouped_heldout(
    X: pd.DataFrame,
    y: pd.Series,
    groups: pd.Series,
    *,
    test_size: float,
    random_state: int,
    model_config: dict[str, Any],
) -> tuple[dict[str, Any], pd.DataFrame]:
    if test_size <= 0 or GroupShuffleSplit is None:
        return {"heldout_enabled": False, "heldout_reason": "disabled"}, pd.DataFrame()

    splitter = GroupShuffleSplit(
        n_splits=50, test_size=test_size, random_state=random_state
    )
    for attempt, (train_idx, test_idx) in enumerate(splitter.split(X, y, groups), start=1):
        y_train = y.iloc[train_idx]
        y_test = y.iloc[test_idx]
        if y_train.nunique() < 2 or y_test.nunique() < 2:
            continue
        model = _model(
            random_state,
            model_type=model_config["model_type"],
            C=model_config["C"],
            l1_ratio=model_config.get("l1_ratio"),
        )
        model, converged = _fit_checked(model, X.iloc[train_idx], y_train)
        if not converged:
            continue
        prob = model.predict_proba(X.iloc[test_idx])[:, 1]
        pred = (prob >= 0.5).astype(int)
        pred_df = pd.DataFrame({
            "row_index": test_idx,
            "group_id": groups.iloc[test_idx].to_numpy(),
            "label": y_test.to_numpy(),
            "heldout_probability": prob,
        })
        return ({
            "heldout_enabled": True,
            "heldout_method": "GroupShuffleSplit",
            "heldout_attempt": attempt,
            "heldout_test_size": test_size,
            "heldout_n_test": int(len(test_idx)),
            "heldout_n_test_groups": int(groups.iloc[test_idx].nunique()),
            "heldout_roc_auc": _safe_metric(roc_auc_score, y_test, prob),
            "heldout_average_precision": _safe_metric(average_precision_score, y_test, prob),
            "heldout_balanced_accuracy": _safe_metric(balanced_accuracy_score, y_test, pred),
        }, pred_df)

    return {
        "heldout_enabled": False,
        "heldout_reason": "could_not_create_two-class_grouped_split",
    }, pd.DataFrame()


def _evidence_tiebreaker(M: pd.DataFrame) -> pd.Series:
    cols = [c for c in SUPERVISED_FEATURES if c in M.columns]
    if not cols:
        return pd.Series(0.0, index=M.index)
    values = M[cols].apply(pd.to_numeric, errors="coerce")
    return values.mean(axis=1, skipna=True).fillna(0.0).clip(0.0, 1.0)


def _recompute_priority_reporting(M: pd.DataFrame) -> pd.DataFrame:
    """Rebuild tiers/display values after supervised probabilities replace score."""
    M = M.copy()
    gate = M.get("priority_gate_pass", pd.Series(False, index=M.index)).fillna(False).astype(bool)
    strict_scores = pd.to_numeric(M.loc[gate, "rank_score"], errors="coerce")
    q75 = strict_scores.quantile(0.75) if len(strict_scores) else np.inf
    relaxed = M.get("dsRNA_case_priority", pd.Series("", index=M.index)).eq(
        "case_supported_missing_RI_or_annotation"
    )
    M["priority_tier"] = np.select(
        [gate & M["rank_score"].ge(q75), gate, relaxed],
        ["tier1_strict_high", "tier2_strict", "tier3_relaxed"],
        default="not_prioritized",
    )

    gate_cols = [
        c for c in [
            "priority_gate_orientation", "priority_gate_annotation",
            "priority_gate_case_TE", "priority_gate_case_editing",
            "priority_gate_case_RI",
        ] if c in M.columns
    ]
    M["strict_gate_count_0_to_5"] = (
        M[gate_cols].fillna(False).astype(bool).sum(axis=1) if gate_cols else 0
    )
    M["evidence_score_raw"] = pd.to_numeric(M["rank_score"], errors="coerce").fillna(0.0)
    M["evidence_score_percentile"] = M["evidence_score_raw"].rank(
        method="average", pct=True
    ).fillna(0.0)
    M["evidence_tiebreaker_score"] = _evidence_tiebreaker(M)
    # Exploratory model-only rank can reveal strong gate-failing candidates,
    # while the primary priority_rank below remains conservative/tier-aware.
    M["supervised_global_rank"] = M["rank_score"].rank(
        method="first", ascending=False
    ).astype(int)
    M["_tier_order"] = M["priority_tier"].map(TIER_ORDER).fillna(99).astype(int)

    M["within_tier_rank"] = 0
    M["within_tier_percentile"] = 0.0
    for tier, idx in M.groupby("priority_tier", sort=False).groups.items():
        ordered = M.loc[list(idx)].sort_values(
            ["rank_score", "evidence_tiebreaker_score", "pair_id"],
            ascending=[False, False, True],
            kind="mergesort",
        )
        n = len(ordered)
        ranks = pd.Series(np.arange(1, n + 1), index=ordered.index)
        pct = pd.Series(1.0 if n == 1 else (n - ranks) / (n - 1), index=ordered.index)
        M.loc[ordered.index, "within_tier_rank"] = ranks.astype(int)
        M.loc[ordered.index, "within_tier_percentile"] = pct.astype(float)

    band_map = {
        "tier1_strict_high": (90.0, 100.0, "90-100"),
        "tier2_strict": (70.0, 89.0, "70-89"),
        "tier3_relaxed": (50.0, 69.0, "50-69"),
        "not_prioritized": (0.0, 49.0, "0-49"),
    }
    display: list[float] = []
    bands: list[str] = []
    for tier, pct in zip(M["priority_tier"].astype(str), M["within_tier_percentile"]):
        lo, hi, label = band_map.get(tier, (0.0, 49.0, "0-49"))
        display.append(lo + float(pct) * (hi - lo))
        bands.append(label)
    M["display_priority_score"] = pd.Series(display, index=M.index).round(2)
    M["tier_display_band"] = bands

    M = M.sort_values(
        ["_tier_order", "rank_score", "evidence_tiebreaker_score", "pair_id"],
        ascending=[True, False, False, True],
        kind="mergesort",
    ).drop(columns="_tier_order")
    M["priority_rank"] = np.arange(1, len(M) + 1)
    return M


def _prepare_training_data(
    M: pd.DataFrame,
    labels_df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, bool]:
    """Return training rows, target rows, and whether training is external."""
    feature_count = sum(c in labels_df.columns for c in [*SUPERVISED_FEATURES, *ROUTE6_FEATURES])
    external = feature_count >= 2

    if external:
        training = labels_df.copy()
        target = M.copy()
        target["supervised_training_label"] = pd.Series(pd.NA, index=target.index, dtype="Float64")
        target["supervised_label_source"] = "external_training_matrix"
        return training, target, True

    if "pair_id" not in labels_df.columns:
        raise ValueError("Target-specific labels require pair_id")

    # Reject contradictory duplicate labels before merge.
    labeled = labels_df.loc[labels_df["label"].notna(), ["pair_id", "label"]]
    conflicting = labeled.groupby("pair_id")["label"].nunique()
    bad = conflicting[conflicting > 1]
    if len(bad):
        raise ValueError(f"Conflicting labels for pair IDs: {list(bad.index[:10])}")
    labels_one = labels_df.drop_duplicates("pair_id", keep="first").copy()

    merge_cols = [
        c for c in [
            "pair_id", "label", "label_source", "matched_positive_symbols",
            "matched_negative_symbols", "matched_truth_symbols",
        ] if c in labels_one.columns
    ]
    target = M.merge(labels_one[merge_cols], on="pair_id", how="left")
    target["supervised_training_label"] = _numeric_label(target["label"])
    target["supervised_label_source"] = target.get(
        "label_source", pd.Series("", index=target.index)
    ).fillna("unlabeled")
    target = target.drop(columns=["label", "label_source"], errors="ignore")
    training = target.loc[target["supervised_training_label"].notna()].copy()
    training["label"] = training["supervised_training_label"].astype(int)
    return training, target, False


def apply_supervised_priority(
    M: pd.DataFrame,
    args,
    *,
    outdir: str | Path,
    case: str,
    control: str,
) -> pd.DataFrame:
    """Fit optional label-dependent reranking and report unbiased evaluation.

    The final model is a deployment model fitted to all labeled training rows.
    When tuning is enabled, performance is estimated separately by true nested
    grouped cross-validation. Fitted-all-label scores are never reported as an
    independent performance estimate.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    M = ensure_route_features(M.copy())
    training_labels = getattr(args, "training_labels", None)
    truth_table = getattr(args, "training_truth_table", None)
    if training_labels:
        labels_df, truth_full, truth_positive = read_pair_labels(training_labels, M)
    elif truth_table:
        labels_df, truth_full, truth_positive = derive_pair_labels_from_truth_table(
            M, truth_table,
            symbol_col=getattr(args, "truth_symbol_col", "Symbol"),
            truth_label_mode=getattr(args, "truth_label_mode", "positive_logfc_padj"),
            truth_label_col=getattr(args, "truth_label_col", None),
            padj_col=getattr(args, "truth_padj_col", "padj"),
            logfc_col=getattr(args, "truth_logfc_col", "log2FoldChange"),
            padj_max=float(getattr(args, "truth_padj_max", 0.05)),
        )
    else:
        raise ValueError(
            "--priority-score-mode supervised requires --training-labels or --training-truth-table"
        )
    labels_df = ensure_route_features(labels_df)
    labels_df.to_csv(
        outdir / f"TEpair_dsRNA_supervised_labels_input.{case}.tsv", sep="\t", index=False
    )
    if truth_full is not None:
        truth_full.to_csv(outdir / f"TEpair_dsRNA_supervised_truth_genes_full_labeled.{case}.csv", index=False)
    if truth_positive is not None:
        truth_positive.to_csv(outdir / f"TEpair_dsRNA_supervised_truth_genes_positive.{case}.csv", index=False)

    training, target, external = _prepare_training_data(M, labels_df)
    training = ensure_route_features(training)
    target = ensure_route_features(target)
    training["label"] = _numeric_label(training["label"])
    training = training.loc[training["label"].notna()].copy().reset_index(drop=True)
    training["label"] = training["label"].astype(int)
    y = training["label"].astype(int).reset_index(drop=True)
    groups = _derive_group_id(training).reset_index(drop=True)
    n_pos, n_neg = int(y.eq(1).sum()), int(y.eq(0).sum())
    if n_pos < 2 or n_neg < 2:
        raise ValueError(
            "Supervised training requires at least 2 positive and 2 confident-negative labeled rows; "
            f"found positives={n_pos}, negatives={n_neg}. Do not turn unlabeled candidates into negatives."
        )
    random_state = int(getattr(args, "supervised_random_state", 1))
    tune = bool(getattr(args, "supervised_tune", False)) or str(
        getattr(args, "supervised_feature_panel", "raw7")
    ) == "auto_prespecified"
    if tune:
        nested_report, oof_df, fold_df = _nested_grouped_cv(
            training, target, y, groups, args
        )
    else:
        nested_report, oof_df, fold_df = ({
            "nested_cv_enabled": False,
            "nested_cv_reason": "fixed_configuration_requested",
        }, pd.DataFrame(), pd.DataFrame())

    family_col = None
    for col in ("study_family", "dataset_family"):
        if col in training.columns:
            family_col = training[col].fillna("").astype(str).reset_index(drop=True)
            break
    config, tuning_df = _select_configuration(
        training, target, y, groups, args, family=family_col
    )
    weights = (
        _family_class_balanced_weights(y, family_col)
        if family_col is not None else _class_balanced_weights(y)
    )
    final_model, config, fit_status = _fit_selected_with_fallback(
        training, target, y, config, tuning_df,
        random_state=random_state, sample_weight=weights,
    )
    features = list(config["features"])
    if not features:
        raise ValueError("Selected supervised configuration contains no usable features")
    score = final_model.predict_proba(
        target[features].apply(pd.to_numeric, errors="coerce")
    )[:, 1]
    # Keep the historical probability column as a compatibility alias, but the
    # scientifically preferred term is ranking score unless calibration is
    # independently demonstrated.
    target["supervised_priority_probability"] = score
    target["supervised_priority_score"] = score
    target["supervised_model_type"] = config["model_type"]
    target["supervised_feature_panel"] = config["feature_panel"]
    target["case_priority_score"] = score
    target["rank_score"] = score

    report: dict[str, Any] = {
        "model": "class_balanced_regularized_logistic_regression",
        "model_type": config["model_type"],
        "feature_panel": config["feature_panel"],
        "selection_status": config["selection_status"],
        "final_fit_status": fit_status,
        "selection_metric": config.get("selection_metric"),
        "inner_cv_method": config.get("inner_cv_method"),
        "inner_cv_folds_requested": config.get("inner_cv_folds_requested"),
        "inner_cv_score": config.get("inner_cv_score"),
        "tuning_enabled": tune, "C": config["C"],
        "l1_ratio": config.get("l1_ratio"),
        "training_source": "external_feature_matrix" if external else "current_candidate_labels",
        "ranking_application": "external_to_target" if external else "same_candidate_universe",
        "n_target_candidates": int(len(target)), "n_training_rows": int(len(training)),
        "n_positive_labels": n_pos, "n_confident_negative_labels": n_neg,
        "n_training_groups": int(groups.nunique()), "features_used": ";".join(features),
        "unlabeled_as_negative": False,
        "group_definition": "provided group_id, otherwise gene pair, otherwise 100-kb locus",
        "sample_weighting": "study-and-class-balanced" if family_col is not None else "class-balanced",
        "interpretation": "ranking score; not a calibrated probability without independent calibration",
        "performance_source": "nested OOF predictions when nested_cv_enabled; never fitted-all-labels scores",
        **nested_report,
    }
    lr = final_model.named_steps["logistic"]
    pd.DataFrame({
        "feature": features, "coefficient_standardized": lr.coef_[0],
        "odds_ratio_per_1SD": np.exp(lr.coef_[0]),
        "model_type": config["model_type"], "feature_panel": config["feature_panel"],
        "C": config["C"], "l1_ratio": config.get("l1_ratio"),
    }).sort_values("coefficient_standardized", ascending=False).to_csv(
        outdir / f"TEpair_dsRNA_supervised_coefficients.{case}.csv", index=False
    )
    tuning_df.to_csv(
        outdir / f"TEpair_dsRNA_supervised_final_model_selection.{case}.csv", index=False
    )
    train_ids = pd.DataFrame({
        "training_row": np.arange(len(training)),
        "pair_id": training.get("pair_id", pd.Series("", index=training.index)).astype(str).to_numpy(),
        "dataset_id": training.get("dataset_id", pd.Series("", index=training.index)).astype(str).to_numpy(),
        "study_family": training.get("study_family", pd.Series("", index=training.index)).astype(str).to_numpy(),
        "group_id": groups.to_numpy(), "label": y.to_numpy(),
    })
    train_ids.to_csv(
        outdir / f"TEpair_dsRNA_supervised_training_rows.{case}.tsv", sep="\t", index=False
    )
    if not oof_df.empty:
        oof_df.merge(train_ids, left_on="row_index", right_on="training_row", how="left").to_csv(
            outdir / f"TEpair_dsRNA_supervised_nested_oof_predictions.{case}.tsv",
            sep="\t", index=False,
        )
    if not fold_df.empty:
        fold_df.to_csv(outdir / f"TEpair_dsRNA_supervised_nested_outer_folds.{case}.csv", index=False)
    pd.DataFrame([report]).to_csv(
        outdir / f"TEpair_dsRNA_supervised_training_report.{case}.csv", index=False
    )
    with open(outdir / f"TEpair_dsRNA_supervised_training_report.{case}.json", "w") as handle:
        json.dump(report, handle, indent=2, default=str)
    try:
        import joblib
        joblib.dump({
            "model": final_model, "features": features, "configuration": config,
            "report": report, "route_definitions": ROUTE_DEFINITIONS,
        }, outdir / f"TEpair_dsRNA_supervised_model.{case}.joblib")
    except Exception as exc:  # pragma: no cover
        report["model_save_warning"] = str(exc)
    return _recompute_priority_reporting(target)

