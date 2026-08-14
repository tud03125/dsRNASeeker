from __future__ import annotations

"""
Label-free ranking-robustness diagnostics for dsRNASeeker.

This module does NOT use RIP/J2/Z22/dsRNA-seq labels and does NOT replace ADPS.
It asks a narrower question that can be answered from an unlabeled candidate
table: how stable is each candidate's rank under several pre-specified,
biologically interpretable evidence views?

The broad profile set mirrors the fixed profiles used in the manuscript
benchmark:
    - integrated ADPS
    - all seven components, equal weight
    - duplex structure only (energy + interface)
    - geometry + annotation
    - condition evidence (expression + editing + RI)

The six route scores mirror modules.supervised.ROUTE_DEFINITIONS.  They are
reported separately and are not allowed to vote six times in the broad-profile
agreement diagnostic.

Agreement is not accuracy.  A candidate can rank consistently across profiles
and still fail experimental validation; conversely, a mechanism-specific true
positive can be profile-sensitive.  The outputs are decision-support
diagnostics, not calibrated probabilities of duplex formation.
"""

from pathlib import Path
from typing import Iterable
import json
import math
import re

import numpy as np
import pandas as pd

from .supervised import ROUTE_DEFINITIONS


COMPONENTS = {
    "orientation": "orientation_adps",
    "annotation": "annotation_adps",
    "expression": "case_expression_adps",
    "energy": "energy_adps",
    "interface": "interface_adps",
    "editing": "case_editing_adps",
    "RI": "RI_adps",
}

WEIGHT_COLUMNS = {
    name: f"adaptive_weight_{col}"
    for name, col in COMPONENTS.items()
}

BROAD_PROFILE_DEFINITIONS = {
    # integrated_adps is read from the stored score, not recomputed here.
    "equal_weight_all": list(COMPONENTS),
    "duplex_structure_only": ["energy", "interface"],
    "geometry_annotation": ["orientation", "annotation"],
    "condition_evidence": ["expression", "editing", "RI"],
}

# Diagnostic cutoffs only.  They describe rank agreement and weight
# concentration; they are NOT accuracy thresholds.
DEFAULT_DOMINANT_WEIGHT_THRESHOLD = 0.80
DEFAULT_TIGHT_IQR = 0.10
DEFAULT_WIDE_IQR = 0.25

# Human-readable diagnostic cutoffs. These are operational interpretation
# thresholds, not accuracy/calibration thresholds and were not fit to RIP labels.
DEFAULT_MIN_REFERENCE_N = 25
DEFAULT_MIN_BACKGROUND_N = 25
DEFAULT_MIN_REFERENCE_FRACTION = 0.05
DEFAULT_MAX_REFERENCE_FRACTION = 0.95

# ADPS numerical-resolution cutoffs. A score is called COARSE when it has very
# few distinct values, a very low unique-score fraction, or a very large tied
# group. HIGH/MODERATE/COARSE describe ranking granularity only.
DEFAULT_COARSE_UNIQUE_N = 20
DEFAULT_MODERATE_UNIQUE_N = 100
DEFAULT_COARSE_UNIQUE_FRACTION = 0.01
DEFAULT_MODERATE_UNIQUE_FRACTION = 0.05
DEFAULT_COARSE_LARGEST_TIE_FRACTION = 0.25
DEFAULT_MODERATE_LARGEST_TIE_FRACTION = 0.10


def _numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def _component_table(df: pd.DataFrame) -> pd.DataFrame:
    out = pd.DataFrame(index=df.index)
    for name, col in COMPONENTS.items():
        if col in df.columns:
            out[name] = _numeric(df[col]).clip(0.0, 1.0)
        else:
            out[name] = np.nan
    return out


def _mean_score(comp: pd.DataFrame, keep: Iterable[str]) -> pd.Series:
    cols = [x for x in keep if x in comp.columns]
    if not cols:
        return pd.Series(np.nan, index=comp.index, dtype=float)
    return comp[cols].mean(axis=1, skipna=True)


def _stored_adps(df: pd.DataFrame) -> pd.Series:
    for col in ("adaptive_priority_score", "case_priority_score", "rank_score"):
        if col in df.columns:
            x = _numeric(df[col])
            if x.notna().any():
                return x
    raise ValueError(
        "No stored ADPS score found. Expected one of: "
        "adaptive_priority_score, case_priority_score, rank_score."
    )


def _rank_percentile(score: pd.Series) -> pd.Series:
    """Higher score -> higher percentile; ties receive their average percentile."""
    x = _numeric(score)
    out = pd.Series(np.nan, index=x.index, dtype=float)
    ok = x.notna()
    if ok.any():
        out.loc[ok] = x.loc[ok].rank(method="average", pct=True)
    return out


def build_profile_scores(df: pd.DataFrame) -> pd.DataFrame:
    """Return fixed, label-free profile scores for every candidate."""
    comp = _component_table(df)
    scores = pd.DataFrame(index=df.index)
    scores["integrated_adps"] = _stored_adps(df)

    for name, parts in BROAD_PROFILE_DEFINITIONS.items():
        scores[name] = _mean_score(comp, parts)

    for route_name, feature_cols in ROUTE_DEFINITIONS.items():
        # ROUTE_DEFINITIONS uses keys such as score_route_substrate_context.
        # Store the internal profile name without the leading "score_" because
        # candidate output columns add that prefix consistently below.
        clean_route_name = (
            route_name[len("score_"):]
            if str(route_name).startswith("score_")
            else str(route_name)
        )
        parts = []
        for feature_col in feature_cols:
            matching = [k for k, v in COMPONENTS.items() if v == feature_col]
            if matching:
                parts.append(matching[0])
        scores[clean_route_name] = _mean_score(comp, parts)

    route_cols = [c for c in scores.columns if c.startswith("route_")]
    if route_cols:
        scores["route_any_max"] = scores[route_cols].max(axis=1, skipna=True)

    return scores


def _weights_from_summary(df: pd.DataFrame) -> tuple[dict[str, float], str]:
    vals: dict[str, float] = {}
    for name, col in WEIGHT_COLUMNS.items():
        if col in df.columns:
            x = _numeric(df[col]).dropna()
            vals[name] = float(x.median()) if len(x) else 0.0
        else:
            vals[name] = 0.0

    source = ""
    if "adaptive_weight_source" in df.columns:
        x = df["adaptive_weight_source"].dropna().astype(str)
        if len(x):
            source = str(x.iloc[0])

    positive = {k: max(0.0, float(v)) for k, v in vals.items()}
    total = sum(positive.values())
    if total <= 0:
        positive = {k: 1.0 / len(COMPONENTS) for k in COMPONENTS}
        if not source:
            source = "uniform_fallback"
    else:
        positive = {k: v / total for k, v in positive.items()}

    return positive, source


def _first_numeric_value(df: pd.DataFrame, col: str, default: float = np.nan) -> float:
    if col not in df.columns:
        return default
    x = _numeric(df[col]).dropna()
    return float(x.iloc[0]) if len(x) else default


def _weight_diagnostics(
    df: pd.DataFrame,
    dominant_weight_threshold: float = DEFAULT_DOMINANT_WEIGHT_THRESHOLD,
) -> dict[str, object]:
    weights, source = _weights_from_summary(df)
    w = np.asarray(list(weights.values()), dtype=float)
    w = np.where(w > 0, w, 0.0)
    if w.sum() > 0:
        w = w / w.sum()

    max_i = int(np.argmax(w)) if len(w) else 0
    names = list(weights)
    max_weight = float(w[max_i]) if len(w) else np.nan
    dominant_block = names[max_i] if len(names) else ""

    effective_n = float(1.0 / np.sum(w * w)) if len(w) and np.sum(w * w) > 0 else np.nan
    nz = w[w > 0]
    if len(nz) <= 1:
        entropy_norm = 0.0
    else:
        entropy = float(-np.sum(nz * np.log(nz)))
        entropy_norm = float(entropy / math.log(len(w))) if len(w) > 1 else 0.0

    ref_n = _first_numeric_value(df, "adaptive_gate_positive_n", np.nan)
    bg_n = _first_numeric_value(df, "adaptive_gate_background_n", np.nan)
    denom = ref_n + bg_n if np.isfinite(ref_n) and np.isfinite(bg_n) else np.nan
    ref_fraction = float(ref_n / denom) if np.isfinite(denom) and denom > 0 else np.nan

    source_low = source.lower()
    if "uniform" in source_low or (np.isfinite(ref_n) and ref_n == 0):
        mode = "equal_weight_fallback"
    elif np.isfinite(max_weight) and max_weight >= float(dominant_weight_threshold):
        mode = "single_component_dominated"
    else:
        mode = "dataset_adaptive"

    out: dict[str, object] = {
        "adps_weighting_mode": mode,
        "adaptive_weight_source": source,
        "adaptive_reference_n": int(ref_n) if np.isfinite(ref_n) else np.nan,
        "adaptive_background_n": int(bg_n) if np.isfinite(bg_n) else np.nan,
        "adaptive_reference_fraction": ref_fraction,
        "adps_max_weight": max_weight,
        "adps_dominant_component": dominant_block,
        "adps_effective_component_n": effective_n,
        "adps_nonzero_weight_n": int(np.sum(w > 1e-12)),
        "adps_weight_entropy_normalized": entropy_norm,
        "dominant_weight_threshold": float(dominant_weight_threshold),
    }
    for name, value in zip(names, w):
        out[f"adps_weight_{name}"] = float(value)
    return out


def _score_granularity(score: pd.Series) -> dict[str, float | int]:
    x = _numeric(score).dropna()
    if len(x) == 0:
        return {
            "adps_scored_candidate_n": 0,
            "adps_unique_score_n": 0,
            "adps_unique_score_fraction": np.nan,
            "adps_largest_tie_n": 0,
            "adps_largest_tie_fraction": np.nan,
        }
    counts = x.value_counts(dropna=True)
    return {
        "adps_scored_candidate_n": int(len(x)),
        "adps_unique_score_n": int(x.nunique(dropna=True)),
        "adps_unique_score_fraction": float(x.nunique(dropna=True) / len(x)),
        "adps_largest_tie_n": int(counts.max()),
        "adps_largest_tie_fraction": float(counts.max() / len(x)),
    }


def _adps_weight_status(
    wdiag: dict[str, object],
    *,
    min_reference_n: int = DEFAULT_MIN_REFERENCE_N,
    min_background_n: int = DEFAULT_MIN_BACKGROUND_N,
    min_reference_fraction: float = DEFAULT_MIN_REFERENCE_FRACTION,
    max_reference_fraction: float = DEFAULT_MAX_REFERENCE_FRACTION,
) -> str:
    """Return a descriptive ADPS provenance/support label.

    This status is intentionally label-free. ADAPTIVE_SUPPORTED means only that
    the non-uniform weight solution was estimated from non-trivial internal
    reference and background groups. It does not mean experimentally validated.
    """
    mode = str(wdiag.get("adps_weighting_mode", ""))
    if mode == "equal_weight_fallback":
        return "EQUAL_WEIGHT_FALLBACK"
    if mode == "single_component_dominated":
        return "SINGLE_COMPONENT_DOMINATED"

    ref_n = pd.to_numeric(
        pd.Series([wdiag.get("adaptive_reference_n")]), errors="coerce"
    ).iloc[0]
    bg_n = pd.to_numeric(
        pd.Series([wdiag.get("adaptive_background_n")]), errors="coerce"
    ).iloc[0]
    frac = pd.to_numeric(
        pd.Series([wdiag.get("adaptive_reference_fraction")]), errors="coerce"
    ).iloc[0]

    supported = (
        np.isfinite(ref_n) and ref_n >= int(min_reference_n) and
        np.isfinite(bg_n) and bg_n >= int(min_background_n) and
        np.isfinite(frac) and
        float(min_reference_fraction) <= frac <= float(max_reference_fraction)
    )
    return "ADAPTIVE_SUPPORTED" if supported else "LIMITED_REFERENCE"


def _adps_score_resolution(granularity: dict[str, float | int]) -> str:
    """Classify ADPS score granularity; never interpret this as accuracy."""
    n_unique = pd.to_numeric(
        pd.Series([granularity.get("adps_unique_score_n")]), errors="coerce"
    ).iloc[0]
    unique_frac = pd.to_numeric(
        pd.Series([granularity.get("adps_unique_score_fraction")]), errors="coerce"
    ).iloc[0]
    tie_frac = pd.to_numeric(
        pd.Series([granularity.get("adps_largest_tie_fraction")]), errors="coerce"
    ).iloc[0]

    if (
        not np.isfinite(n_unique) or n_unique < DEFAULT_COARSE_UNIQUE_N or
        not np.isfinite(unique_frac) or unique_frac < DEFAULT_COARSE_UNIQUE_FRACTION or
        (np.isfinite(tie_frac) and tie_frac >= DEFAULT_COARSE_LARGEST_TIE_FRACTION)
    ):
        return "COARSE"
    if (
        n_unique < DEFAULT_MODERATE_UNIQUE_N or
        unique_frac < DEFAULT_MODERATE_UNIQUE_FRACTION or
        (np.isfinite(tie_frac) and tie_frac >= DEFAULT_MODERATE_LARGEST_TIE_FRACTION)
    ):
        return "MODERATE"
    return "HIGH"


def _candidate_user_interpretation(
    candidate: pd.DataFrame,
    *,
    adps_weight_status: str,
    tight_iqr: float,
    wide_iqr: float,
) -> pd.Series:
    """Return a non-probabilistic, user-facing candidate interpretation."""
    iqr = pd.to_numeric(
        candidate["broad_profile_percentile_iqr"], errors="coerce"
    )
    top10 = pd.to_numeric(
        candidate["broad_profile_top10_count"], errors="coerce"
    ).fillna(0)

    # Strongest candidate-level statement we can make without labels: high rank
    # under multiple broad views with tight dispersion. This remains meaningful
    # even when ADPS used a fallback because it describes cross-profile stability.
    robust = (top10 >= 3) & (iqr <= float(tight_iqr))
    out = pd.Series("", index=candidate.index, dtype=object)
    out.loc[robust] = "ROBUST_MULTI_PROFILE_PRIORITY"

    remaining = ~robust
    if adps_weight_status == "SINGLE_COMPONENT_DOMINATED":
        out.loc[remaining] = "SINGLE_COMPONENT_DRIVEN"
    elif adps_weight_status == "EQUAL_WEIGHT_FALLBACK":
        out.loc[remaining] = "FALLBACK_RANKING"
    elif adps_weight_status == "LIMITED_REFERENCE":
        out.loc[remaining] = "LIMITED_REFERENCE_RANKING"
    else:
        wide = iqr > float(wide_iqr)
        out.loc[remaining & wide] = "ADPS_SUPPORTED_PROFILE_SENSITIVE"
        out.loc[remaining & ~wide] = "ADPS_SUPPORTED_GENERAL"
    return out


def _dataset_guidance_label(weight_status: str, score_resolution: str) -> str:
    if weight_status == "EQUAL_WEIGHT_FALLBACK":
        return "INSPECT_MULTIPLE_PROFILES_FALLBACK"
    if weight_status == "SINGLE_COMPONENT_DOMINATED":
        return "INSPECT_MULTIPLE_PROFILES_SINGLE_COMPONENT"
    if weight_status == "LIMITED_REFERENCE":
        return "INSPECT_MULTIPLE_PROFILES_LIMITED_REFERENCE"
    if score_resolution == "COARSE":
        return "ADPS_REFERENCE_SUPPORTED_BUT_COARSE"
    return "ADPS_REFERENCE_SUPPORTED"


def _pairwise_spearman(scores: pd.DataFrame) -> pd.DataFrame:
    corr = scores.corr(method="spearman", min_periods=3)
    rows = []
    cols = list(corr.columns)
    for i, a in enumerate(cols):
        for j, b in enumerate(cols):
            if j <= i:
                continue
            rows.append({
                "profile_a": a,
                "profile_b": b,
                "spearman_rho": corr.loc[a, b],
            })
    return pd.DataFrame(rows, columns=["profile_a", "profile_b", "spearman_rho"])


def add_ranking_robustness(
    df: pd.DataFrame,
    *,
    dominant_weight_threshold: float = DEFAULT_DOMINANT_WEIGHT_THRESHOLD,
    tight_iqr: float = DEFAULT_TIGHT_IQR,
    wide_iqr: float = DEFAULT_WIDE_IQR,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Build candidate-, dataset-, and profile-level label-free diagnostics.

    Returns:
        candidate_table, dataset_summary, profile_correlations
    """
    if "pair_id" not in df.columns:
        raise ValueError("Summary must contain pair_id.")

    scores = build_profile_scores(df)

    # Five broad views only.  The six biological routes are reported separately
    # so overlapping route definitions do not overweight particular components.
    broad = [
        "integrated_adps",
        "equal_weight_all",
        "duplex_structure_only",
        "geometry_annotation",
        "condition_evidence",
    ]
    broad = [c for c in broad if c in scores.columns]

    pct = pd.DataFrame(index=df.index)
    for col in scores.columns:
        pct[f"pct_{col}"] = _rank_percentile(scores[col])

    broad_pct_cols = [f"pct_{c}" for c in broad]
    B = pct[broad_pct_cols]

    candidate = pd.DataFrame(index=df.index)
    keep_identity = [
        "pair_id", "priority_tier", "priority_gate_pass", "dsRNA_case_priority",
        "A_TE_id", "B_TE_id", "A_SYMBOL", "B_SYMBOL",
        "A_repFamily", "A_repName", "B_repFamily", "B_repName",
        "A_annotation", "B_annotation",
    ]
    for col in keep_identity:
        if col in df.columns:
            candidate[col] = df[col]

    for col in scores.columns:
        candidate[f"score_{col}"] = scores[col]
    candidate = pd.concat([candidate, pct], axis=1)

    candidate["broad_profile_percentile_median"] = B.median(axis=1, skipna=True)
    candidate["broad_profile_percentile_q25"] = B.quantile(0.25, axis=1)
    candidate["broad_profile_percentile_q75"] = B.quantile(0.75, axis=1)
    candidate["broad_profile_percentile_iqr"] = (
        candidate["broad_profile_percentile_q75"] -
        candidate["broad_profile_percentile_q25"]
    )
    candidate["broad_profile_percentile_min"] = B.min(axis=1, skipna=True)
    candidate["broad_profile_percentile_max"] = B.max(axis=1, skipna=True)
    candidate["broad_profile_percentile_range"] = (
        candidate["broad_profile_percentile_max"] -
        candidate["broad_profile_percentile_min"]
    )
    candidate["broad_profile_top10_count"] = (B >= 0.90).sum(axis=1)
    candidate["broad_profile_top25_count"] = (B >= 0.75).sum(axis=1)
    candidate["broad_profile_bottom25_count"] = (B <= 0.25).sum(axis=1)
    candidate["adps_vs_broad_median_abs_gap"] = (
        candidate["pct_integrated_adps"] -
        candidate["broad_profile_percentile_median"]
    ).abs()

    # Descriptive agreement bands only; these are not experimentally calibrated
    # confidence categories.
    iqr = candidate["broad_profile_percentile_iqr"]
    candidate["broad_profile_agreement_band"] = np.select(
        [iqr <= float(tight_iqr), iqr > float(wide_iqr)],
        ["tight", "wide"],
        default="moderate",
    )
    candidate["broad_profile_agreement_note"] = np.select(
        [
            (candidate["broad_profile_top10_count"] >= 3) & (iqr <= float(tight_iqr)),
            (candidate["broad_profile_top10_count"] >= 1) & (iqr > float(wide_iqr)),
        ],
        [
            "high_rank_across_multiple_broad_profiles",
            "high_in_at_least_one_profile_but_profile_sensitive",
        ],
        default="descriptive_only",
    )

    # Route diagnostics are kept separate from broad-profile agreement.
    route_pct_cols = [
        c for c in pct.columns
        if c.startswith("pct_route_") and c != "pct_route_any_max"
    ]
    if route_pct_cols:
        R = pct[route_pct_cols]
        candidate["route_percentile_median"] = R.median(axis=1, skipna=True)
        candidate["route_percentile_iqr"] = R.quantile(0.75, axis=1) - R.quantile(0.25, axis=1)
        candidate["route_top10_count"] = (R >= 0.90).sum(axis=1)
        candidate["route_top25_count"] = (R >= 0.75).sum(axis=1)

    # Attach the same dataset-level ADPS provenance to every candidate so the
    # table is self-describing when subsetted/exported.
    wdiag = _weight_diagnostics(df, dominant_weight_threshold)
    for key in (
        "adps_weighting_mode", "adaptive_reference_n",
        "adaptive_reference_fraction", "adps_max_weight",
        "adps_dominant_component", "adps_effective_component_n",
    ):
        candidate[key] = wdiag.get(key)

    # Human-readable status fields. These preserve all numerical diagnostics and
    # add plain-language labels; they are not probabilities or validation calls.
    adps_weight_status = _adps_weight_status(wdiag)
    granularity = _score_granularity(scores["integrated_adps"])
    adps_score_resolution = _adps_score_resolution(granularity)
    candidate["ADPS_WEIGHT_STATUS"] = adps_weight_status
    candidate["ADPS_SCORE_RESOLUTION"] = adps_score_resolution
    candidate["CANDIDATE_RANK_STABILITY"] = (
        candidate["broad_profile_agreement_band"].astype(str).str.upper()
    )
    candidate["USER_INTERPRETATION"] = _candidate_user_interpretation(
        candidate,
        adps_weight_status=adps_weight_status,
        tight_iqr=tight_iqr,
        wide_iqr=wide_iqr,
    )

    # Dataset-level summary.
    summary_row: dict[str, object] = {
        "candidate_n": int(len(df)),
        "broad_profile_n": int(len(broad)),
        "tight_iqr_threshold": float(tight_iqr),
        "wide_iqr_threshold": float(wide_iqr),
        "candidate_tight_agreement_fraction": float((iqr <= float(tight_iqr)).mean()) if len(iqr) else np.nan,
        "candidate_wide_agreement_fraction": float((iqr > float(wide_iqr)).mean()) if len(iqr) else np.nan,
        "candidate_high_multi_profile_fraction": float(
            ((candidate["broad_profile_top10_count"] >= 3) & (iqr <= float(tight_iqr))).mean()
        ) if len(candidate) else np.nan,
    }
    summary_row.update(wdiag)
    summary_row.update(granularity)
    summary_row["ADPS_WEIGHT_STATUS"] = adps_weight_status
    summary_row["ADPS_SCORE_RESOLUTION"] = adps_score_resolution
    summary_row["DATASET_USER_GUIDANCE"] = _dataset_guidance_label(
        adps_weight_status, adps_score_resolution
    )
    summary_row["robust_multi_profile_priority_fraction"] = float(
        candidate["USER_INTERPRETATION"].eq("ROBUST_MULTI_PROFILE_PRIORITY").mean()
    ) if len(candidate) else np.nan
    summary_row["profile_sensitive_fraction"] = float(
        candidate["CANDIDATE_RANK_STABILITY"].eq("WIDE").mean()
    ) if len(candidate) else np.nan

    broad_corr = _pairwise_spearman(scores[broad])
    finite_corr = pd.to_numeric(broad_corr.get("spearman_rho"), errors="coerce").dropna()
    summary_row["broad_profile_pairwise_spearman_median"] = float(finite_corr.median()) if len(finite_corr) else np.nan
    summary_row["broad_profile_pairwise_spearman_min"] = float(finite_corr.min()) if len(finite_corr) else np.nan

    adps_corr = broad_corr[
        broad_corr["profile_a"].eq("integrated_adps") |
        broad_corr["profile_b"].eq("integrated_adps")
    ] if len(broad_corr) else broad_corr
    finite_adps_corr = pd.to_numeric(adps_corr.get("spearman_rho"), errors="coerce").dropna()
    summary_row["adps_to_other_broad_profiles_spearman_median"] = float(finite_adps_corr.median()) if len(finite_adps_corr) else np.nan
    summary_row["adps_to_other_broad_profiles_spearman_min"] = float(finite_adps_corr.min()) if len(finite_adps_corr) else np.nan

    dataset_summary = pd.DataFrame([summary_row])

    all_corr = _pairwise_spearman(scores)
    return candidate, dataset_summary, all_corr


def _infer_case_label(summary_path: Path) -> str:
    parent = summary_path.parent
    patterns = [
        "TEpair_dsRNA_adaptive_weights_long.*.csv",
        "TEpair_dsRNA_adaptive_weights.*.csv",
    ]
    labels = []
    rx = re.compile(r"TEpair_dsRNA_adaptive_weights(?:_long)?\.(.+)\.csv$")
    for pattern in patterns:
        for p in parent.glob(pattern):
            m = rx.match(p.name)
            if m:
                labels.append(m.group(1))
    labels = sorted(set(labels))
    if len(labels) == 1:
        return labels[0]
    return "CASE"


def _guidance_markdown(
    dataset_summary: pd.DataFrame,
    *,
    case_label: str,
) -> str:
    r = dataset_summary.iloc[0].to_dict()
    mode = str(r.get("adps_weighting_mode", ""))
    ref_n = r.get("adaptive_reference_n", np.nan)
    ref_frac = r.get("adaptive_reference_fraction", np.nan)
    max_w = r.get("adps_max_weight", np.nan)
    dom = str(r.get("adps_dominant_component", ""))
    neff = r.get("adps_effective_component_n", np.nan)
    uniq = r.get("adps_unique_score_n", np.nan)
    n = r.get("adps_scored_candidate_n", np.nan)
    corr = r.get("adps_to_other_broad_profiles_spearman_median", np.nan)
    weight_status = str(r.get("ADPS_WEIGHT_STATUS", ""))
    score_resolution = str(r.get("ADPS_SCORE_RESOLUTION", ""))
    dataset_guidance = str(r.get("DATASET_USER_GUIDANCE", ""))

    def fmt(x, digits=3):
        try:
            if pd.isna(x):
                return "NA"
            return f"{float(x):.{digits}f}"
        except Exception:
            return str(x)

    mode_explanation = {
        "equal_weight_fallback":
            "The dataset-specific median-separation rule could not produce a positive weight solution; ADPS is an equal-weight fallback rather than a dataset-specific adaptive solution.",
        "single_component_dominated":
            f"One evidence block ({dom}) carries at least the configured dominant-weight threshold; ADPS is formally adaptive but is effectively dominated by that block.",
        "dataset_adaptive":
            "ADPS used a non-uniform, dataset-specific median-separation weight solution.",
    }.get(mode, "ADPS weighting mode could not be classified.")

    return f"""# dsRNASeeker label-free ranking-robustness report ({case_label})

## What this report can and cannot tell you

No RIP/J2/Z22/dsRNA-seq labels are used here. These diagnostics therefore do **not**
identify the evidence route with the highest true AP/ROC-AUC and do not estimate
the probability that a candidate forms a duplex. They quantify **ranking
stability and evidence-profile agreement** under pre-specified label-free views.

## Quick-read status

- **ADPS_WEIGHT_STATUS:** `{weight_status}`
- **ADPS_SCORE_RESOLUTION:** `{score_resolution}`
- **DATASET_USER_GUIDANCE:** `{dataset_guidance}`

These labels summarize provenance and ranking granularity only. They are not
experimental confidence classes and were not calibrated against RIP labels.

## ADPS provenance

- weighting mode: **{mode}**
- internal adaptive-reference candidates: **{fmt(ref_n, 0)}**
- internal adaptive-reference fraction: **{fmt(ref_frac)}**
- largest ADPS component weight: **{fmt(max_w)}** ({dom})
- effective number of weighted components: **{fmt(neff)}**
- distinct ADPS score values: **{fmt(uniq, 0)} / {fmt(n, 0)}** scored candidates
- median Spearman correlation of ADPS with the other broad profiles: **{fmt(corr)}**

{mode_explanation}

## Practical interpretation for an unlabeled study

1. Treat ADPS as the pre-specified label-independent reference ranking, not as a
   calibrated probability.
2. Start with the plain-language columns in the candidate table:
   `CANDIDATE_RANK_STABILITY` and `USER_INTERPRETATION`.
   `ROBUST_MULTI_PROFILE_PRIORITY` means the candidate is in the top 10% under
   at least three broad evidence views with tight rank dispersion. It does **not**
   mean experimentally confirmed dsRNA.
3. For details, inspect `broad_profile_top10_count`,
   `broad_profile_percentile_median`, and `broad_profile_percentile_iqr`.
   If a candidate is high in one profile but has a wide cross-profile rank
   spread, treat it as **profile-sensitive**. Inspect the component and route
   columns rather than assuming the ADPS route is correct.
4. If the biological question was specified in advance (for example,
   structure-focused or editing-marked dsRNA), the corresponding route score may
   be used as a hypothesis-specific view **alongside** ADPS. Do not select a route
   retrospectively because it would have produced a better unknown validation
   result.
5. When candidate-level experimental labels are available, use the supervised
   benchmarking path with grouped held-out evaluation. When labels are absent,
   orthogonal validation is especially important for profile-sensitive
   candidates.

The `tight`, `moderate`, and `wide` agreement bands are descriptive cutoffs for
rank dispersion only. They are not experimentally calibrated confidence classes.
"""


def write_ranking_robustness(
    df: pd.DataFrame,
    *,
    outdir: str | Path,
    case_label: str,
    dominant_weight_threshold: float = DEFAULT_DOMINANT_WEIGHT_THRESHOLD,
    tight_iqr: float = DEFAULT_TIGHT_IQR,
    wide_iqr: float = DEFAULT_WIDE_IQR,
) -> dict[str, Path]:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)

    candidate, summary, corr = add_ranking_robustness(
        df,
        dominant_weight_threshold=dominant_weight_threshold,
        tight_iqr=tight_iqr,
        wide_iqr=wide_iqr,
    )

    paths = {
        "candidate": out / f"TEpair_dsRNA_ranking_robustness.{case_label}.csv",
        "summary": out / f"TEpair_dsRNA_ranking_robustness_summary.{case_label}.csv",
        "correlations": out / f"TEpair_dsRNA_profile_rank_correlations.{case_label}.csv",
        "guidance": out / f"TEpair_dsRNA_user_guidance.{case_label}.md",
    }
    candidate.to_csv(paths["candidate"], index=False)
    summary.to_csv(paths["summary"], index=False)
    corr.to_csv(paths["correlations"], index=False)
    paths["guidance"].write_text(
        _guidance_markdown(summary, case_label=case_label),
        encoding="utf-8",
    )

    # Machine-readable provenance/definitions for reproducibility.
    definitions = {
        "label_free": True,
        "replaces_adps": False,
        "broad_profile_definitions": {
            "integrated_adps": ["stored adaptive_priority_score"],
            **BROAD_PROFILE_DEFINITIONS,
        },
        "route_definitions": ROUTE_DEFINITIONS,
        "dominant_weight_threshold": dominant_weight_threshold,
        "tight_iqr": tight_iqr,
        "wide_iqr": wide_iqr,
        "human_readable_status_cutoffs": {
            "min_reference_n": DEFAULT_MIN_REFERENCE_N,
            "min_background_n": DEFAULT_MIN_BACKGROUND_N,
            "min_reference_fraction": DEFAULT_MIN_REFERENCE_FRACTION,
            "max_reference_fraction": DEFAULT_MAX_REFERENCE_FRACTION,
            "coarse_unique_n": DEFAULT_COARSE_UNIQUE_N,
            "moderate_unique_n": DEFAULT_MODERATE_UNIQUE_N,
            "coarse_unique_fraction": DEFAULT_COARSE_UNIQUE_FRACTION,
            "moderate_unique_fraction": DEFAULT_MODERATE_UNIQUE_FRACTION,
            "coarse_largest_tie_fraction": DEFAULT_COARSE_LARGEST_TIE_FRACTION,
            "moderate_largest_tie_fraction": DEFAULT_MODERATE_LARGEST_TIE_FRACTION,
            "note": (
                "Operational, label-free interpretation cutoffs only; not fitted "
                "to RIP labels and not accuracy/confidence thresholds."
            ),
        },
        "status_definitions": {
            "ADPS_WEIGHT_STATUS": {
                "ADAPTIVE_SUPPORTED": "Non-uniform ADPS weights with non-trivial internal reference and background groups.",
                "LIMITED_REFERENCE": "Non-uniform ADPS weights, but the internal reference/background support is small or highly imbalanced.",
                "EQUAL_WEIGHT_FALLBACK": "Dataset-specific adaptive weights were unavailable; equal weights were used.",
                "SINGLE_COMPONENT_DOMINATED": "At least one ADPS component carries the configured dominant-weight threshold.",
            },
            "ADPS_SCORE_RESOLUTION": {
                "HIGH": "Many distinct ADPS values with limited tying.",
                "MODERATE": "Intermediate ADPS score granularity.",
                "COARSE": "Few distinct ADPS values, very low unique-score fraction, or a large tied-score group.",
            },
            "CANDIDATE_RANK_STABILITY": {
                "TIGHT": "Small cross-profile rank IQR.",
                "MODERATE": "Intermediate cross-profile rank IQR.",
                "WIDE": "Large cross-profile rank IQR; profile-sensitive candidate.",
            },
            "USER_INTERPRETATION": {
                "ROBUST_MULTI_PROFILE_PRIORITY": "Top-10% rank in at least three broad profiles with tight dispersion; robust computational priority only.",
                "ADPS_SUPPORTED_PROFILE_SENSITIVE": "ADPS provenance is supported but the candidate rank changes widely across profiles.",
                "ADPS_SUPPORTED_GENERAL": "ADPS provenance is supported; candidate is not in the strict multi-profile priority class and is not widely profile-sensitive.",
                "LIMITED_REFERENCE_RANKING": "Adaptive ADPS used a limited internal reference/background basis.",
                "FALLBACK_RANKING": "ADPS is the equal-weight fallback in this dataset.",
                "SINGLE_COMPONENT_DRIVEN": "ADPS is dominated by one evidence component in this dataset.",
            },
        },
        "warning": (
            "Agreement is not accuracy. These diagnostics do not identify the "
            "unknown RIP-optimal evidence route in an unlabeled dataset."
        ),
    }
    paths["definitions"] = out / f"TEpair_dsRNA_ranking_robustness_definitions.{case_label}.json"
    paths["definitions"].write_text(
        json.dumps(definitions, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return paths


def run_robustness(args) -> None:
    summary_path = Path(args.summary_in)
    if not summary_path.is_file():
        raise FileNotFoundError(summary_path)
    df = pd.read_csv(summary_path)

    case = str(getattr(args, "case_label", "") or "").strip()
    if not case:
        case = _infer_case_label(summary_path)

    outdir = Path(getattr(args, "output_dir", "") or summary_path.parent)
    paths = write_ranking_robustness(
        df,
        outdir=outdir,
        case_label=case,
        dominant_weight_threshold=float(
            getattr(args, "dominant_weight_threshold", DEFAULT_DOMINANT_WEIGHT_THRESHOLD)
        ),
        tight_iqr=float(getattr(args, "agreement_tight_iqr", DEFAULT_TIGHT_IQR)),
        wide_iqr=float(getattr(args, "agreement_wide_iqr", DEFAULT_WIDE_IQR)),
    )
    for key, path in paths.items():
        print(f"[ROBUSTNESS] wrote {key}: {path}")
