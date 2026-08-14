from __future__ import annotations
# DSRNASEEKER_T2T_RI_GENEID_FALLBACK_V3
# DSRNASEEKER_PAIR_ID_RESOLUTION_FIX_V3
import os
import json
import re
from pathlib import Path
import tempfile
import subprocess
import numpy as np
import pandas as pd
from modules.priority import add_priority_columns, priority_front_columns
from modules.supervised import apply_supervised_priority
from modules.robustness import write_ranking_robustness


# DSRNASEEKER_PAIR_ID_RESOLUTION_FIX_V2

def _nonempty_text(value) -> str:
    """Normalize an identifier-like value without turning missing data into text."""
    if pd.isna(value):
        return ""
    text = str(value)
    return "" if text.strip().lower() in {"", "nan", "na", "none", "null"} else text


def _pair_id_splits(pair_id: str):
    """Yield every possible split around an overlapping '__' delimiter."""
    start = 0
    while True:
        index = pair_id.find("__", start)
        if index < 0:
            return
        yield pair_id[:index], pair_id[index + 2:]
        # Advance by one, not two: an arm ID ending in '_' can create an
        # overlapping delimiter at the true A/B boundary.
        start = index + 1


def _resolve_pair_arm_ids(
    pair_ids: pd.Series,
    known_te_ids,
    explicit_a: pd.Series | None = None,
    explicit_b: pd.Series | None = None,
) -> pd.DataFrame:
    """Resolve exact A/B TE IDs without assuming a particular ID naming scheme.

    pair_id is constructed as A_TE_ID + '__' + B_TE_ID.  Repeat identifiers may
    themselves contain underscores, may end in an underscore, and are not
    guaranteed to end in a genomic strand sign.  The only unambiguous parser is
    therefore to test possible delimiters against the exact TE identifiers from
    the input TE table.

    Future per-condition master files carry explicit A_TE_id/B_TE_id columns;
    those values are preferred and validated when available.
    """
    known = {_nonempty_text(value) for value in known_te_ids}
    known.discard("")

    pair_ids = pair_ids.astype(str)
    if explicit_a is None:
        explicit_a = pd.Series(pd.NA, index=pair_ids.index, dtype="object")
    if explicit_b is None:
        explicit_b = pd.Series(pd.NA, index=pair_ids.index, dtype="object")

    resolved_a = []
    resolved_b = []
    unmatched = []
    ambiguous = []

    for row_index, pair_id, a_value, b_value in zip(
        pair_ids.index, pair_ids, explicit_a, explicit_b
    ):
        a_text = _nonempty_text(a_value)
        b_text = _nonempty_text(b_value)

        if a_text and b_text:
            if a_text in known and b_text in known:
                resolved_a.append(a_text)
                resolved_b.append(b_text)
                continue
            unmatched.append(
                {
                    "row": row_index,
                    "pair_id": pair_id,
                    "reason": "explicit arm IDs are absent from the TE table",
                    "A_TE_id": a_text,
                    "B_TE_id": b_text,
                }
            )
            resolved_a.append(pd.NA)
            resolved_b.append(pd.NA)
            continue

        candidates = {
            (left, right)
            for left, right in _pair_id_splits(pair_id)
            if left in known and right in known
        }

        if len(candidates) == 1:
            left, right = next(iter(candidates))
            resolved_a.append(left)
            resolved_b.append(right)
        elif len(candidates) == 0:
            unmatched.append(
                {
                    "row": row_index,
                    "pair_id": pair_id,
                    "reason": "no delimiter produced two exact TE-table IDs",
                }
            )
            resolved_a.append(pd.NA)
            resolved_b.append(pd.NA)
        else:
            ambiguous.append(
                {
                    "row": row_index,
                    "pair_id": pair_id,
                    "candidate_splits": sorted(candidates),
                }
            )
            resolved_a.append(pd.NA)
            resolved_b.append(pd.NA)

    if unmatched or ambiguous:
        details = []
        if unmatched:
            details.append(f"unmatched examples={unmatched[:5]}")
        if ambiguous:
            details.append(f"ambiguous examples={ambiguous[:5]}")
        raise ValueError(
            "Could not resolve exact TE-pair arm identifiers using the input "
            "TE table. " + "; ".join(details)
        )

    return pd.DataFrame(
        {"A_TE_id": resolved_a, "B_TE_id": resolved_b},
        index=pair_ids.index,
    )


def run_summary(args) -> None:
    work = Path(args.output_dir)
    tag = args.analyze_subset
    case = args.case_label
    control = args.control_label

    case_dedup = work / case / tag / case / f"TEpair_dsRNA_master.{case}.dedup.tsv"
    ctrl_dedup = work / control / tag / control / f"TEpair_dsRNA_master.{control}.dedup.tsv"
    if not case_dedup.exists() or not ctrl_dedup.exists():
        raise FileNotFoundError(
            f"Missing per-condition dedup masters\n  {case_dedup}\n  {ctrl_dedup}"
        )

    case_df = pd.read_csv(case_dedup, sep="\t")
    ctrl_df = pd.read_csv(ctrl_dedup, sep="\t")
    H = case_df.add_suffix("_H")
    F = ctrl_df.add_suffix("_F")
    master = H.merge(F, left_on="pair_id_H", right_on="pair_id_F", how="outer")

    def coalesce(colH: str, colF: str):
        vH = master[colH] if colH in master.columns else pd.Series([np.nan] * len(master), index=master.index)
        vF = master[colF] if colF in master.columns else pd.Series([np.nan] * len(master), index=master.index)
        return vH.combine_first(vF)

    M = pd.DataFrame(index=master.index.copy())
    M["pair_id"] = coalesce("pair_id_H", "pair_id_F")
    for base in [
        "A_TE_id", "B_TE_id",
        "A_SYMBOL", "B_SYMBOL", "A_annotation", "B_annotation",
        "A_repFamily", "A_repName", "B_repFamily", "B_repName",
        "genomic_orientation", "transcript_orientation",
    ]:
        M[base] = coalesce(f"{base}_H", f"{base}_F")

    for base in [
        "RNAcofold_MFE_kcalmol", "MFE_norm_kcalpermkb",
        "RNAfold_A_MFE_kcalmol", "RNAfold_B_MFE_kcalmol",
        "ddG_interaction_kcalmol", "ddG_norm_kcalpermkb", "ddG_Z",
        "ddG_mu_null", "ddG_sd_null", "null_n_requested", "null_n_effective",
        "null_shuffle_exact_dinuc", "null_shuffle_method",
        "interface_bpp_sum", "interface_bpp_max", "interface_bpp_n",
        "interface_bpp_n_ge_1e5", "interface_bpp_expected_fraction_shorter",
        "interface_bpp_mean_arm_fraction",
    ]:
        M[base] = coalesce(f"{base}_H", f"{base}_F")

    for base in ["total", "fwd_frac", "both_strands", "arm_opposite", "arms_both_cov"]:
        M[f"{case}_{base}"] = coalesce(f"{case}_{base}_H", f"{case}_{base}_F")
        M[f"{control}_{base}"] = coalesce(f"{control}_{base}_H", f"{control}_{base}_F")

    M[f"{case}_AtoI_hits_window"] = master["AtoI_hits_window_H"] if "AtoI_hits_window_H" in master.columns else np.nan
    M[f"{control}_AtoI_hits_window"] = master["AtoI_hits_window_F"] if "AtoI_hits_window_F" in master.columns else np.nan
    M[f"{case}_REDI_hits_window"] = master["REDI_hits_window_H"] if "REDI_hits_window_H" in master.columns else np.nan
    M[f"{control}_REDI_hits_window"] = master["REDI_hits_window_F"] if "REDI_hits_window_F" in master.columns else np.nan

    M["AtoI_hits_window"] = (
        pd.to_numeric(M[f"{case}_AtoI_hits_window"], errors="coerce").fillna(0)
        + pd.to_numeric(M[f"{control}_AtoI_hits_window"], errors="coerce").fillna(0)
    )
    M["REDI_hits_window"] = (
        pd.to_numeric(M[f"{case}_REDI_hits_window"], errors="coerce").fillna(0)
        + pd.to_numeric(M[f"{control}_REDI_hits_window"], errors="coerce").fillna(0)
    )

    for col in ["bias_penalty", "expr_points", "energy_points", "editing_points", "rank_score"]:
        M[col] = coalesce(f"{col}_H", f"{col}_F")

    M[f"dsRNA_confidence_{case}"] = master["dsRNA_confidence_H"] if "dsRNA_confidence_H" in master.columns else np.nan
    M[f"dsRNA_confidence_{control}"] = master["dsRNA_confidence_F"] if "dsRNA_confidence_F" in master.columns else np.nan

    order = {"high": 3, "probable": 2, "possible": 1, "uncertain": 0, None: -1, np.nan: -1}

    def combine_conf(h, f):
        rH = order.get(h, -1)
        rF = order.get(f, -1)
        if rH < 0 and rF < 0:
            return np.nan
        score = min(max(rH, 0), max(rF, 0))
        return {0: "uncertain", 1: "possible", 2: "probable", 3: "high"}[score]

    M["dsRNA_confidence"] = [
        combine_conf(h, f) for h, f in zip(M[f"dsRNA_confidence_{case}"], M[f"dsRNA_confidence_{control}"])
    ]

    csv = args.csv_in
    if csv and Path(csv).exists():
        te = pd.read_csv(csv)
        id_col = "Row.names" if "Row.names" in te.columns else te.columns[0]
        cols = [id_col] + [
            c for c in [
                "log2FoldChange", "padj", "repClass",
                "SYMBOL", "annotation", "geneId", "GENENAME", "ENSEMBL",
            ] if c in te.columns
        ]
        te = te[cols].drop_duplicates(id_col).rename(columns={id_col: "TE_id"})
        if "repClass" not in te.columns:
            te["repClass"] = np.nan

        # Resolve exact arm identifiers against the actual TE table instead of
        # assuming a naming suffix. atena IDs such as AluJb_dup128444,
        # B1F1_dup11564 and 7SK_dup613 do not end in genomic strand signs.
        # Future/current per-condition master files already carry explicit
        # A_TE_id/B_TE_id columns; these are preferred and validated. The
        # delimiter search is a compatibility fallback for older masters.
        explicit_a = M["A_TE_id"] if "A_TE_id" in M.columns else None
        explicit_b = M["B_TE_id"] if "B_TE_id" in M.columns else None
        ab = _resolve_pair_arm_ids(
            M["pair_id"],
            te["TE_id"],
            explicit_a=explicit_a,
            explicit_b=explicit_b,
        )
        M["A_TE_id"] = ab["A_TE_id"]
        M["B_TE_id"] = ab["B_TE_id"]

        A = te.rename(columns={
            "TE_id": "A_TE_id",
            "log2FoldChange": "A_log2FC",
            "padj": "A_padj",
            "repClass": "A_repClass",
            "SYMBOL": "A_SYMBOL_from_TE",
            "annotation": "A_annotation_from_TE",
            "geneId": "A_geneId_from_TE",
            "GENENAME": "A_GENENAME_from_TE",
            "ENSEMBL": "A_ENSEMBL_from_TE",
        })
        B = te.rename(columns={
            "TE_id": "B_TE_id",
            "log2FoldChange": "B_log2FC",
            "padj": "B_padj",
            "repClass": "B_repClass",
            "SYMBOL": "B_SYMBOL_from_TE",
            "annotation": "B_annotation_from_TE",
            "geneId": "B_geneId_from_TE",
            "GENENAME": "B_GENENAME_from_TE",
            "ENSEMBL": "B_ENSEMBL_from_TE",
        })
        M = M.merge(A, on="A_TE_id", how="left").merge(B, on="B_TE_id", how="left")

        def _clean_missing_text(series: pd.Series) -> pd.Series:
            out = series.copy()
            text = out.astype(str).str.strip()
            bad = out.isna() | text.str.lower().isin({"", "nan", "na", "n/a", "none", "null"})
            return out.mask(bad)

        # A corrected custom-genome TE annotation can be applied during summary
        # regeneration without rerunning coverage, energetics, editing, or the
        # per-condition pair construction.
        for arm in ["A", "B"]:
            for base in ["SYMBOL", "annotation"]:
                source_col = f"{arm}_{base}_from_TE"
                target_col = f"{arm}_{base}"
                if source_col not in M.columns:
                    continue
                source = _clean_missing_text(M[source_col])
                if target_col in M.columns:
                    target = _clean_missing_text(M[target_col])
                    M[target_col] = target.combine_first(source)
                else:
                    M[target_col] = source

        def pair_lfc(row):
            vals = [row.get("A_log2FC"), row.get("B_log2FC")]
            vals = [v for v in vals if pd.notna(v)]
            return np.nan if not vals else float(np.mean(vals))

        def pair_padj(row):
            vals = [row.get("A_padj"), row.get("B_padj")]
            vals = [v for v in vals if pd.notna(v)]
            return np.nan if not vals else float(np.min(vals))

        M["log2FoldChange"] = M.apply(pair_lfc, axis=1)
        M["padj"] = M.apply(pair_padj, axis=1)

        def side_call(fc, pj, padj_thr=0.05, lfc_thr=0.5):
            try:
                fc = float(fc)
                pj = float(pj)
            except (TypeError, ValueError):
                return "not_tested"
            if np.isnan(fc) or np.isnan(pj):
                return "not_tested"
            if pj > padj_thr:
                return "ns"
            if fc >= lfc_thr:
                return f"{case}_up"
            if fc <= -lfc_thr:
                return f"{control}_up"
            return "weak"

        M["A_side_call"] = [side_call(fc, pj) for fc, pj in zip(M.get("A_log2FC"), M.get("A_padj"))]
        M["B_side_call"] = [side_call(fc, pj) for fc, pj in zip(M.get("B_log2FC"), M.get("B_padj"))]
        primary = {f"{case}_up", f"{control}_up"}

        def summary_side(a, b):
            if a in primary and b in primary:
                if a == b:
                    return f"both_{a}"
                return "discordant"
            if a in primary and b not in primary:
                return f"A_{a}"
            if b in primary and a not in primary:
                return f"B_{b}"
            if a == b:
                return a
            return "mixed"

        M["summary_side"] = [summary_side(a, b) for a, b in zip(M["A_side_call"], M["B_side_call"])]

    # initialize RI columns.  The historical RI_* columns are kept, but now they
    # are slop-aware when --rmats-overlap-slop > 0.  Exact-overlap and nearest-RI
    # diagnostics are added separately so we can debug whether candidates are
    # inside, near, or only gene-associated with rMATS RI events.
    ri_numeric_defaults = {
        "RI_overlap_any": 0, "RI_overlap_W": 0, "RI_overlap_A": 0, "RI_overlap_B": 0, "RI_overlap_both_arms": 0,
        "RI_event_count_W": 0, "RI_event_count_A": 0, "RI_event_count_B": 0,
        "RI_exact_overlap_any": 0, "RI_exact_overlap_W": 0, "RI_exact_overlap_A": 0, "RI_exact_overlap_B": 0, "RI_exact_overlap_both_arms": 0,
        "RI_nearby_any": 0, "RI_nearby_W": 0, "RI_nearby_A": 0, "RI_nearby_B": 0, "RI_nearby_both_arms": 0,
        "RI_gene_match_any": 0, "RI_gene_match_count": 0,
        "RI_gene_symbol_match_any": 0, "RI_gene_id_match_any": 0,
        "RI_gene_symbol_match_count": 0, "RI_gene_id_match_count": 0,
        "RI_search_slop_bp": int(getattr(args, "rmats_overlap_slop", 0)),
    }
    ri_float_defaults = {
        "RI_min_FDR_W": np.nan, "RI_min_FDR_A": np.nan, "RI_min_FDR_B": np.nan,
        "RI_max_abs_dPSI_W": np.nan, "RI_max_abs_dPSI_A": np.nan, "RI_max_abs_dPSI_B": np.nan,
        "RI_gene_min_FDR": np.nan, "RI_gene_max_abs_dPSI": np.nan,
        "RI_nearest_distance_W": np.nan, "RI_nearest_FDR_W": np.nan, "RI_nearest_dPSI_W": np.nan,
    }
    ri_text_defaults = {
        "RI_direction_majority_W": "", "RI_direction_majority_A": "", "RI_direction_majority_B": "",
        "RI_gene_direction_majority": "",
        "RI_gene_match_type": "",
        "RI_gene_matched_symbols": "",
        "RI_gene_matched_ids": "",
        "RI_nearest_geneSymbol_W": "", "RI_nearest_ID_W": "",
        "RI_interval_mode": str(getattr(args, "rmats_interval_mode", "intron_body")),
    }
    for c, v in {**ri_numeric_defaults, **ri_float_defaults, **ri_text_defaults}.items():
        if c not in M.columns:
            M[c] = v

    ri_integration_diagnostics = {
        "rmats_dir": str(getattr(args, "rmats_dir", "") or ""),
        "rmats_track": str(getattr(args, "rmats_track", "") or ""),
        "rmats_fdr_max": float(getattr(args, "rmats_fdr_max", 0.05)),
        "rmats_interval_mode": str(getattr(args, "rmats_interval_mode", "intron_body")),
        "rmats_overlap_slop": int(getattr(args, "rmats_overlap_slop", 0)),
        "status": "not_requested",
        "total_RI_rows": 0,
        "numeric_FDR_dPSI_rows": 0,
        "significant_RI_rows": 0,
        "valid_significant_RI_intervals": 0,
        "candidate_pairs": int(len(M)),
        "candidate_pairs_with_coordinate_window_overlap": 0,
        "candidate_pairs_with_gene_symbol_match": 0,
        "candidate_pairs_with_gene_id_match": 0,
        "candidate_pairs_with_any_gene_match": 0,
        "note": (
            "Coordinate RI columns require shared sequence names plus overlapping/proximal numeric intervals. "
            "GeneID fallback is diagnostic/contextual and does not by itself pass priority_gate_case_RI."
        ),
    }

    if args.rmats_dir:
        rmats_file = Path(args.rmats_dir) / f"RI.MATS.{args.rmats_track}.txt"
        ri_integration_diagnostics["rmats_file"] = str(rmats_file)
        if rmats_file.exists() and csv and Path(csv).exists() and "A_TE_id" in M.columns:
            ri = pd.read_csv(rmats_file, sep="\t", dtype=str)
            ri["FDR_num"] = pd.to_numeric(ri.get("FDR", np.nan), errors="coerce")
            ri["dPSI_num"] = pd.to_numeric(ri.get("IncLevelDifference", np.nan), errors="coerce")
            # Optional sign correction for legacy rMATS runs where --b1/--b2 were
            # accidentally control/case but the summary should be interpreted as
            # case-minus-control. After this flip, positive dPSI means
            # rmats_group1_label has higher RI and negative dPSI means
            # rmats_group2_label has higher RI.
            if getattr(args, "rmats_flip_dpsi", False):
                ri["dPSI_num"] = -ri["dPSI_num"]
            ri = ri.dropna(subset=["FDR_num", "dPSI_num"])
            ri_integration_diagnostics["numeric_FDR_dPSI_rows"] = int(len(ri))
            ri_filt = ri[ri["FDR_num"] <= float(args.rmats_fdr_max)].copy()
            ri_integration_diagnostics["significant_RI_rows"] = int(len(ri_filt))
            if not ri_filt.empty:
                te_full = pd.read_csv(csv, sep=None, engine="python")
                id_col = "Row.names" if "Row.names" in te_full.columns else te_full.columns[0]
                chr_col = next((c for c in ["seqnames", "chr", "chrom", "Chromosome", "chromosome"] if c in te_full.columns), None)
                start_col = "start" if "start" in te_full.columns else ("Start" if "Start" in te_full.columns else None)
                end_col = "end" if "end" in te_full.columns else ("End" if "End" in te_full.columns else None)
                if chr_col and start_col and end_col:
                    te_coords = te_full[[id_col, chr_col, start_col, end_col]].copy().dropna()
                    te_coords[start_col] = pd.to_numeric(te_coords[start_col], errors="coerce")
                    te_coords[end_col] = pd.to_numeric(te_coords[end_col], errors="coerce")
                    te_coords = te_coords.dropna(subset=[start_col, end_col]).set_index(id_col)

                    def clean_symbol(x) -> str:
                        if pd.isna(x):
                            return ""
                        return str(x).strip().strip('"').strip("'")

                    def get_ri_interval(r):
                        """Return 0-based BED interval for a rMATS RI event.

                        Default uses the actual intron body upstreamEE -> downstreamES.
                        Fallback uses the full rMATS RI event span riExonStart_0base -> riExonEnd.
                        """
                        mode = str(getattr(args, "rmats_interval_mode", "intron_body"))
                        s = e = np.nan
                        if mode == "intron_body":
                            s = pd.to_numeric(r.get("upstreamEE"), errors="coerce")
                            e = pd.to_numeric(r.get("downstreamES"), errors="coerce")
                        if pd.isna(s) or pd.isna(e) or int(e) <= int(s):
                            s = pd.to_numeric(r.get("riExonStart_0base"), errors="coerce")
                            e = pd.to_numeric(r.get("riExonEnd"), errors="coerce")
                        if pd.isna(s) or pd.isna(e) or int(e) <= int(s):
                            return None
                        return int(s), int(e)

                    def clean_bed_field(x) -> str:
                        """Return a non-empty, tab-safe field for BEDTools I/O."""
                        y = clean_symbol(x)
                        y = y.replace("\t", " ").replace("\n", " ").replace("\r", " ").strip()
                        return y if y else "."

                    def normalize_gene_id(x) -> str:
                        y = clean_symbol(x)
                        if not y:
                            return ""
                        # Ensembl/GENCODE style version suffix; keep exact ID columns too.
                        return re.sub(r"\.\d+$", "", y)

                    ri_records = []
                    for _, r in ri_filt.iterrows():
                        iv = get_ri_interval(r)
                        if iv is None:
                            continue
                        gene = clean_bed_field(r.get("geneSymbol", ""))
                        gene_id = clean_bed_field(r.get("GeneID", ""))
                        ri_records.append({
                            "chr": clean_bed_field(r.get("chr", "")),
                            "start": iv[0],
                            "end": iv[1],
                            "ri_id": clean_bed_field(r.get("ID", "")),
                            "ri_fdr": float(r["FDR_num"]),
                            "ri_dpsi": float(r["dPSI_num"]),
                            "ri_gene": gene,
                            "ri_gene_id": gene_id,
                            "ri_gene_id_norm": normalize_gene_id(gene_id),
                        })
                    ri_integration_diagnostics["valid_significant_RI_intervals"] = int(len(ri_records))
                    if ri_records:
                        ri_df = pd.DataFrame(ri_records)

                        with tempfile.TemporaryDirectory() as td:
                            pair_bed = Path(td) / "pairs_window.slop.bed"
                            arm_bed = Path(td) / "pairs_arms.slop.bed"
                            pair_bed_exact = Path(td) / "pairs_window.exact.bed"
                            arm_bed_exact = Path(td) / "pairs_arms.exact.bed"
                            ri_bed = Path(td) / "rmats_RI.bed"

                            with ri_bed.open("w") as f:
                                for _, r in ri_df.iterrows():
                                    f.write(
                                        f"{r['chr']}\t{int(r['start'])}\t{int(r['end'])}\t"
                                        f"{r['ri_id']}\t{r['ri_fdr']}\t{r['ri_dpsi']}\t{r['ri_gene']}\t{r['ri_gene_id']}\n"
                                    )

                            def get_coord(te_id):
                                if te_id in te_coords.index:
                                    row = te_coords.loc[te_id]
                                    # Some TE ids can be duplicated in unusual annotations. Keep the first row.
                                    if isinstance(row, pd.DataFrame):
                                        row = row.iloc[0]
                                    return str(row[chr_col]), int(float(row[start_col]) - 1), int(float(row[end_col]))
                                return None

                            slop = max(0, int(getattr(args, "rmats_overlap_slop", 0)))

                            def pad_interval(s, e, pad):
                                return max(0, int(s) - pad), int(e) + pad

                            with pair_bed.open("w") as fw, arm_bed.open("w") as fa, pair_bed_exact.open("w") as fwe, arm_bed_exact.open("w") as fae:
                                for pid, a_id, b_id in zip(M["pair_id"].astype(str), M["A_TE_id"].astype(str), M["B_TE_id"].astype(str)):
                                    ca = get_coord(a_id)
                                    cb = get_coord(b_id)
                                    if ca is None or cb is None:
                                        continue
                                    chrA, sA, eA = ca
                                    chrB, sB, eB = cb
                                    if chrA != chrB:
                                        continue
                                    # Exact windows/arms.
                                    wS, wE = min(sA, sB), max(eA, eB)
                                    fwe.write(f"{chrA}\t{wS}\t{wE}\t{pid}\n")
                                    fae.write(f"{chrA}\t{sA}\t{eA}\t{pid}|A\n")
                                    fae.write(f"{chrB}\t{sB}\t{eB}\t{pid}|B\n")
                                    # Slop-aware windows/arms used for the main RI_* evidence columns.
                                    sA2, eA2 = pad_interval(sA, eA, slop)
                                    sB2, eB2 = pad_interval(sB, eB, slop)
                                    fw.write(f"{chrA}\t{min(sA2, sB2)}\t{max(eA2, eB2)}\t{pid}\n")
                                    fa.write(f"{chrA}\t{sA2}\t{eA2}\t{pid}|A\n")
                                    fa.write(f"{chrB}\t{sB2}\t{eB2}\t{pid}|B\n")

                            bedtools = getattr(args, "bedtools_exe", "bedtools")

                            def sort_bed(in_path: Path) -> Path:
                                """Return a coordinate-sorted copy of a BED file for bedtools.

                                bedtools closest/intersect can fail with "out of order record"
                                when a temporary BED is written in candidate-table order.
                                Sorting all temporary BEDs with the same key keeps chromosome
                                order consistent between -a and -b.
                                """
                                out_path = in_path.with_name(in_path.stem + ".sorted" + in_path.suffix)
                                with out_path.open("w") as out:
                                    p_sort = subprocess.run(
                                        ["sort", "-k1,1", "-k2,2n", str(in_path)],
                                        stdout=out,
                                        stderr=subprocess.PIPE,
                                        text=True,
                                        env={**os.environ, "LC_ALL": "C"},
                                    )
                                if p_sort.returncode != 0:
                                    raise RuntimeError(f"Failed to sort BED {in_path}:\n{p_sort.stderr}")
                                return out_path

                            # BEDTools closest requires coordinate-sorted BED input.
                            # Sort every temporary BED used by intersect/closest.
                            ri_bed_sorted = sort_bed(ri_bed)
                            pair_bed_sorted = sort_bed(pair_bed)
                            arm_bed_sorted = sort_bed(arm_bed)
                            pair_bed_exact_sorted = sort_bed(pair_bed_exact)
                            arm_bed_exact_sorted = sort_bed(arm_bed_exact)

                            def intersect(a_path: Path) -> pd.DataFrame:
                                p = subprocess.run([bedtools, "intersect", "-wa", "-wb", "-a", str(a_path), "-b", str(ri_bed_sorted)], capture_output=True, text=True)
                                if p.returncode != 0:
                                    raise RuntimeError(p.stderr)
                                if not p.stdout.strip():
                                    return pd.DataFrame()
                                rows = [ln.split("\t") for ln in p.stdout.strip().split("\n")]
                                df = pd.DataFrame(rows, columns=[
                                    "a_chr", "a_start", "a_end", "a_name",
                                    "b_chr", "b_start", "b_end", "ri_id", "ri_fdr", "ri_dpsi", "ri_gene", "ri_gene_id",
                                ])
                                df["ri_fdr"] = pd.to_numeric(df["ri_fdr"], errors="coerce")
                                df["ri_dpsi"] = pd.to_numeric(df["ri_dpsi"], errors="coerce")
                                return df

                            def closest(a_path: Path) -> pd.DataFrame:
                                p = subprocess.run([bedtools, "closest", "-d", "-t", "first", "-a", str(a_path), "-b", str(ri_bed_sorted)], capture_output=True, text=True)
                                if p.returncode != 0:
                                    raise RuntimeError(p.stderr)
                                if not p.stdout.strip():
                                    return pd.DataFrame()
                                rows = [ln.split("\t") for ln in p.stdout.strip().split("\n")]
                                # bedtools closest emits b fields as -1 / . if there is no match on the chromosome.
                                df = pd.DataFrame(rows, columns=[
                                    "a_chr", "a_start", "a_end", "a_name",
                                    "b_chr", "b_start", "b_end", "ri_id", "ri_fdr", "ri_dpsi", "ri_gene", "ri_gene_id", "distance",
                                ])
                                df["distance"] = pd.to_numeric(df["distance"], errors="coerce")
                                df["ri_fdr"] = pd.to_numeric(df["ri_fdr"], errors="coerce")
                                df["ri_dpsi"] = pd.to_numeric(df["ri_dpsi"], errors="coerce")
                                df = df[df["distance"].ge(0)]
                                return df

                            def majority_direction(series):
                                pos = (series > 0).sum()
                                neg = (series < 0).sum()
                                if pos == 0 and neg == 0:
                                    return ""
                                if pos > neg:
                                    return f"{args.rmats_group1_label}_high_RI"
                                if neg > pos:
                                    return f"{args.rmats_group2_label}_high_RI"
                                return "mixed"

                            def populate_hits(hit_df: pd.DataFrame, prefix: str, arm_specific: bool = False):
                                if hit_df.empty:
                                    return
                                if not arm_specific:
                                    g = hit_df.groupby("a_name")
                                    size = g.size()
                                    idx = size.index
                                    M.loc[M["pair_id"].isin(idx), f"RI_{prefix}_W"] = 1
                                    M[f"RI_event_count_W"] = M["pair_id"].map(size).combine_first(M[f"RI_event_count_W"]).fillna(M[f"RI_event_count_W"])
                                    M["RI_min_FDR_W"] = M["pair_id"].map(g["ri_fdr"].min()).combine_first(M["RI_min_FDR_W"])
                                    M["RI_max_abs_dPSI_W"] = M["pair_id"].map(g["ri_dpsi"].apply(lambda s: np.nanmax(np.abs(s.values)))).combine_first(M["RI_max_abs_dPSI_W"])
                                    M["RI_direction_majority_W"] = M["pair_id"].map(g["ri_dpsi"].apply(majority_direction)).combine_first(M["RI_direction_majority_W"])
                                else:
                                    hit_df = hit_df.copy()
                                    hit_df["pair_id"] = hit_df["a_name"].str.replace(r"\|[AB]$", "", regex=True)
                                    hit_df["arm"] = hit_df["a_name"].str.extract(r"\|([AB])$", expand=False)
                                    for arm in ["A", "B"]:
                                        sub = hit_df[hit_df["arm"] == arm]
                                        if sub.empty:
                                            continue
                                        g = sub.groupby("pair_id")
                                        size = g.size()
                                        idx = size.index
                                        M.loc[M["pair_id"].isin(idx), f"RI_{prefix}_{arm}"] = 1
                                        M[f"RI_event_count_{arm}"] = M["pair_id"].map(size).combine_first(M[f"RI_event_count_{arm}"]).fillna(M[f"RI_event_count_{arm}"])
                                        M[f"RI_min_FDR_{arm}"] = M["pair_id"].map(g["ri_fdr"].min()).combine_first(M[f"RI_min_FDR_{arm}"])
                                        M[f"RI_max_abs_dPSI_{arm}"] = M["pair_id"].map(g["ri_dpsi"].apply(lambda s: np.nanmax(np.abs(s.values)))).combine_first(M[f"RI_max_abs_dPSI_{arm}"])
                                        M[f"RI_direction_majority_{arm}"] = M["pair_id"].map(g["ri_dpsi"].apply(majority_direction)).combine_first(M[f"RI_direction_majority_{arm}"])

                            # Main RI_* evidence: slop-aware. When slop=0 this is exact behavior.
                            hitW = intersect(pair_bed_sorted)
                            hitA = intersect(arm_bed_sorted)
                            populate_hits(hitW, "overlap", arm_specific=False)
                            populate_hits(hitA, "overlap", arm_specific=True)

                            # Exact-overlap diagnostics are independent of slop and not used directly by priority.py.
                            hitW_exact = intersect(pair_bed_exact_sorted)
                            hitA_exact = intersect(arm_bed_exact_sorted)
                            if not hitW_exact.empty:
                                idx = hitW_exact.groupby("a_name").size().index
                                M.loc[M["pair_id"].isin(idx), "RI_exact_overlap_W"] = 1
                            if not hitA_exact.empty:
                                hitA_exact = hitA_exact.copy()
                                hitA_exact["pair_id"] = hitA_exact["a_name"].str.replace(r"\|[AB]$", "", regex=True)
                                hitA_exact["arm"] = hitA_exact["a_name"].str.extract(r"\|([AB])$", expand=False)
                                for arm in ["A", "B"]:
                                    idx = hitA_exact.loc[hitA_exact["arm"].eq(arm), "pair_id"].unique()
                                    M.loc[M["pair_id"].isin(idx), f"RI_exact_overlap_{arm}"] = 1

                            M["RI_exact_overlap_both_arms"] = ((M["RI_exact_overlap_A"] == 1) & (M["RI_exact_overlap_B"] == 1)).astype(int)
                            M["RI_exact_overlap_any"] = ((M["RI_exact_overlap_W"] == 1) | (M["RI_exact_overlap_A"] == 1) | (M["RI_exact_overlap_B"] == 1)).astype(int)

                            M["RI_overlap_both_arms"] = ((M["RI_overlap_A"] == 1) & (M["RI_overlap_B"] == 1)).astype(int)
                            M["RI_overlap_any"] = ((M["RI_overlap_W"] == 1) | (M["RI_overlap_A"] == 1) | (M["RI_overlap_B"] == 1)).astype(int)
                            ri_integration_diagnostics["candidate_pairs_with_coordinate_window_overlap"] = int((M["RI_overlap_W"] == 1).sum())

                            # Nearby means it is found by the slop-aware search but not exact.
                            for part in ["W", "A", "B"]:
                                M[f"RI_nearby_{part}"] = ((M[f"RI_overlap_{part}"] == 1) & (M[f"RI_exact_overlap_{part}"] == 0)).astype(int)
                            M["RI_nearby_both_arms"] = ((M["RI_nearby_A"] == 1) & (M["RI_nearby_B"] == 1)).astype(int)
                            M["RI_nearby_any"] = ((M["RI_nearby_W"] == 1) | (M["RI_nearby_A"] == 1) | (M["RI_nearby_B"] == 1)).astype(int)

                            # Nearest significant RI event to the exact pair window.
                            nearW = closest(pair_bed_exact_sorted)
                            if not nearW.empty:
                                nearW = nearW.sort_values(["a_name", "distance", "ri_fdr"]).drop_duplicates("a_name")
                                M["RI_nearest_distance_W"] = M["pair_id"].map(nearW.set_index("a_name")["distance"]).combine_first(M["RI_nearest_distance_W"])
                                M["RI_nearest_FDR_W"] = M["pair_id"].map(nearW.set_index("a_name")["ri_fdr"]).combine_first(M["RI_nearest_FDR_W"])
                                M["RI_nearest_dPSI_W"] = M["pair_id"].map(nearW.set_index("a_name")["ri_dpsi"]).combine_first(M["RI_nearest_dPSI_W"])
                                M["RI_nearest_geneSymbol_W"] = M["pair_id"].map(nearW.set_index("a_name")["ri_gene"]).combine_first(M["RI_nearest_geneSymbol_W"])
                                M["RI_nearest_ID_W"] = M["pair_id"].map(nearW.set_index("a_name")["ri_id"]).combine_first(M["RI_nearest_ID_W"])

                            # Gene-level fallback, kept separate from coordinate/nearby RI evidence.
                            # Standard assemblies usually have symbols; custom/T2T assemblies often
                            # preserve ENSMUSG/ENSG gene IDs while symbols are blank. We therefore
                            # match both symbol and GeneID, and record which route succeeded. Gene-only
                            # evidence is diagnostic/contextual and does not directly pass RI_adps.
                            if "ri_gene" in ri_df.columns:
                                ri_symbol_df = ri_df[ri_df["ri_gene"].ne(".")].copy()
                                ri_geneid_df = ri_df[ri_df["ri_gene_id"].ne(".")].copy()
                                ri_by_symbol = {gene: sub for gene, sub in ri_symbol_df.groupby("ri_gene") if gene and gene != "."}
                                ri_by_geneid = {gid: sub for gid, sub in ri_geneid_df.groupby("ri_gene_id") if gid and gid != "."}
                                ri_by_geneid_norm = {gid: sub for gid, sub in ri_geneid_df.groupby("ri_gene_id_norm") if gid}

                                gene_count = []
                                gene_symbol_count = []
                                gene_id_count = []
                                gene_min_fdr = []
                                gene_max_abs = []
                                gene_dir = []
                                gene_any = []
                                gene_symbol_any = []
                                gene_id_any = []
                                gene_match_type = []
                                matched_symbols = []
                                matched_ids = []

                                for _, row in M.iterrows():
                                    symbols = {clean_symbol(row.get("A_SYMBOL", "")), clean_symbol(row.get("B_SYMBOL", ""))}
                                    symbols.discard("")
                                    gene_ids_exact = {clean_symbol(row.get("A_geneId_from_TE", "")), clean_symbol(row.get("B_geneId_from_TE", ""))}
                                    gene_ids_exact.discard("")
                                    gene_ids_norm = {normalize_gene_id(x) for x in gene_ids_exact if normalize_gene_id(x)}

                                    symbol_hits = []
                                    symbol_hit_names = []
                                    for sym in sorted(symbols):
                                        if sym in ri_by_symbol:
                                            symbol_hits.append(ri_by_symbol[sym])
                                            symbol_hit_names.append(sym)

                                    id_hits = []
                                    id_hit_names = []
                                    id_norm_hit_names = []
                                    for gid in sorted(gene_ids_exact):
                                        if gid in ri_by_geneid:
                                            id_hits.append(ri_by_geneid[gid])
                                            id_hit_names.append(gid)
                                    for gid in sorted(gene_ids_norm):
                                        if gid in ri_by_geneid_norm and gid not in id_hit_names:
                                            id_hits.append(ri_by_geneid_norm[gid])
                                            id_norm_hit_names.append(gid)

                                    sub_hits = symbol_hits + id_hits
                                    if not sub_hits:
                                        gene_count.append(0); gene_symbol_count.append(0); gene_id_count.append(0)
                                        gene_min_fdr.append(np.nan); gene_max_abs.append(np.nan); gene_dir.append("")
                                        gene_any.append(0); gene_symbol_any.append(0); gene_id_any.append(0)
                                        gene_match_type.append("none"); matched_symbols.append(""); matched_ids.append("")
                                        continue

                                    sub = pd.concat(sub_hits, ignore_index=True).drop_duplicates(["ri_id", "ri_gene", "ri_gene_id", "ri_fdr", "ri_dpsi"])
                                    gene_count.append(int(len(sub)))
                                    gene_symbol_count.append(int(sum(len(x) for x in symbol_hits)))
                                    gene_id_count.append(int(sum(len(x) for x in id_hits)))
                                    gene_min_fdr.append(float(np.nanmin(sub["ri_fdr"].values)))
                                    gene_max_abs.append(float(np.nanmax(np.abs(sub["ri_dpsi"].values))))
                                    gene_dir.append(majority_direction(sub["ri_dpsi"]))
                                    gene_any.append(1)
                                    gene_symbol_any.append(1 if symbol_hits else 0)
                                    gene_id_any.append(1 if id_hits else 0)
                                    if symbol_hits and id_hits:
                                        gene_match_type.append("symbol_and_gene_id")
                                    elif symbol_hits:
                                        gene_match_type.append("symbol_exact")
                                    elif id_hit_names and not id_norm_hit_names:
                                        gene_match_type.append("gene_id_exact")
                                    else:
                                        gene_match_type.append("gene_id_version_normalized")
                                    matched_symbols.append(";".join(symbol_hit_names))
                                    matched_ids.append(";".join(id_hit_names + id_norm_hit_names))

                                M["RI_gene_match_any"] = gene_any
                                M["RI_gene_symbol_match_any"] = gene_symbol_any
                                M["RI_gene_id_match_any"] = gene_id_any
                                M["RI_gene_match_count"] = gene_count
                                M["RI_gene_symbol_match_count"] = gene_symbol_count
                                M["RI_gene_id_match_count"] = gene_id_count
                                M["RI_gene_match_type"] = gene_match_type
                                M["RI_gene_matched_symbols"] = matched_symbols
                                M["RI_gene_matched_ids"] = matched_ids
                                M["RI_gene_min_FDR"] = gene_min_fdr
                                M["RI_gene_max_abs_dPSI"] = gene_max_abs
                                M["RI_gene_direction_majority"] = gene_dir
                                ri_integration_diagnostics["candidate_pairs_with_gene_symbol_match"] = int(sum(gene_symbol_any))
                                ri_integration_diagnostics["candidate_pairs_with_gene_id_match"] = int(sum(gene_id_any))
                                ri_integration_diagnostics["candidate_pairs_with_any_gene_match"] = int(sum(gene_any))

    outdir = work / tag / "summary"
    outdir.mkdir(parents=True, exist_ok=True)
    try:
        ri_integration_diagnostics["final_RI_overlap_any_pairs"] = int((M.get("RI_overlap_any", 0) == 1).sum())
        ri_integration_diagnostics["final_RI_nearest_significant_pairs"] = int(pd.to_numeric(M.get("RI_nearest_distance_W", pd.Series(dtype=float)), errors="coerce").notna().sum())
        ri_integration_diagnostics["final_RI_gene_match_pairs"] = int((M.get("RI_gene_match_any", 0) == 1).sum())
        ri_integration_diagnostics["final_RI_gene_id_match_pairs"] = int((M.get("RI_gene_id_match_any", 0) == 1).sum())
        if ri_integration_diagnostics.get("status") == "loaded_rmats_RI":
            if ri_integration_diagnostics.get("significant_RI_rows", 0) == 0:
                ri_integration_diagnostics["status"] = "no_significant_RI_at_requested_FDR"
            elif ri_integration_diagnostics.get("valid_significant_RI_intervals", 0) == 0:
                ri_integration_diagnostics["status"] = "no_valid_significant_RI_intervals"
            elif ri_integration_diagnostics.get("final_RI_overlap_any_pairs", 0) == 0 and ri_integration_diagnostics.get("final_RI_nearest_significant_pairs", 0) == 0:
                ri_integration_diagnostics["status"] = "significant_RI_loaded_but_no_coordinate_matches"
            else:
                ri_integration_diagnostics["status"] = "RI_integration_completed"
        with (outdir / f"TEpair_dsRNA_RI_integration_diagnostics.{case}.json").open("w") as fh:
            json.dump(ri_integration_diagnostics, fh, indent=2, sort_keys=True)
            fh.write("\n")
    except Exception as e:
        print(f"[WARN] failed to write RI integration diagnostics: {e}")

    require_case_editing = bool(getattr(args, "require_case_editing", True))
    require_case_ri = bool(getattr(args, "require_case_ri", True))
    if getattr(args, "priority_mode", "strict") == "relaxed":
        require_case_editing = False
        require_case_ri = False

    M = M.dropna(subset=["pair_id"])
    score_mode = getattr(args, "priority_score_mode", "adaptive")
    M = add_priority_columns(
        M,
        case=case,
        control=control,
        require_case_editing=require_case_editing,
        require_case_ri=require_case_ri,
        score_mode=score_mode,
    )

    if score_mode == "supervised":
        M = apply_supervised_priority(M, args, outdir=outdir, case=case, control=control)

    # Drop ambiguous historical total-burden column names from final summary outputs.
    # Public-facing summaries expose explicit total and case-minus-control columns:
    #   SPRINT_total_hits_window, SPRINT_delta_case_minus_control
    #   REDI_total_hits_window,   REDI_delta_case_minus_control
    M = M.drop(columns=["AtoI_hits_window", "REDI_hits_window"], errors="ignore")

    front = [c for c in priority_front_columns(case, control) if c in M.columns]
    rest = [c for c in M.columns if c not in front]
    M = M[front + rest]

    ri_cols = [c for c in M.columns if c.startswith("RI_")]

    # Keep the default summary filenames reserved for non-supervised ranking
    # (adaptive/expert). In supervised mode, write separate files so an ML
    # benchmark run does not overwrite the adaptive summary produced earlier.
    # Examples:
    #   adaptive:   TEpair_dsRNA_master.summary.with_RI.csv
    #   supervised: TEpair_dsRNA_master_supervised.summary.with_RI.csv
    output_mode = str(getattr(args, "priority_score_mode", "adaptive")).lower()
    mode_suffix = "_supervised" if output_mode == "supervised" else ""

    no_ri = outdir / f"TEpair_dsRNA_master{mode_suffix}.summary.csv"
    with_ri = outdir / f"TEpair_dsRNA_master{mode_suffix}.summary.with_RI.csv"
    strict_path = outdir / f"TEpair_dsRNA_high_priority{mode_suffix}.{case}.strict.csv"
    topn_path = outdir / f"TEpair_dsRNA_high_priority{mode_suffix}.{case}.top{int(getattr(args, 'priority_top_n', 20))}.csv"
    relaxed_path = outdir / f"TEpair_dsRNA_high_priority{mode_suffix}.{case}.relaxed.csv"
    # Raw ADPS weight diagnostics: original metric/value layout.
    weights_path = outdir / f"TEpair_dsRNA_adaptive_weights.{case}.csv"
    # Human-readable ADPS weight diagnostics: one row per evidence block.
    weights_long_path = outdir / f"TEpair_dsRNA_adaptive_weights_long.{case}.csv"

    M.drop(columns=ri_cols, errors="ignore").to_csv(no_ri, index=False)
    M.to_csv(with_ri, index=False)

    adps_blocks = [
        "orientation_adps",
        "annotation_adps",
        "case_expression_adps",
        "energy_adps",
        "interface_adps",
        "case_editing_adps",
        "RI_adps",
    ]
    adps_labels = {
        "orientation_adps": "orientation",
        "annotation_adps": "annotation",
        "case_expression_adps": "case_expression",
        "energy_adps": "energy",
        "interface_adps": "interface",
        "case_editing_adps": "case_editing",
        "RI_adps": "RI",
    }
    meta_cols = [
        c for c in [
            "adaptive_weight_source",
            "adaptive_gate_positive_n",
            "adaptive_gate_background_n",
        ]
        if c in M.columns
    ]
    weight_cols = [f"adaptive_weight_{b}" for b in adps_blocks if f"adaptive_weight_{b}" in M.columns]
    if weight_cols:
        first = M.head(1).iloc[0]

        # Original raw metric/value output. This is useful for scripts that expect
        # a flat diagnostic listing of every ADPS weight/statistic.
        raw_cols = meta_cols + weight_cols
        raw_cols += [f"adaptive_separation_{b}" for b in adps_blocks if f"adaptive_separation_{b}" in M.columns]
        raw_cols += [f"adaptive_pos_median_{b}" for b in adps_blocks if f"adaptive_pos_median_{b}" in M.columns]
        raw_cols += [f"adaptive_bg_median_{b}" for b in adps_blocks if f"adaptive_bg_median_{b}" in M.columns]
        weight_summary = M[raw_cols].head(1).T.reset_index()
        weight_summary.columns = ["metric", "value"]
        weight_summary.to_csv(weights_path, index=False)

        # Human-readable table output. This has the same values as the raw file,
        # but organized as one row per evidence block.
        rows = []
        for block in adps_blocks:
            rows.append({
                "evidence_block": adps_labels[block],
                "positive_median": first.get(f"adaptive_pos_median_{block}", 0.0),
                "background_median": first.get(f"adaptive_bg_median_{block}", 0.0),
                "separation": first.get(f"adaptive_separation_{block}", 0.0),
                "adaptive_weight": first.get(f"adaptive_weight_{block}", 0.0),
                "adps_feature_column": block,
                "adaptive_weight_source": first.get("adaptive_weight_source", ""),
                "adaptive_gate_positive_n": first.get("adaptive_gate_positive_n", ""),
                "adaptive_gate_background_n": first.get("adaptive_gate_background_n", ""),
            })
        weight_table_columns = [
            "evidence_block",
            "positive_median",
            "background_median",
            "separation",
            "adaptive_weight",
            "adps_feature_column",
            "adaptive_weight_source",
            "adaptive_gate_positive_n",
            "adaptive_gate_background_n",
        ]
        pd.DataFrame(rows)[weight_table_columns].to_csv(weights_long_path, index=False)

    strict = M[M["priority_gate_pass"]].sort_values("rank_score", ascending=False)
    strict.to_csv(strict_path, index=False)
    strict.head(int(getattr(args, "priority_top_n", 20))).to_csv(topn_path, index=False)

    relaxed = M[M["dsRNA_case_priority"].isin([
        "case_high_priority",
        "case_supported_missing_RI_or_annotation",
        "case_TE_only",
    ])].sort_values("rank_score", ascending=False)
    relaxed.to_csv(relaxed_path, index=False)

    # Label-free ranking-robustness diagnostics are derived entirely from the
    # finished summary table. They do not rerun Step 5 coverage/strand calls,
    # RNAfold/RNAcofold/null-Z/interface calculations, editing, or rMATS.
    # Keep them in separate files so the historical primary summary remains
    # byte/column compatible with existing benchmark scripts.
    try:
        robustness_paths = write_ranking_robustness(
            M,
            outdir=outdir,
            case_label=case,
        )
    except Exception as e:
        robustness_paths = {}
        print(f"[WARN] ranking-robustness diagnostics were not written: {e}")

    print(f"[SUMMARY] wrote: {no_ri}")
    print(f"[SUMMARY] wrote: {with_ri}")
    print(f"[SUMMARY] wrote: {strict_path} ({len(strict)} rows)")
    print(f"[SUMMARY] wrote: {topn_path}")
    if weight_cols:
        print(f"[SUMMARY] wrote: {weights_path}")
        print(f"[SUMMARY] wrote: {weights_long_path}")
    print(f"[SUMMARY] wrote: {relaxed_path} ({len(relaxed)} rows)")
    for robustness_name, robustness_path in robustness_paths.items():
        print(f"[SUMMARY] wrote robustness {robustness_name}: {robustness_path}")
    if getattr(args, "priority_score_mode", "adaptive") == "supervised":
        print(f"[SUMMARY] supervised outputs written under: {outdir}")
