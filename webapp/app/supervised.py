from __future__ import annotations
import numpy as np
import pandas as pd
from .bridge import load_dsrnaseeker

REQ=["pair_id","label","group_id","orientation_adps","annotation_adps","case_expression_adps","energy_adps","interface_adps","case_editing_adps","RI_adps"]

def same_study_rerank(df: pd.DataFrame, outer_folds: int=5, seed: int=20260925) -> tuple[pd.DataFrame,pd.DataFrame]:
    missing=[c for c in REQ if c not in df.columns]
    if missing: raise ValueError(f"Missing columns: {missing}")
    d=df.copy(); d["label"]=pd.to_numeric(d["label"],errors="coerce"); d=d[d.label.isin([0,1])].copy(); d["label"]=d.label.astype(int)
    _, _, fixed_grouped_oof=load_dsrnaseeker()
    score, folds=fixed_grouped_oof(d,outer_folds=outer_folds,seed=seed)
    d["supervised_score"]=score
    return d.sort_values("supervised_score",ascending=False),folds
