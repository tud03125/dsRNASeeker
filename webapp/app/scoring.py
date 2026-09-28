from __future__ import annotations
import pandas as pd
from .bridge import load_dsrnaseeker

def calculate_adps(df: pd.DataFrame, case: str, control: str, annotation_policy: str="conservative") -> pd.DataFrame:
    _, add_priority_columns, _ = load_dsrnaseeker()
    out=add_priority_columns(df,case=case,control=control,score_mode="adaptive",annotation_policy=annotation_policy)
    front=[c for c in ["pair_id","adaptive_priority_score","rank_score","priority_tier","orientation_adps","annotation_adps","case_expression_adps","energy_adps","interface_adps","case_editing_adps","RI_adps"] if c in out.columns]
    rest=[c for c in out.columns if c not in front]
    return out[front+rest].sort_values("rank_score",ascending=False)
