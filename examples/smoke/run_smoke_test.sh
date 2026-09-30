#!/usr/bin/env bash
set -euo pipefail
REPO="${1:-$PWD}"
cd "$REPO"
IN=examples/smoke/adps_input.tsv
EXP=examples/smoke/expected/adps_ranked.tsv
[[ -s "$IN" && -s "$EXP" ]] || { echo "[ERROR] smoke fixture/expected output missing"; exit 2; }

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

python3 - "$IN" "$tmp/observed.tsv" <<'PY'
import sys
import pandas as pd
from modules.priority import add_priority_columns
inp,out=sys.argv[1:]
d=pd.read_csv(inp,sep="\t")
r=add_priority_columns(d,case="CASE",control="CONTROL",
                       score_mode="adaptive",annotation_policy="conservative")
keep=[c for c in ["pair_id","adaptive_priority_score","rank_score","priority_tier",
                  "orientation_adps","annotation_adps","case_expression_adps",
                  "energy_adps","interface_adps","case_editing_adps","RI_adps"] if c in r.columns]
r[keep].sort_values(["rank_score","pair_id"],ascending=[False,True]).to_csv(
    out,sep="\t",index=False,float_format="%.12g")
PY

diff -u "$EXP" "$tmp/observed.tsv"
echo "[PASS] deterministic ADPS smoke test"

# Supervised branch smoke test on the bundled example.
python3 - <<'PY'
import pandas as pd
from webapp.app.supervised import same_study_rerank
p="webapp/examples/example_supervised_audit.tsv"
d=pd.read_csv(p,sep="\t")
ranked,folds=same_study_rerank(d,outer_folds=3,seed=20260925)
assert len(ranked)>0
assert ranked["supervised_score"].notna().all()
assert len(folds)>0
print("[PASS] grouped same-study supervised smoke test", len(ranked), "rows")
PY

python3 main.py --help >/dev/null
python3 main.py workflow --help >/dev/null
echo "[PASS] CLI import/help smoke test"
