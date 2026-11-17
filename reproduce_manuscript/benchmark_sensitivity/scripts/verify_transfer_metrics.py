#!/usr/bin/env python3
"""Verify archived cross-study performance without refitting models."""

import argparse
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

p = argparse.ArgumentParser()
p.add_argument("--transfer", type=Path, required=True)
a = p.parse_args()

root = a.transfer

meta = json.loads((root / "COMPLETE.json").read_text())

sha = hashlib.sha256(
    (root / "all_predictions.csv").read_bytes()
).hexdigest()

assert sha == meta["prediction_file_sha256"]

pred = pd.read_csv(root / "all_predictions.csv")

pred = pred.loc[
    pred.evaluation.eq("directed") &
    pred.panel.eq("raw7")
].copy()

screen = pd.read_csv(
    root / "TRANSFER_SEQUENCE_SCREEN_PER_PAIR.csv"
)

summary = pd.read_csv(
    root / "TRANSFER_SEQUENCE_SCREEN_SUMMARY.csv"
)

keys = ["train", "test", "pair_id"]

d = pred.merge(
    screen,
    on=keys,
    how="outer",
    validate="one_to_one",
    indicator=True
)

assert d["_merge"].eq("both").all()

assert len(summary) == 20
assert not summary.duplicated(["train", "test"]).any()
assert d.groupby(["train", "test"]).ngroups == 20

for r in summary.itertuples(index=False):

    full = d.loc[
        d.train.eq(r.train) &
        d.test.eq(r.test)
    ]

    retained = full.loc[
        full.screen.eq("NO_MATCH_SCREENED")
    ]

    assert len(full) == r.n_test
    assert len(retained) == r.retained_n
    assert retained.label.sum() == r.retained_positive
    assert len(retained) - retained.label.sum() == r.retained_negative

    assert full.screen.eq("EXACT_ARM_MATCH").sum() == r.exact_arm_match
    assert full.screen.eq("NEAR_ARM_MATCH").sum() == r.near_only_arm_match
    assert full.screen.eq("UNRESOLVED").sum() == r.unresolved

    observed = [
        average_precision_score(full.label, full.prediction),
        roc_auc_score(full.label, full.prediction),
        average_precision_score(retained.label, retained.prediction),
        roc_auc_score(retained.label, retained.prediction)
    ]

    expected = [
        r.full_AP, r.full_AUC,
        r.retained_AP, r.retained_AUC
    ]

    assert np.allclose(
        observed, expected,
        rtol=0, atol=1e-10
    )

    print("PASS:", r.train, "->", r.test)

print("\nPASS: All 20 directed transfers verified.")
