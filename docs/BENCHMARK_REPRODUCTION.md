# Benchmark and manuscript reproduction

This document describes the frozen manuscript-facing reproduction bundle.

## Principal labeled benchmarks

- GSE308488/9, FLAG-ZBP1, 1-kb candidate radius.
- GSE184962/4, Z22, 5-kb expanded quantitative benchmark.
- GSE244103/4, J2, 1-kb candidate radius.

Exact positive/confident-negative/unlabeled definitions are in:

`reproduce_manuscript/data/label_definitions_exact.tsv`

Verification of the reconstructed rules is in:

`reproduce_manuscript/data/label_definition_verification.tsv`

## Figure 2

The source candidate-level inputs are staged under:

`reproduce_manuscript/input_root/`

Reproduce Figure 2 from repository root:

```bash
python3 reproduce_manuscript/scripts/plot_Figure2_AugustStyle_5fold_v3_frozenfolds.py \
  --root reproduce_manuscript/input_root \
  --fold-assignments reproduce_manuscript/data/Figure2_outer_fold_assignments.tsv \
  --outdir reproduce_manuscript/reproduced/Figure2
```

The frozen fold assignments, mean/SD metric table, and original rendered figure are
also distributed under `reproduce_manuscript/data/` and
`reproduce_manuscript/figures/`.

## Figure 4 and sensitivity analyses

`reproduce_manuscript/data/Table_supervised_generalization.csv` contains the
same-study and held-out-family supervised generalization results.

`reproduce_manuscript/data/mapping_confidence_canonical_ADPS_metrics.csv`
contains the mapping-confidence sensitivity metrics.

`reproduce_manuscript/data/evidence_direction_metrics_GSE184962_5kb.csv`
contains the Z22 evidence-direction diagnostic.

Final rendered figures are supplied in `reproduce_manuscript/figures/`.

## Public-data provenance

Public sequence data are identified by GEO accessions in the manuscript and benchmark
registry. The in-house mouse-liver RNA-seq accession must be inserted here and in the
manuscript Data Availability Statement once deposition is complete.

## Scope of the smoke test

`examples/smoke/run_smoke_test.sh` is deliberately small. It verifies deterministic
ADPS scoring, the grouped supervised branch, and CLI imports. It is not a substitute
for the full biological benchmark reproduction above.
