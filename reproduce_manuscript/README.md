# Manuscript reproduction package

## Current manuscript benchmark materials

The November 2026 manuscript uses the updated statistical comparisons,
individual-channel baselines, fixed-cohort radius analysis, sequence-screened
cross-study transfers and genomic-orientation sensitivity results.

The current reproduction instructions are in:

`reproduce_manuscript/benchmark_sensitivity/README.md`

Current Figure 2 and Figure S1 files are named `Figure2_FINAL.*` and
`FigureS1_FINAL.*` under `reproduce_manuscript/figures/`.

Extended Data S1-S5 are in `reproduce_manuscript/extended_data/`.

The earlier grouped-fold Figure 2 analysis is retained for
historical reproducibility and is not the current manuscript figure.

This directory contains the public-safe manuscript-facing inputs, derived
benchmark tables, scripts, and frozen figure outputs used for the revised
dsRNASeeker analyses.

Detailed reproduction commands, benchmark definitions, and provenance are
documented in:

    docs/BENCHMARK_REPRODUCTION.md

from the repository root.

Canonical manuscript-facing rendered figures are stored under:

    reproduce_manuscript/figures/

Newly regenerated outputs should be written under:

    reproduce_manuscript/reproduced/

That directory is intentionally excluded from version control so that
reproduction runs do not modify the frozen source bundle.

Public sequencing datasets are identified by accession in the manuscript
and benchmark documentation. The public mouse-liver FPC application uses
GSE162876, part of SuperSeries GSE162878, with frozen manuscript-facing
materials archived under:

    reproduce_manuscript/applications/GSE162876_FPC/
