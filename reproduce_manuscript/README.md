# Manuscript reproduction package

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
and benchmark documentation. The in-house mouse-liver sequencing accession
will be added to the manuscript and reproduction documentation once its
public deposition is finalized.
