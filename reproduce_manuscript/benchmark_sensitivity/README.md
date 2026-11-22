# Benchmark sensitivity and statistical comparisons

These materials support the November 2026 dsRNASeeker manuscript.
The original principal benchmark was produced for v1.0.4; later radius-sensitivity and execution-reliability revisions are archived here.

## Contents

- 01_benchmarks: candidate-level scores, individual-channel
  benchmarks, grouped bootstrap comparisons, precision at k,
  and genomic-orientation sensitivity.
- 02_radius: search-radius statistics and shared-pair rank stability.
- 03_transfer: archived cross-study predictions, reference-arm
  sequences, sequence-screening records and performance summaries.
- 04_resources: measured FASTQ-workflow computational resources.
- 05_shared_locus_cv: Z22 shared-TE-locus validation results.

## Reproduction

Run the following commands from the repository root using
the documented dsRNASeeker Python environment.

    B=reproduce_manuscript/benchmark_sensitivity
    OUT=reproduce_manuscript/reproduced

    python3 "$B/scripts/01_analyze.py" \
      --inputs "$B/01_benchmarks/inputs" \
      --out "$OUT/benchmarks" \
      --bootstrap 5000 \
      --seed 20261114

    python3 "$B/scripts/plot_Figure2.py" \
      --support "$B" \
      --out "$OUT/Figure2"

    python3 "$B/scripts/plot_FigureS1.py" \
      --support "$B" \
      --out "$OUT/FigureS1"

    python3 "$B/scripts/verify_transfer_metrics.py" \
      --transfer "$B/03_transfer"

## Scope

The scripts reproduce performance statistics and figures from
archived candidate-level scores and processed results.

They do not repeat FASTQ processing, reference-sequence extraction,
supervised training or the original BLAST searches.

The strictly inverted analysis evaluates existing scores and
predictions without refitting the models.

The J2 radius series uses the matched no-shuffle scoring
configuration with radius-specific positive/negative labels:
642/575 at 1 kb, 1018/1036 at 2 kb, and 1737/1877 at 5 kb.
FLAG and Z22 radius-series values retain their reported precision.

Historical source-path prefixes have been anonymized in public
metadata. Historical source checksums refer to the original files,
not to the edited public metadata.

The complete software workflow and previous reproduction materials
are retained elsewhere in this repository.
