# Manuscript benchmark reproduction

Principal FLAG-ZBP1 (1 kb), J2 (1 kb) and condition-specific Z22 (5 kb) include the sequence-null statistic. FLAG retains the verified principal source. J2 and Z22 use the completed N=50 calculations. Candidate and label counts are unchanged. Radius sensitivity uses separate no-Null-Z scores and radius-specific cohorts. Transfer uses five study-specific scoring sources, not a uniform Null-Z experiment.

From the repository root:

```bash
B=reproduce_manuscript/benchmark_sensitivity
OUT=reproduce_manuscript/reproduced
python3 "$B/scripts/01_analyze.py" --inputs "$B/01_benchmarks/inputs" --out "$OUT/benchmarks" --bootstrap 5000 --seed 20261114
python3 "$B/scripts/plot_Figure2.py" --support "$B" --out "$OUT/Figure2"
python3 "$B/scripts/plot_FigureS1.py" --support "$B" --out "$OUT/FigureS1"
python3 "$B/scripts/verify_transfer_metrics.py" --transfer "$B/03_transfer"
python3 main.py supervised-transfer --manifest "$B/03_transfer/manifest.tsv" --output-dir "$OUT/transfer" --mode both --bootstrap 0 --seed 20261011
```

The statistical scripts reproduce fixed-score statistics and the curve data; frozen publication artwork is supplied separately. Transfer training inputs are projected from verified source files, with provenance and checksums. The within-study Full-7 held-out predictions and per-fold results are archived alongside the seven channel scores. Z22 shared-locus folds and predictions are in 05_shared_locus_cv. The code does not repeat upstream FASTQ alignment or BLAST in these reproduction commands.
