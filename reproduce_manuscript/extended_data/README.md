# dsRNASeeker — Bioinformatics Extended Data, manuscript supporting data (November 2026)

Data S1 contains **96 rows**, including 36 corrected primary-study top-k rows (three current studies, six methods, K=10/20) and 60 preserved earlier secondary-study rows (six separate historical analysis variants, five methods, K=10/20). The `source_lineage` column distinguishes these. A difference in label definitions or method panels across cohorts must not be mistaken for a unified independent validation. Primary ADPS FLAG P@10/P@20 = 0.90/0.70.

Data S2 retains **all nine** original rows byte-for-byte: the three principal records agree with the October corrected/rechecked source for matched counts, class counts, prevalence, odds ratio and Fisher p. The other six are historically archived and were not independently rechecked by the newer three-study source.

Data S3 (eight Z22 5-kb case-specific evidence rows) and Data S4 (three FLAG candidate-generation radius rows) are byte-identical to the original upload and validated against anchor values.

The primary source CSVs are preserved byte-for-byte outside extended_data in 90_SOURCE_PROVENANCE/. No source-only or LOSO benchmarks are included. The separate current Figure S1/Table S2 historical noNullZ radius sensitivity and 20-edge directed transfers already have their own supplementary figure/table locations.

Primary and historical secondary benchmark cohorts are distinguished by their documented source lineage. Historical secondary results should not be interpreted as the current principal benchmark results.


## Data S5 — Benchmark sensitivity and statistical comparisons
Data S5 (`DataS5_benchmark_sensitivity_and_statistical_comparisons.xlsx`) contains named worksheets. Data S1–S4 retain their original bytes. Detailed candidate-level inputs, sequence records and scripts are retained separately in the reproducibility archive. The radius-series J2 values use the same principal label cohort, but distinct no-shuffle scores. FLAG and Z22 radius-series performance is reproduced from the earlier Table S2 at its reported precision.
