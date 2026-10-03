# dsRNASeeker extended benchmark data

These files contain secondary benchmark diagnostics retained for
reproducibility but not presented as formatted Supplementary Tables
in the reader-facing manuscript.

## Data S1
`DataS1_tie_aware_precision_at_k.csv`

Tie-aware expected Precision@10 and Precision@20 for primary
dsRNASeeker and projected comparator scores. The expected value
averages over possible orderings when the top-k boundary intersects
a score tie.

## Data S2
`DataS2_dsRNAscan_binary_enrichment.csv`

Matched-versus-unmatched enrichment analysis for the binary
dsRNAscan two-arm comparator within each frozen labeled benchmark
universe.

## Data S3
`DataS3_GSE184962_evidence_direction.csv`

Feature-wise AP, ROC-AUC, grouped-bootstrap confidence intervals,
and association direction for the GSE184962/4 5-kb Z22 benchmark.
These values underlie Figure S3.

## Data S4
`DataS4_FLAG_candidate_generation_capture.tsv`

Candidate-generation capture of the 158 published GSE308488
FLAG-ZBP1-positive source clusters at 1-, 2-, and 5-kb TE-pair
search radii, evaluated before ADPS ranking.

These analyses are secondary diagnostics and are provided to support
complete reproduction and audit of the manuscript analyses.
