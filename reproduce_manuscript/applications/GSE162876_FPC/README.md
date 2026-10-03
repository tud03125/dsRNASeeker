# GSE162876 FPC mouse-liver application

This directory contains the frozen manuscript-facing materials for the
public GSE162876 mouse-liver RNA-seq application, part of SuperSeries
GSE162878.

Contrasts:
- FPC Mild versus FPC Control
- FPC Advanced versus FPC Control

The analysis uses public HEP-INTACT hepatocyte-nuclei RNA-seq samples.
dsRNASeeker candidate discovery used the manuscript 1-kb TE-pair search
radius on mm39.

Hallmark gene-set enrichment was performed on the full ranked
gene-expression universe. dsRNASeeker candidate-associated genes were
subsequently annotated within Hallmark leading-edge subsets; Hallmark
enrichment itself was not calculated only from dsRNASeeker candidate genes.

Contents:
- GSE162876_FPC_sample_manifest.tsv: frozen sample/accession mapping
- candidate_summaries/: dsRNASeeker candidate summaries
- outputs/: DESeq2, ranked-gene, fgsea, metadata, and plot outputs
- scripts/: public-safe analysis scripts/wrapper

The figure outputs provide biological-context analyses rather than
candidate-level validation benchmarks.

To rerun the application, provide MATRIX, GTF, and DSRROOT as environment
variables. PYTHON_EXE and RSCRIPT_EXE may optionally be overridden.
