# dsRNASeeker CLI reference

Generated from commit `e9e6868826683d07c9f922321269ad480df7f571` on 2026-09-29T15:02:14-04:00.

```text
usage: dsRNASeeker [-h]
                   {workflow,run,summary,delta,check,zrna,molecule-model,supervised-benchmark,robustness}
                   ...

Condition-agnostic TE-pair dsRNA discovery pipeline.

positional arguments:
  {workflow,run,summary,delta,check,zrna,molecule-model,supervised-benchmark,robustness}
    workflow            Run end-to-end workflow from FASTQ/BAM through
                        alignment, TE analysis, rMATS, RNA editing,
                        dsRNASeeker, delta, and zrna.
    run                 Run one condition through module-driven execution.
    summary             Build fused summary across case/control using per-
                        condition outputs.
    delta               Build delta table across case/control outputs.
    check               Check runtime dependencies, files, and samplesheet.
    zrna                Annotate inverted TE-pair dsRNA candidates with A-form
                        support and Z-RNA propensity.
    molecule-model      Annotate candidates with conservative
                        intramolecular/intermolecular compatibility.
    supervised-benchmark
                        Run manifest-driven nested grouped and leave-one-
                        study-family-out supervised evaluation.
    robustness          Build label-free ranking-robustness diagnostics from
                        an existing summary; does not rerun alignment,
                        coverage, energetics, editing, or pairing.

options:
  -h, --help            show this help message and exit
```

## `workflow`

```text
usage: dsRNASeeker workflow [-h] --output-dir OUTPUT_DIR --case-label
                            CASE_LABEL --control-label CONTROL_LABEL
                            --samplesheet SAMPLESHEET
                            [--input-mode {fastq,bam}] --fasta FASTA --gtf GTF
                            [--te-gtf TE_GTF]
                            [--strandedness {auto,forward,reverse,unstranded,fr-firststrand,fr-secondstrand}]
                            [--paired] [--single-end]
                            [--read-length READ_LENGTH] [--threads THREADS]
                            [--force] [--quiet] [--verbose]
                            [--infer-fastq-records INFER_FASTQ_RECORDS]
                            [--infer-experiment-exe INFER_EXPERIMENT_EXE]
                            [--strandedness-bed STRANDEDNESS_BED]
                            [--strandedness-sample-size STRANDEDNESS_SAMPLE_SIZE]
                            [--stranded-threshold STRANDED_THRESHOLD]
                            [--unstranded-threshold UNSTRANDED_THRESHOLD]
                            [--strandedness-fallback {unstranded,forward,reverse,error}]
                            [--star-index STAR_INDEX] [--build-star-index]
                            [--sjdb-overhang SJDB_OVERHANG]
                            [--skip-te-analysis] [--te-mode {advanced}]
                            [--te-genome TE_GENOME]
                            [--te-rmsk-rds TE_RMSK_RDS]
                            [--te-force-rebuild-rmsk] [--te-use-strand]
                            [--te-ignore-strand]
                            [--te-yield-size TE_YIELD_SIZE]
                            [--te-min-max-count TE_MIN_MAX_COUNT]
                            [--te-shrink-type {ashr,none}]
                            [--te-txdb-package TE_TXDB_PACKAGE]
                            [--te-txdb-gtf TE_TXDB_GTF]
                            [--te-txdb-rds TE_TXDB_RDS]
                            [--te-orgdb-package TE_ORGDB_PACKAGE]
                            [--te-feature-type TE_FEATURE_TYPE]
                            [--te-attribute TE_ATTRIBUTE]
                            [--te-padj-max TE_PADJ_MAX]
                            [--te-lfc-min TE_LFC_MIN]
                            [--te-candidate-mode {strict,expressed}]
                            [--te-candidate-min-mean TE_CANDIDATE_MIN_MEAN]
                            [--te-candidate-padj-max TE_CANDIDATE_PADJ_MAX]
                            [--te-candidate-lfc-min TE_CANDIDATE_LFC_MIN]
                            [--featurecounts-exe FEATURECOUNTS_EXE]
                            [--skip-rmats] [--rmats-exe RMATS_EXE]
                            [--rmats-track {JC,JCEC}]
                            [--rmats-fdr-max RMATS_FDR_MAX]
                            [--rmats-overlap-slop RMATS_OVERLAP_SLOP]
                            [--rmats-interval-mode {intron_body,event_span}]
                            [--rmats-cstat RMATS_CSTAT]
                            [--rmats-min-intron-length RMATS_MIN_INTRON_LENGTH]
                            [--rmats-max-exon-length INT|auto]
                            [--rmats-libtype {None,fr-unstranded,fr-firststrand,fr-secondstrand}]
                            [--rmats-novel-ss] [--no-rmats-novel-ss]
                            [--skip-reditools] [--reditools-exe REDITOOLS_EXE]
                            [--reditools-strand {auto,0,1,2}]
                            [--reditools-extra REDITOOLS_EXTRA]
                            [--reditools-post-rscript REDITOOLS_POST_RSCRIPT]
                            [--reditools-min-meanq REDITOOLS_MIN_MEANQ]
                            [--reditools-min-coverage REDITOOLS_MIN_COVERAGE]
                            [--reditools-min-frequency REDITOOLS_MIN_FREQUENCY]
                            [--run-sprint] [--sprint-exe SPRINT_EXE]
                            [--sprint-repeat-bed SPRINT_REPEAT_BED]
                            [--sprint-geta2i SPRINT_GETA2I]
                            [--sprint-strand-specific {auto,0,1}]
                            [--bwa-exe BWA_EXE] [--sprint-extra SPRINT_EXTRA]
                            [--sprint-auto-decompress]
                            [--no-sprint-auto-decompress]
                            [--precomputed-csv-in PRECOMPUTED_CSV_IN]
                            [--precomputed-rmats-dir PRECOMPUTED_RMATS_DIR]
                            [--precomputed-redit-dir PRECOMPUTED_REDIT_DIR]
                            [--precomputed-sprint-dir PRECOMPUTED_SPRINT_DIR]
                            [--analyze-subset {inverted,hairpin,allpairs}]
                            [--window-w WINDOW_W] [--arm-aware]
                            [--no-arm-aware] [--arm-pad ARM_PAD]
                            [--arm-min-cov ARM_MIN_COV]
                            [--min-selected-candidates MIN_SELECTED_CANDIDATES]
                            [--do-ddg] [--no-ddg] [--do-pf-interface]
                            [--do-null-z] [--null-n NULL_N]
                            [--null-seed NULL_SEED] [--do-intarna]
                            [--cofold-strong COFOLD_STRONG]
                            [--cofold-moderate COFOLD_MODERATE]
                            [--transcript-mapping-rscript TRANSCRIPT_MAPPING_RSCRIPT]
                            [--priority-top-n PRIORITY_TOP_N]
                            [--priority-mode {strict,relaxed}]
                            [--no-require-case-editing] [--no-require-case-ri]
                            [--priority-score-mode {expert,adaptive,balanced,supervised}]
                            [--annotation-policy {conservative,zrna_permissive}]
                            [--training-truth-table TRAINING_TRUTH_TABLE]
                            [--training-labels TRAINING_LABELS]
                            [--truth-symbol-col TRUTH_SYMBOL_COL]
                            [--truth-label-mode {positive_logfc_padj,padj_only,all_table_rows,explicit_label_col}]
                            [--truth-label-col TRUTH_LABEL_COL]
                            [--truth-padj-col TRUTH_PADJ_COL]
                            [--truth-logfc-col TRUTH_LOGFC_COL]
                            [--truth-padj-max TRUTH_PADJ_MAX]
                            [--supervised-test-size SUPERVISED_TEST_SIZE]
                            [--cv-folds CV_FOLDS]
                            [--supervised-random-state SUPERVISED_RANDOM_STATE]
                            [--supervised-model {legacy_l2,elasticnet}]
                            [--supervised-feature-panel {raw7,routes6,orientation_only,annotation_only,orientation_annotation,structure_only,condition_only,structure_condition,no_orientation_annotation,annotation_structure_condition,compact_v1,all_current,auto_prespecified}]
                            [--supervised-tune] [--no-supervised-tune]
                            [--supervised-inner-cv-folds SUPERVISED_INNER_CV_FOLDS]
                            [--supervised-selection-metric {average_precision,roc_auc,balanced_accuracy}]
                            [--supervised-c SUPERVISED_C]
                            [--supervised-l1-ratio SUPERVISED_L1_RATIO]
                            [--supervised-c-grid SUPERVISED_C_GRID]
                            [--supervised-l1-ratio-grid SUPERVISED_L1_RATIO_GRID]
                            [--zrna-score-mode {pc1,sequence_pc1,consensus}]
                            [--zrna-class-mode {quantile,fixed}]
                            [--zrna-moderate-threshold ZRNA_MODERATE_THRESHOLD]
                            [--zrna-high-threshold ZRNA_HIGH_THRESHOLD]
                            [--star-exe STAR_EXE] [--python-exe PYTHON_EXE]
                            [--rscript-exe RSCRIPT_EXE]
                            [--bedtools-exe BEDTOOLS_EXE]
                            [--samtools-exe SAMTOOLS_EXE]
                            [--bamcoverage-exe BAMCOVERAGE_EXE]
                            [--multibigwigsummary-exe MULTIBIGWIGSUMMARY_EXE]
                            [--rnacofold-exe RNACOFOLD_EXE]
                            [--rnafold-exe RNAFOLD_EXE]
                            [--intarna-exe INTARNA_EXE]

options:
  -h, --help            show this help message and exit
  --output-dir OUTPUT_DIR
  --case-label CASE_LABEL
  --control-label CONTROL_LABEL
  --samplesheet SAMPLESHEET
                        FASTQ TSV/CSV with sample_id, condition, fastq_1,
                        fastq_2, strandedness; or BAM TSV with bam_path if
                        --input-mode bam
  --input-mode {fastq,bam}
  --fasta FASTA
  --gtf GTF
  --te-gtf TE_GTF       RepeatMasker/RMSK TE GTF; required only for --te-mode
                        simple. Advanced atena mode uses --te-genome/--te-
                        rmsk-rds.
  --strandedness {auto,forward,reverse,unstranded,fr-firststrand,fr-secondstrand}
                        Default auto: infer after alignment with RSeQC and
                        propagate to TE, rMATS, REDItools2, and SPRINT.
  --paired              Optional override. By default, infer paired/single
                        layout from the samplesheet and input files.
  --single-end          Optional override. By default, infer paired/single
                        layout automatically.
  --read-length READ_LENGTH
                        Optional override. By default, infer the maximum
                        observed read length from FASTQ/BAM input.
  --threads THREADS
  --force               Re-run steps even when expected outputs already exist
  --quiet               Redirect noisy tool stdout/stderr to
                        pipeline_info/logs (default)
  --verbose             Print internal STAR/rMATS/REDItools/SPRINT output to
                        terminal
  --infer-fastq-records INFER_FASTQ_RECORDS
                        FASTQ records inspected per mate/sample for layout and
                        read-length inference.
  --infer-experiment-exe INFER_EXPERIMENT_EXE
                        RSeQC strandedness inference executable.
  --strandedness-bed STRANDEDNESS_BED
                        Optional BED12 transcript annotation for RSeQC. If
                        absent, generated from --gtf.
  --strandedness-sample-size STRANDEDNESS_SAMPLE_SIZE
  --stranded-threshold STRANDED_THRESHOLD
  --unstranded-threshold UNSTRANDED_THRESHOLD
  --strandedness-fallback {unstranded,forward,reverse,error}
  --star-index STAR_INDEX
  --build-star-index    Legacy compatibility flag. A missing STAR index is now
                        built automatically.
  --sjdb-overhang SJDB_OVERHANG
                        Default: read_length - 1
  --skip-te-analysis    Skip TE analysis and reuse --precomputed-csv-in or
                        existing internal TE CSV
  --te-mode {advanced}  advanced = atena/qtex + DESeq2 + ChIPseeker (public
                        v1.x mode)
  --te-genome TE_GENOME
                        Genome key used by the advanced atena/ChIPseeker
                        module. Standard values include hg38, mm39, and mm10.
                        For a custom assembly such as C57BL_6J_T2T_v1, use
                        custom together with --te-rmsk-rds and --te-txdb-gtf.
  --te-rmsk-rds TE_RMSK_RDS
                        Optional cached atena RepeatMasker GRanges RDS. If
                        absent, built under 02_te/rmsk_<genome>_used.rds.
  --te-force-rebuild-rmsk
  --te-use-strand       Optional override. By default, strand-aware TE
                        counting follows inferred library strandedness.
  --te-ignore-strand
  --te-yield-size TE_YIELD_SIZE
  --te-min-max-count TE_MIN_MAX_COUNT
  --te-shrink-type {ashr,none}
  --te-txdb-package TE_TXDB_PACKAGE
                        Optional packaged TxDb override, e.g.
                        TxDb.Hsapiens.UCSC.hg38.knownGene. Do not use this for
                        a custom/T2T assembly when --te-txdb-gtf is supplied.
  --te-txdb-gtf TE_TXDB_GTF
                        Gene GTF used to build a custom TxDb for ChIPseeker.
                        Use this for assemblies without a packaged TxDb, such
                        as C57BL_6J_T2T_v1. The GTF must match --fasta and
                        --te-rmsk-rds coordinates.
  --te-txdb-rds TE_TXDB_RDS
                        Optional path for a cached custom TxDb RDS. If it
                        exists, the R module reuses it; otherwise it is built
                        from --te-txdb-gtf and saved here.
  --te-orgdb-package TE_ORGDB_PACKAGE
                        Optional OrgDb package override, e.g. org.Hs.eg.db or
                        org.Mm.eg.db. For mouse T2T, use org.Mm.eg.db.
  --te-feature-type TE_FEATURE_TYPE
                        Legacy internal fallback option; ignored by public
                        advanced TE mode
  --te-attribute TE_ATTRIBUTE
                        Legacy internal fallback option; ignored by public
                        advanced TE mode
  --te-padj-max TE_PADJ_MAX
  --te-lfc-min TE_LFC_MIN
                        For advanced mode, mirrors your old |log2FC| > 1
                        significant TE threshold by default.
  --te-candidate-mode {strict,expressed}
                        strict: dsRNASeeker pairs are built only from
                        padj/|LFC|-significant TE loci. expressed: pairs are
                        built from observed/expressed TE loci, while strict DE
                        status is kept as evidence columns. Use expressed for
                        RIP/dsRNA-seq validation datasets.
  --te-candidate-min-mean TE_CANDIDATE_MIN_MEAN
                        For --te-candidate-mode expressed, minimum max(mean
                        normalized count in case/control) for a TE locus to
                        enter the pair universe.
  --te-candidate-padj-max TE_CANDIDATE_PADJ_MAX
                        For --te-candidate-mode expressed, optional padj
                        cutoff for the candidate universe. 1.0 keeps all
                        expressed tested rows regardless of padj.
  --te-candidate-lfc-min TE_CANDIDATE_LFC_MIN
                        For --te-candidate-mode expressed, optional
                        abs(log2FC) cutoff for the candidate universe. 0 keeps
                        all expressed tested rows.
  --featurecounts-exe FEATURECOUNTS_EXE
                        Legacy internal fallback option; ignored by public
                        advanced TE mode
  --skip-rmats
  --rmats-exe RMATS_EXE
  --rmats-track {JC,JCEC}
  --rmats-fdr-max RMATS_FDR_MAX
  --rmats-overlap-slop RMATS_OVERLAP_SLOP
                        Pad TE-pair/arm intervals by this many bp when
                        intersecting rMATS RI events. 0 preserves exact-
                        overlap behavior.
  --rmats-interval-mode {intron_body,event_span}
                        For RI matching, use upstreamEE→downstreamES intron
                        body when available, or full
                        riExonStart_0base→riExonEnd event span.
  --rmats-cstat RMATS_CSTAT
  --rmats-min-intron-length RMATS_MIN_INTRON_LENGTH
                        Pass rMATS --mil. This only changes --novelSS event
                        detection. Default: 1, the most permissive positive
                        minimum intron length.
  --rmats-max-exon-length INT|auto
                        Pass rMATS --mel (maximum EXON length, not intron
                        length). Default: auto, which scans exon records in
                        --gtf and uses the largest annotated exon length
                        (never below the rMATS default 500). A positive
                        integer may be supplied to override the automatic
                        value.
  --rmats-libtype {None,fr-unstranded,fr-firststrand,fr-secondstrand}
  --rmats-novel-ss
  --no-rmats-novel-ss
  --skip-reditools
  --reditools-exe REDITOOLS_EXE
                        REDItools2 executable/script. Default assumes
                        src/cineca/reditools.py-style arguments
  --reditools-strand {auto,0,1,2}
                        Default auto: derive REDItools2 -s from inferred
                        library orientation.
  --reditools-extra REDITOOLS_EXTRA
                        Extra raw REDItools2 arguments, quoted as one string
  --reditools-post-rscript REDITOOLS_POST_RSCRIPT
                        Optional generic R postprocessor; default uses
                        r/reditools_filter_a2i.R
  --reditools-min-meanq REDITOOLS_MIN_MEANQ
  --reditools-min-coverage REDITOOLS_MIN_COVERAGE
  --reditools-min-frequency REDITOOLS_MIN_FREQUENCY
  --run-sprint
  --sprint-exe SPRINT_EXE
  --sprint-repeat-bed SPRINT_REPEAT_BED
                        RepeatMasker/repeat BED required by SPRINT -rp
  --sprint-geta2i SPRINT_GETA2I
                        Path to SPRINT/utilities/getA2I.py
  --sprint-strand-specific {auto,0,1}
                        Default auto: 0 for inferred unstranded data,
                        otherwise 1.
  --bwa-exe BWA_EXE
  --sprint-extra SPRINT_EXTRA
  --sprint-auto-decompress
                        Automatically materialize gzipped FASTQs for SPRINT
                        and reuse the decompressed cache.
  --no-sprint-auto-decompress
  --precomputed-csv-in PRECOMPUTED_CSV_IN
  --precomputed-rmats-dir PRECOMPUTED_RMATS_DIR
  --precomputed-redit-dir PRECOMPUTED_REDIT_DIR
  --precomputed-sprint-dir PRECOMPUTED_SPRINT_DIR
  --analyze-subset {inverted,hairpin,allpairs}
  --window-w WINDOW_W
  --arm-aware
  --no-arm-aware
  --arm-pad ARM_PAD
  --arm-min-cov ARM_MIN_COV
  --min-selected-candidates MIN_SELECTED_CANDIDATES
  --do-ddg
  --no-ddg
  --do-pf-interface
  --do-null-z
  --null-n NULL_N
  --null-seed NULL_SEED
  --do-intarna
  --cofold-strong COFOLD_STRONG
  --cofold-moderate COFOLD_MODERATE
  --transcript-mapping-rscript TRANSCRIPT_MAPPING_RSCRIPT
  --priority-top-n PRIORITY_TOP_N
  --priority-mode {strict,relaxed}
  --no-require-case-editing
  --no-require-case-ri
  --priority-score-mode {expert,adaptive,balanced,supervised}
  --annotation-policy {conservative,zrna_permissive}
  --training-truth-table TRAINING_TRUTH_TABLE
  --training-labels TRAINING_LABELS
  --truth-symbol-col TRUTH_SYMBOL_COL
  --truth-label-mode {positive_logfc_padj,padj_only,all_table_rows,explicit_label_col}
  --truth-label-col TRUTH_LABEL_COL
  --truth-padj-col TRUTH_PADJ_COL
  --truth-logfc-col TRUTH_LOGFC_COL
  --truth-padj-max TRUTH_PADJ_MAX
  --supervised-test-size SUPERVISED_TEST_SIZE
  --cv-folds CV_FOLDS
  --supervised-random-state SUPERVISED_RANDOM_STATE
  --supervised-model {legacy_l2,elasticnet}
  --supervised-feature-panel {raw7,routes6,orientation_only,annotation_only,orientation_annotation,structure_only,condition_only,structure_condition,no_orientation_annotation,annotation_structure_condition,compact_v1,all_current,auto_prespecified}
  --supervised-tune
  --no-supervised-tune
  --supervised-inner-cv-folds SUPERVISED_INNER_CV_FOLDS
  --supervised-selection-metric {average_precision,roc_auc,balanced_accuracy}
  --supervised-c SUPERVISED_C
  --supervised-l1-ratio SUPERVISED_L1_RATIO
  --supervised-c-grid SUPERVISED_C_GRID
  --supervised-l1-ratio-grid SUPERVISED_L1_RATIO_GRID
  --zrna-score-mode {pc1,sequence_pc1,consensus}
  --zrna-class-mode {quantile,fixed}
  --zrna-moderate-threshold ZRNA_MODERATE_THRESHOLD
  --zrna-high-threshold ZRNA_HIGH_THRESHOLD
  --star-exe STAR_EXE
  --python-exe PYTHON_EXE
  --rscript-exe RSCRIPT_EXE
  --bedtools-exe BEDTOOLS_EXE
  --samtools-exe SAMTOOLS_EXE
  --bamcoverage-exe BAMCOVERAGE_EXE
  --multibigwigsummary-exe MULTIBIGWIGSUMMARY_EXE
  --rnacofold-exe RNACOFOLD_EXE
  --rnafold-exe RNAFOLD_EXE
  --intarna-exe INTARNA_EXE
```

## `run`

```text
usage: dsRNASeeker run [-h] --output-dir OUTPUT_DIR --case-label CASE_LABEL
                       --control-label CONTROL_LABEL --condition CONDITION
                       --samplesheet SAMPLESHEET --csv-in CSV_IN --fasta FASTA
                       --gtf GTF [--sprint-a2i-dir SPRINT_A2I_DIR]
                       [--redit-dirs [REDIT_DIRS ...]]
                       [--analyze-subset {inverted,hairpin,allpairs}]
                       [--window-w WINDOW_W] [--arm-aware] [--no-arm-aware]
                       [--arm-pad ARM_PAD] [--arm-min-cov ARM_MIN_COV]
                       [--do-ddg] [--no-ddg] [--do-pf-interface] [--do-null-z]
                       [--null-n NULL_N] [--null-seed NULL_SEED]
                       [--do-intarna] [--cofold-strong COFOLD_STRONG]
                       [--cofold-moderate COFOLD_MODERATE]
                       [--transcript-mapping-rscript TRANSCRIPT_MAPPING_RSCRIPT]
                       [--python-exe PYTHON_EXE] [--rscript-exe RSCRIPT_EXE]
                       [--bedtools-exe BEDTOOLS_EXE]
                       [--samtools-exe SAMTOOLS_EXE]
                       [--bamcoverage-exe BAMCOVERAGE_EXE]
                       [--multibigwigsummary-exe MULTIBIGWIGSUMMARY_EXE]
                       [--rnacofold-exe RNACOFOLD_EXE]
                       [--rnafold-exe RNAFOLD_EXE] [--intarna-exe INTARNA_EXE]

options:
  -h, --help            show this help message and exit
  --output-dir OUTPUT_DIR
  --case-label CASE_LABEL
  --control-label CONTROL_LABEL
  --condition CONDITION
  --samplesheet SAMPLESHEET
                        TSV with columns: sample_id, condition, bam_path
  --csv-in CSV_IN
  --fasta FASTA
  --gtf GTF
  --sprint-a2i-dir SPRINT_A2I_DIR
  --redit-dirs [REDIT_DIRS ...]
  --analyze-subset {inverted,hairpin,allpairs}
  --window-w WINDOW_W
  --arm-aware
  --no-arm-aware
  --arm-pad ARM_PAD
  --arm-min-cov ARM_MIN_COV
  --do-ddg
  --no-ddg
  --do-pf-interface
  --do-null-z
  --null-n NULL_N
  --null-seed NULL_SEED
  --do-intarna
  --cofold-strong COFOLD_STRONG
  --cofold-moderate COFOLD_MODERATE
  --transcript-mapping-rscript TRANSCRIPT_MAPPING_RSCRIPT
  --python-exe PYTHON_EXE
  --rscript-exe RSCRIPT_EXE
  --bedtools-exe BEDTOOLS_EXE
  --samtools-exe SAMTOOLS_EXE
  --bamcoverage-exe BAMCOVERAGE_EXE
  --multibigwigsummary-exe MULTIBIGWIGSUMMARY_EXE
  --rnacofold-exe RNACOFOLD_EXE
  --rnafold-exe RNAFOLD_EXE
  --intarna-exe INTARNA_EXE
```

## `summary`

```text
usage: dsRNASeeker summary [-h] --output-dir OUTPUT_DIR --case-label
                           CASE_LABEL --control-label CONTROL_LABEL --csv-in
                           CSV_IN
                           [--analyze-subset {inverted,hairpin,allpairs}]
                           [--rmats-dir RMATS_DIR] [--rmats-track {JC,JCEC}]
                           [--rmats-fdr-max RMATS_FDR_MAX]
                           [--rmats-overlap-slop RMATS_OVERLAP_SLOP]
                           [--rmats-interval-mode {intron_body,event_span}]
                           [--rmats-group1-label RMATS_GROUP1_LABEL]
                           [--rmats-group2-label RMATS_GROUP2_LABEL]
                           [--rmats-flip-dpsi] [--bedtools-exe BEDTOOLS_EXE]
                           [--priority-top-n PRIORITY_TOP_N]
                           [--priority-mode {strict,relaxed}]
                           [--no-require-case-editing] [--no-require-case-ri]
                           [--priority-score-mode {expert,adaptive,balanced,supervised}]
                           [--annotation-policy {conservative,zrna_permissive}]
                           [--training-truth-table TRAINING_TRUTH_TABLE]
                           [--training-labels TRAINING_LABELS]
                           [--truth-symbol-col TRUTH_SYMBOL_COL]
                           [--truth-label-mode {positive_logfc_padj,padj_only,all_table_rows,explicit_label_col}]
                           [--truth-label-col TRUTH_LABEL_COL]
                           [--truth-padj-col TRUTH_PADJ_COL]
                           [--truth-logfc-col TRUTH_LOGFC_COL]
                           [--truth-padj-max TRUTH_PADJ_MAX]
                           [--supervised-test-size SUPERVISED_TEST_SIZE]
                           [--cv-folds CV_FOLDS]
                           [--supervised-random-state SUPERVISED_RANDOM_STATE]
                           [--supervised-model {legacy_l2,elasticnet}]
                           [--supervised-feature-panel {raw7,routes6,orientation_only,annotation_only,orientation_annotation,structure_only,condition_only,structure_condition,no_orientation_annotation,annotation_structure_condition,compact_v1,all_current,auto_prespecified}]
                           [--supervised-tune] [--no-supervised-tune]
                           [--supervised-inner-cv-folds SUPERVISED_INNER_CV_FOLDS]
                           [--supervised-selection-metric {average_precision,roc_auc,balanced_accuracy}]
                           [--supervised-c SUPERVISED_C]
                           [--supervised-l1-ratio SUPERVISED_L1_RATIO]
                           [--supervised-c-grid SUPERVISED_C_GRID]
                           [--supervised-l1-ratio-grid SUPERVISED_L1_RATIO_GRID]

options:
  -h, --help            show this help message and exit
  --output-dir OUTPUT_DIR
  --case-label CASE_LABEL
  --control-label CONTROL_LABEL
  --csv-in CSV_IN
  --analyze-subset {inverted,hairpin,allpairs}
  --rmats-dir RMATS_DIR
  --rmats-track {JC,JCEC}
  --rmats-fdr-max RMATS_FDR_MAX
  --rmats-overlap-slop RMATS_OVERLAP_SLOP
                        Pad TE-pair/arm intervals by this many bp when
                        intersecting rMATS RI events. 0 preserves exact-
                        overlap behavior.
  --rmats-interval-mode {intron_body,event_span}
                        For RI matching, use upstreamEE→downstreamES intron
                        body when available, or full
                        riExonStart_0base→riExonEnd event span.
  --rmats-group1-label RMATS_GROUP1_LABEL
  --rmats-group2-label RMATS_GROUP2_LABEL
  --rmats-flip-dpsi     Flip rMATS IncLevelDifference sign before RI direction
                        calls. Use when rMATS was run as control-minus-case
                        but summary should report case-minus-control.
  --bedtools-exe BEDTOOLS_EXE
  --priority-top-n PRIORITY_TOP_N
                        Number of strict high-priority candidates to export
                        separately.
  --priority-mode {strict,relaxed}
                        strict requires case editing and case RI; relaxed
                        keeps gates but labels incomplete candidates.
  --no-require-case-editing
                        Do not require case-enriched SPRINT/REDI editing for
                        strict priority_gate_pass.
  --no-require-case-ri  Do not require case-high rMATS RI for strict
                        priority_gate_pass.
  --priority-score-mode {expert,adaptive,balanced,supervised}
                        adaptive uses gate-derived ADPS weights; balanced uses
                        equal weights across non-gate evidence blocks;
                        supervised uses label-driven logistic ranking.
  --annotation-policy {conservative,zrna_permissive}
                        zrna_permissive differs only by allowing known
                        3UTR-3UTR pairs; unknown annotations remain rejected.
  --training-truth-table TRAINING_TRUTH_TABLE
                        Gene-level truth table used to derive pair labels by
                        matching truth symbols to A_SYMBOL/B_SYMBOL.
  --training-labels TRAINING_LABELS
                        Optional precomputed pair-level labels file with
                        columns pair_id and label.
  --truth-symbol-col TRUTH_SYMBOL_COL
                        Gene-symbol column in --training-truth-table.
  --truth-label-mode {positive_logfc_padj,padj_only,all_table_rows,explicit_label_col}
                        How positives are defined from --training-truth-table.
  --truth-label-col TRUTH_LABEL_COL
                        0/1 label column used when --truth-label-mode
                        explicit_label_col.
  --truth-padj-col TRUTH_PADJ_COL
                        Adjusted-P/FDR column used by positive_logfc_padj or
                        padj_only modes.
  --truth-logfc-col TRUTH_LOGFC_COL
                        Log2 fold-change column used by positive_logfc_padj
                        mode.
  --truth-padj-max TRUTH_PADJ_MAX
                        Adjusted-P/FDR threshold used by positive_logfc_padj
                        or padj_only modes.
  --supervised-test-size SUPERVISED_TEST_SIZE
                        Held-out test fraction for supervised diagnostics. Use
                        0 to disable held-out testing.
  --cv-folds CV_FOLDS   Number of stratified CV folds for supervised
                        diagnostics. Use 0 to disable CV.
  --supervised-random-state SUPERVISED_RANDOM_STATE
                        Random seed for supervised train/test and cross-
                        validation splits.
  --supervised-model {legacy_l2,elasticnet}
                        legacy_l2 is a fixed L2 comparator; elasticnet permits
                        L2/L1/Elastic-Net selection, using LBFGS for
                        l1_ratio=0 and SAGA otherwise.
  --supervised-feature-panel {raw7,routes6,orientation_only,annotation_only,orientation_annotation,structure_only,condition_only,structure_condition,no_orientation_annotation,annotation_structure_condition,compact_v1,all_current,auto_prespecified}
                        Use raw7 evidence components, routes6 pre-specified
                        composites, a backward-compatible panel, or select
                        raw7 versus routes6 inside grouped training CV.
  --supervised-tune     Tune regularization, and panel when auto_prespecified,
                        using grouped CV within the training matrix only.
  --no-supervised-tune
  --supervised-inner-cv-folds SUPERVISED_INNER_CV_FOLDS
  --supervised-selection-metric {average_precision,roc_auc,balanced_accuracy}
  --supervised-c SUPERVISED_C
                        Fixed inverse regularization strength when
                        --supervised-tune is not used.
  --supervised-l1-ratio SUPERVISED_L1_RATIO
                        Fixed elastic-net L1 ratio when --supervised-tune is
                        not used.
  --supervised-c-grid SUPERVISED_C_GRID
  --supervised-l1-ratio-grid SUPERVISED_L1_RATIO_GRID
```

## `delta`

```text
usage: dsRNASeeker delta [-h] --output-dir OUTPUT_DIR --case-label CASE_LABEL
                         --control-label CONTROL_LABEL
                         [--analyze-subset {inverted,hairpin,allpairs}]

options:
  -h, --help            show this help message and exit
  --output-dir OUTPUT_DIR
  --case-label CASE_LABEL
  --control-label CONTROL_LABEL
  --analyze-subset {inverted,hairpin,allpairs}
```

## `zrna`

```text
usage: dsRNASeeker zrna [-h] --output-dir OUTPUT_DIR --case-label CASE_LABEL
                        --control-label CONTROL_LABEL
                        [--analyze-subset {inverted,hairpin,allpairs}]
                        [--summary-in SUMMARY_IN] [--case-fasta CASE_FASTA]
                        [--control-fasta CONTROL_FASTA]
                        [--zrna-score-mode {pc1,sequence_pc1,consensus}]
                        [--zrna-class-mode {quantile,fixed}]
                        [--zrna-moderate-threshold ZRNA_MODERATE_THRESHOLD]
                        [--zrna-high-threshold ZRNA_HIGH_THRESHOLD]

options:
  -h, --help            show this help message and exit
  --output-dir OUTPUT_DIR
  --case-label CASE_LABEL
  --control-label CONTROL_LABEL
  --analyze-subset {inverted,hairpin,allpairs}
  --summary-in SUMMARY_IN
  --case-fasta CASE_FASTA
  --control-fasta CONTROL_FASTA
  --zrna-score-mode {pc1,sequence_pc1,consensus}
  --zrna-class-mode {quantile,fixed}
  --zrna-moderate-threshold ZRNA_MODERATE_THRESHOLD
  --zrna-high-threshold ZRNA_HIGH_THRESHOLD
```

## `molecule-model`

```text
usage: dsRNASeeker molecule-model [-h] --output-dir OUTPUT_DIR --case-label
                                  CASE_LABEL --control-label CONTROL_LABEL
                                  --gtf GTF
                                  [--analyze-subset {inverted,hairpin,allpairs}]
                                  [--summary-in SUMMARY_IN]
                                  [--zrna-summary ZRNA_SUMMARY]
                                  [--include-default-zrna]
                                  [--no-include-default-zrna]
                                  [--transcript-overlap-fraction TRANSCRIPT_OVERLAP_FRACTION]
                                  [--transcript-containment-slop TRANSCRIPT_CONTAINMENT_SLOP]
                                  [--max-transcript-ids MAX_TRANSCRIPT_IDS]

options:
  -h, --help            show this help message and exit
  --output-dir OUTPUT_DIR
  --case-label CASE_LABEL
  --control-label CONTROL_LABEL
  --gtf GTF
  --analyze-subset {inverted,hairpin,allpairs}
  --summary-in SUMMARY_IN
  --zrna-summary ZRNA_SUMMARY
  --include-default-zrna
  --no-include-default-zrna
  --transcript-overlap-fraction TRANSCRIPT_OVERLAP_FRACTION
  --transcript-containment-slop TRANSCRIPT_CONTAINMENT_SLOP
  --max-transcript-ids MAX_TRANSCRIPT_IDS
```

## `robustness`

```text
usage: dsRNASeeker robustness [-h] --summary-in SUMMARY_IN
                              [--output-dir OUTPUT_DIR]
                              [--case-label CASE_LABEL]
                              [--dominant-weight-threshold DOMINANT_WEIGHT_THRESHOLD]
                              [--agreement-tight-iqr AGREEMENT_TIGHT_IQR]
                              [--agreement-wide-iqr AGREEMENT_WIDE_IQR]

options:
  -h, --help            show this help message and exit
  --summary-in SUMMARY_IN
                        Existing TEpair_dsRNA_master.summary.with_RI.csv. No
                        upstream stages are rerun.
  --output-dir OUTPUT_DIR
                        Output directory for diagnostics. Default: directory
                        containing --summary-in.
  --case-label CASE_LABEL
                        Case label used in output filenames. If omitted, infer
                        from nearby adaptive-weight files when possible.
  --dominant-weight-threshold DOMINANT_WEIGHT_THRESHOLD
                        Descriptive threshold for flagging a single-component-
                        dominated ADPS solution; not an accuracy threshold.
  --agreement-tight-iqr AGREEMENT_TIGHT_IQR
                        Descriptive broad-profile percentile-IQR cutoff for
                        tight rank agreement.
  --agreement-wide-iqr AGREEMENT_WIDE_IQR
                        Descriptive broad-profile percentile-IQR cutoff above
                        which rank agreement is called wide/profile-sensitive.
```

## `supervised-benchmark`

```text
usage: dsRNASeeker supervised-benchmark [-h] --manifest MANIFEST
                                        [--input-root INPUT_ROOT] --output-dir
                                        OUTPUT_DIR [--outer-folds OUTER_FOLDS]
                                        [--inner-folds INNER_FOLDS]
                                        [--selection-metric {average_precision,roc_auc,balanced_accuracy}]
                                        [--c-grid C_GRID]
                                        [--l1-ratio-grid L1_RATIO_GRID]
                                        [--same-study-min-positive SAME_STUDY_MIN_POSITIVE]
                                        [--same-study-min-negative SAME_STUDY_MIN_NEGATIVE]
                                        [--bootstrap BOOTSTRAP] [--seed SEED]

options:
  -h, --help            show this help message and exit
  --manifest MANIFEST   TSV/CSV with dataset_id, analysis_variant,
                        study_family, audit_table, same_study_eligible,
                        loso_representative, and loso_target.
  --input-root INPUT_ROOT
                        Base directory for relative audit_table paths in the
                        manifest.
  --output-dir OUTPUT_DIR
  --outer-folds OUTER_FOLDS
  --inner-folds INNER_FOLDS
  --selection-metric {average_precision,roc_auc,balanced_accuracy}
  --c-grid C_GRID
  --l1-ratio-grid L1_RATIO_GRID
  --same-study-min-positive SAME_STUDY_MIN_POSITIVE
  --same-study-min-negative SAME_STUDY_MIN_NEGATIVE
  --bootstrap BOOTSTRAP
  --seed SEED
```

## `check`

```text
usage: dsRNASeeker check [-h] --output-dir OUTPUT_DIR --case-label CASE_LABEL
                         --control-label CONTROL_LABEL --condition CONDITION
                         --samplesheet SAMPLESHEET --csv-in CSV_IN --fasta
                         FASTA --gtf GTF [--sprint-a2i-dir SPRINT_A2I_DIR]
                         [--redit-dirs [REDIT_DIRS ...]]
                         [--analyze-subset {inverted,hairpin,allpairs}]
                         [--window-w WINDOW_W] [--arm-aware] [--no-arm-aware]
                         [--arm-pad ARM_PAD] [--arm-min-cov ARM_MIN_COV]
                         [--do-ddg] [--no-ddg] [--do-pf-interface]
                         [--do-null-z] [--null-n NULL_N]
                         [--null-seed NULL_SEED] [--do-intarna]
                         [--cofold-strong COFOLD_STRONG]
                         [--cofold-moderate COFOLD_MODERATE]
                         [--transcript-mapping-rscript TRANSCRIPT_MAPPING_RSCRIPT]
                         [--python-exe PYTHON_EXE] [--rscript-exe RSCRIPT_EXE]
                         [--bedtools-exe BEDTOOLS_EXE]
                         [--samtools-exe SAMTOOLS_EXE]
                         [--bamcoverage-exe BAMCOVERAGE_EXE]
                         [--multibigwigsummary-exe MULTIBIGWIGSUMMARY_EXE]
                         [--rnacofold-exe RNACOFOLD_EXE]
                         [--rnafold-exe RNAFOLD_EXE]
                         [--intarna-exe INTARNA_EXE]

options:
  -h, --help            show this help message and exit
  --output-dir OUTPUT_DIR
  --case-label CASE_LABEL
  --control-label CONTROL_LABEL
  --condition CONDITION
  --samplesheet SAMPLESHEET
                        TSV with columns: sample_id, condition, bam_path
  --csv-in CSV_IN
  --fasta FASTA
  --gtf GTF
  --sprint-a2i-dir SPRINT_A2I_DIR
  --redit-dirs [REDIT_DIRS ...]
  --analyze-subset {inverted,hairpin,allpairs}
  --window-w WINDOW_W
  --arm-aware
  --no-arm-aware
  --arm-pad ARM_PAD
  --arm-min-cov ARM_MIN_COV
  --do-ddg
  --no-ddg
  --do-pf-interface
  --do-null-z
  --null-n NULL_N
  --null-seed NULL_SEED
  --do-intarna
  --cofold-strong COFOLD_STRONG
  --cofold-moderate COFOLD_MODERATE
  --transcript-mapping-rscript TRANSCRIPT_MAPPING_RSCRIPT
  --python-exe PYTHON_EXE
  --rscript-exe RSCRIPT_EXE
  --bedtools-exe BEDTOOLS_EXE
  --samtools-exe SAMTOOLS_EXE
  --bamcoverage-exe BAMCOVERAGE_EXE
  --multibigwigsummary-exe MULTIBIGWIGSUMMARY_EXE
  --rnacofold-exe RNACOFOLD_EXE
  --rnafold-exe RNAFOLD_EXE
  --intarna-exe INTARNA_EXE
```
