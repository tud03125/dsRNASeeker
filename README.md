# dsRNASeeker

**dsRNASeeker** is a condition-aware workflow for discovering and prioritizing putative **inverted transposable-element (TE) pairs** that may contribute to endogenous double-stranded RNA (dsRNA) formation from conventional short-read RNA-seq.

> **Scope.** dsRNASeeker targets an important inverted-TE-associated subclass of dsRNA. It is not a universal detector of every biological dsRNA source, and its structural outputs are computational predictions rather than direct proof of in-vivo duplex formation.

## What dsRNASeeker integrates

The end-to-end workflow can combine:

1. read alignment and library/strandedness inference;
2. locus-specific TE expression and annotation;
3. retained-intron/splicing evidence;
4. A-to-I RNA-editing evidence;
5. nearby inverted TE-pair discovery;
6. duplex energetics and pairing-interface support;
7. optional null-model, IntaRNA, and Z-RNA annotations; and
8. evidence integration with the **Adaptive dsRNA Prioritization Score (ADPS)**.

ADPS is a continuous, label-independent candidate-ranking score. Conservative priority tiers are reported separately. Optional supervised reranking is available when independent candidate labels are supplied.

## Repository layout

```text
dsRNASeeker/
├── main.py
├── modules/                 # core Python implementation
├── r/                       # R helpers for TE annotation/editing/transcript mapping
├── wrappers/                # portable wrappers for bundled optional tools
├── tools/
│   ├── REDItools2/          # Git submodule
│   └── SPRINT/              # Git submodule; optional
├── environment.yml          # authoritative Conda environment
├── samplesheet.example.tsv
├── LICENSE
└── README.md
```

The public release should contain one README (`README.md`) and one Conda environment definition (`environment.yml`).

## Platform

The manuscript implementation was developed and tested on 64-bit Linux/HPC systems. The environment uses packages from `conda-forge` and `bioconda`.

## Installation

### 1. Clone the repository and bundled submodules

```bash
git clone --recurse-submodules https://github.com/tud03125/dsRNASeeker.git
cd dsRNASeeker
```

If the repository was already cloned without submodules:

```bash
git submodule update --init --recursive
```

### 2. Configure Conda/Bioconda channels

```bash
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```

### 3. Create and activate the environment

```bash
conda env create -f environment.yml
conda activate dsRNASeeker
```

### 4. Make bundled wrappers executable

```bash
chmod +x wrappers/*.sh
```

### 5. Verify the installation

```bash
python3 main.py --help
python3 main.py workflow --help
```

Optional quick dependency checks:

```bash
python3 - <<'PY'
import numpy, pandas, sklearn, joblib, pysam, pyBigWig, RNA
print("Python imports: OK")
PY

for x in STAR samtools bedtools featureCounts rmats.py infer_experiment.py \
         bamCoverage multiBigwigSummary RNAfold RNAcofold IntaRNA bwa; do
    command -v "$x" >/dev/null && echo "[OK] $x" || echo "[MISSING] $x"
done
```

## Input samplesheets

### FASTQ input

Tab-delimited or CSV:

```text
sample_id	condition	fastq_1	fastq_2
CASE_1	CASE	/path/CASE_1_R1.fastq.gz	/path/CASE_1_R2.fastq.gz
CASE_2	CASE	/path/CASE_2_R1.fastq.gz	/path/CASE_2_R2.fastq.gz
CTRL_1	CONTROL	/path/CTRL_1_R1.fastq.gz	/path/CTRL_1_R2.fastq.gz
CTRL_2	CONTROL	/path/CTRL_2_R1.fastq.gz	/path/CTRL_2_R2.fastq.gz
```

For single-end FASTQ, leave `fastq_2` empty. Paired/single-end layout and read length can be inferred automatically unless explicitly overridden.

### BAM input

```text
sample_id	condition	bam_path
CASE_1	CASE	/path/CASE_1.bam
CASE_2	CASE	/path/CASE_2.bam
CTRL_1	CONTROL	/path/CTRL_1.bam
CTRL_2	CONTROL	/path/CTRL_2.bam
```

## Quick start: end-to-end workflow

The `workflow` command is the recommended public interface.

### Standard genome, FASTQ input

Example for hg38:

```bash
python3 main.py workflow \
  --output-dir results/example_hg38 \
  --case-label CASE \
  --control-label CONTROL \
  --samplesheet samplesheet.fastq.tsv \
  --input-mode fastq \
  --fasta /path/hg38.fa \
  --gtf /path/hg38.gtf \
  --te-mode advanced \
  --te-genome hg38 \
  --te-candidate-mode strict \
  --analyze-subset inverted \
  --window-w 1000 \
  --do-pf-interface \
  --do-null-z \
  --do-intarna \
  --reditools-exe "$PWD/wrappers/reditools2.sh" \
  --threads 16
```

For standard mouse analyses, use `--te-genome mm10` or `--te-genome mm39` with matching FASTA/GTF/BAM coordinates.

### BAM input

```bash
python3 main.py workflow \
  --output-dir results/example_bam \
  --case-label CASE \
  --control-label CONTROL \
  --samplesheet samplesheet.bam.tsv \
  --input-mode bam \
  --fasta /path/reference.fa \
  --gtf /path/genes.gtf \
  --te-mode advanced \
  --te-genome mm39 \
  --te-candidate-mode strict \
  --analyze-subset inverted \
  --window-w 1000 \
  --do-pf-interface \
  --reditools-exe "$PWD/wrappers/reditools2.sh" \
  --threads 16
```

### Validation-oriented / broader candidate universe

For RIP/J2/Z22/FLAG-style validation analyses where the candidate universe should not be restricted to differentially expressed TE loci:

```bash
--te-candidate-mode expressed
```

The strict/expressed choice changes the candidate universe; it is not an experimental ground-truth label.

### Custom or T2T assembly

For an assembly without a packaged TxDb, supply assembly-matched repeat and transcript annotations:

```bash
python3 main.py workflow \
  ... \
  --te-mode advanced \
  --te-genome custom \
  --te-rmsk-rds /path/custom_repeatmasker_granges.rds \
  --te-txdb-gtf /path/custom_genes.gtf
```

The FASTA, GTF, BAM coordinates, RepeatMasker RDS, and custom TxDb GTF must use the same assembly/coordinate system.

## Reusing precomputed components

The workflow can reuse existing products instead of regenerating them:

```text
--precomputed-csv-in
--precomputed-rmats-dir
--precomputed-redit-dir
--precomputed-sprint-dir
```

Corresponding skip flags include:

```text
--skip-te-analysis
--skip-rmats
--skip-reditools
```

This is useful when integrating pre-existing TE, splicing, or editing analyses.

## Main commands

| Command | Purpose |
|---|---|
| `workflow` | Recommended end-to-end FASTQ/BAM workflow |
| `run` | Run one condition through the core dsRNASeeker modules |
| `summary` | Fuse case/control evidence and calculate prioritization outputs |
| `delta` | Build case-control delta tables |
| `zrna` | Add A-form/Z-RNA propensity annotations |
| `molecule-model` | Add conservative intramolecular/intermolecular compatibility annotations |
| `robustness` | Label-free ranking-robustness diagnostics from an existing summary |
| `supervised-benchmark` | Nested grouped and leave-one-study-family-out supervised evaluation |
| `check` | Check core runtime dependencies and input files |

Use the CLI help as the authoritative parameter reference:

```bash
python3 main.py <command> --help
```

## Common workflow options

### Candidate discovery

```text
--analyze-subset {inverted,hairpin,allpairs}
--window-w INT
--arm-aware / --no-arm-aware
--arm-pad INT
--arm-min-cov FLOAT
```

The manuscript uses inverted TE pairs as the primary target class. A 1-kb search radius is the manuscript default; it should be interpreted as a parsimonious computational setting, not a universal biological cutoff.

### TE analysis

```text
--te-mode advanced
--te-genome {hg38,mm10,mm39,custom,...}
--te-candidate-mode {strict,expressed}
--te-padj-max FLOAT
--te-lfc-min FLOAT
```

### Structure support

```text
--do-ddg / --no-ddg
--do-pf-interface
--do-null-z
--null-n INT
--do-intarna
```

### RNA editing and splicing

```text
--skip-rmats
--rmats-track {JC,JCEC}
--rmats-fdr-max FLOAT
--skip-reditools
--reditools-exe PATH
--run-sprint
```

REDItools2 is the default editing route in the end-to-end workflow. SPRINT is optional and may require installation/configuration specific to the local system; use `--sprint-exe`, `--sprint-geta2i`, and related options when enabling it.

### Prioritization

```text
--priority-score-mode {expert,adaptive,balanced,supervised}
--priority-mode {strict,relaxed}
--priority-top-n INT
```

The manuscript-primary canonical workflow uses label-independent adaptive prioritization. Supervised reranking is optional and should not be interpreted as the default ADPS.

## Output structure

An end-to-end run creates a structured output directory containing, as applicable:

```text
00_reference/
01_alignment/
02_te/
03_splicing/
04_editing/
pipeline_info/
```

The core dsRNASeeker stages additionally write condition-specific candidate/evidence outputs, fused case-control summaries, delta tables, ranking/priority columns, and Z-RNA annotations.

`pipeline_info/software_paths.tsv` and `pipeline_info/dsRNASeeker_params.json` record resolved software paths and run parameters for reproducibility.

## Interpretation

dsRNASeeker prioritizes **putative** inverted-TE-pair dsRNA candidates.

- Genomic proximity and opposite orientation do not prove that both TE arms occur on the same RNA molecule.
- Predicted energetics and pairing-interface support do not prove an in-vivo duplex.
- J2, Z22, FLAG-ZBP1, A-to-I editing, and related assays identify overlapping but non-identical dsRNA-associated target classes.
- ADPS is a ranking score, not a calibrated probability of physical dsRNA formation.
- Priority tiers are conservative categorical outputs and are separate from continuous-score ROC/PR benchmarking.

## Reproducibility and release freezing

For a manuscript release, record both the Git tag and exact commit SHA. The human-maintained `environment.yml` should remain the installation source of truth.

For an exact Linux reproduction snapshot of the tested Conda environment, create an explicit specification after activating the release environment:

```bash
conda list --explicit > environment-linux-64-explicit.txt
```

This platform-specific file can be attached to the GitHub release or archived with the manuscript reproducibility materials.

## Citation

A manuscript citation and permanent software DOI will be added to the tagged manuscript release. Until then, cite the repository and the exact commit used.

Repository: https://github.com/tud03125/dsRNASeeker

## License

See `LICENSE`.
