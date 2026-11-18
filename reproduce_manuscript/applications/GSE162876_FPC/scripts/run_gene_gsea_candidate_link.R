#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(fgsea)
  library(msigdbr)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx)) return(default)
  if (idx == length(args)) stop("Missing value after ", flag)
  args[[idx + 1]]
}
required <- c("--counts", "--samples", "--summary", "--case", "--control", "--species", "--dataset", "--outdir")
missing <- required[vapply(required, function(x) is.null(get_arg(x)), logical(1))]
if (length(missing)) stop("Missing arguments: ", paste(missing, collapse = ", "))

counts_file <- get_arg("--counts")
samples_file <- get_arg("--samples")
summary_file <- get_arg("--summary")
case_label <- get_arg("--case")
control_label <- get_arg("--control")
species_key <- tolower(get_arg("--species"))
dataset_label <- get_arg("--dataset")
outdir <- get_arg("--outdir")
min_count <- as.integer(get_arg("--min-count", "10"))
min_size <- as.integer(get_arg("--min-size", "15"))
max_size <- as.integer(get_arg("--max-size", "500"))
top_n <- as.integer(get_arg("--top-n", "15"))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

if (!species_key %in% c("human", "mouse")) stop("--species must be human or mouse")
msig_species <- if (species_key == "human") "Homo sapiens" else "Mus musculus"
db_species <- if (species_key == "human") "HS" else "MM"

message("[INFO] Reading gene counts and sample metadata")
counts_df <- fread(counts_file, data.table = FALSE, check.names = FALSE)
if (!"gene_symbol" %in% names(counts_df)) stop("Counts table must contain gene_symbol")
samples <- fread(samples_file, data.table = FALSE)
samples <- samples %>% filter(condition %in% c(control_label, case_label))
if (!all(c(control_label, case_label) %in% samples$condition)) stop("Both case and control are required")
rownames(samples) <- samples$sample_id
sample_ids <- samples$sample_id
missing_cols <- setdiff(sample_ids, names(counts_df))
if (length(missing_cols)) stop("Counts missing samples: ", paste(missing_cols, collapse = ", "))

cts <- as.matrix(counts_df[, sample_ids, drop = FALSE])
storage.mode(cts) <- "integer"
rownames(cts) <- make.unique(as.character(counts_df$gene_symbol))
samples$condition <- factor(samples$condition, levels = c(control_label, case_label))

# For bulk RNA-seq, remove only extremely low-count genes. Do not preselect
# differentially expressed genes before GSEA.
smallest_group <- min(table(samples$condition))
keep <- rowSums(cts >= min_count) >= smallest_group
cts <- cts[keep, , drop = FALSE]
message("[INFO] Genes retained after low-count filter: ", nrow(cts))

dds <- DESeqDataSetFromMatrix(countData = cts, colData = samples, design = ~ condition)
dds <- DESeq(dds, quiet = TRUE)
res <- results(dds, contrast = c("condition", case_label, control_label), independentFiltering = TRUE)
res_df <- as.data.frame(res) %>%
  tibble::rownames_to_column("gene_symbol") %>%
  mutate(dataset = dataset_label, case = case_label, control = control_label)
write_csv(res_df, file.path(outdir, "gene_DESeq2_results.csv"))

# Use the unshrunk DESeq2 Wald statistic as the pre-ranked metric. This uses
# effect direction and uncertainty and avoids filtering genes by padj.
rank_df <- res_df %>%
  filter(!is.na(stat), is.finite(stat), !is.na(gene_symbol), gene_symbol != "") %>%
  group_by(gene_symbol) %>%
  slice_max(order_by = abs(stat), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(gene_symbol) %>%
  mutate(stat_for_gsea = stat + seq_along(stat) * 1e-12) %>%
  arrange(desc(stat_for_gsea))
gene_ranks <- rank_df$stat_for_gsea
names(gene_ranks) <- rank_df$gene_symbol
write_tsv(rank_df %>% select(gene_symbol, stat_for_gsea), file.path(outdir, "gene_ranked_list.rnk"), col_names = FALSE)

message("[INFO] Reading dsRNASeeker candidate-associated genes")
summary_df <- fread(summary_file, data.table = FALSE, check.names = FALSE)
extract_symbols <- function(df) {
  cols <- intersect(c("A_SYMBOL", "B_SYMBOL", "A_SYMBOL_from_TE", "B_SYMBOL_from_TE"), names(df))
  values <- unlist(df[, cols, drop = FALSE], use.names = FALSE)
  values <- values[!is.na(values)]
  tokens <- unlist(strsplit(as.character(values), "[;,|]"), use.names = FALSE)
  unique(trimws(tokens[trimws(tokens) != "" & !tolower(trimws(tokens)) %in% c("na", "nan", "none", ".")]))
}
all_candidate_genes <- extract_symbols(summary_df)
strict_candidate_genes <- character(0)
if ("priority_gate_pass" %in% names(summary_df)) {
  gate <- tolower(as.character(summary_df$priority_gate_pass)) %in% c("true", "t", "1", "yes", "y")
  strict_candidate_genes <- extract_symbols(summary_df[gate, , drop = FALSE])
}
top20_candidate_genes <- character(0)
if ("priority_rank" %in% names(summary_df)) {
  top_rows <- summary_df %>% mutate(.rank = suppressWarnings(as.numeric(priority_rank))) %>% arrange(.rank) %>% slice_head(n = 20)
  top20_candidate_genes <- extract_symbols(top_rows)
}
structure_candidate_genes <- character(0)
if (all(c("energy_adps", "interface_adps") %in% names(summary_df))) {
  structure_rows <- summary_df %>%
    mutate(.structure = rowMeans(cbind(suppressWarnings(as.numeric(energy_adps)), suppressWarnings(as.numeric(interface_adps))), na.rm = TRUE)) %>%
    arrange(desc(.structure)) %>% slice_head(n = 20)
  structure_candidate_genes <- extract_symbols(structure_rows)
}

candidate_sets <- list(
  dsRNASeeker_all_candidate_genes = intersect(all_candidate_genes, names(gene_ranks)),
  dsRNASeeker_strict_gate_genes = intersect(strict_candidate_genes, names(gene_ranks)),
  dsRNASeeker_priority_top20_genes = intersect(top20_candidate_genes, names(gene_ranks)),
  dsRNASeeker_structure_top20_genes = intersect(structure_candidate_genes, names(gene_ranks))
)
candidate_sets <- candidate_sets[lengths(candidate_sets) >= 5]
write_tsv(tibble(
  candidate_set = names(candidate_sets),
  n_genes_in_ranked_universe = lengths(candidate_sets),
  genes = vapply(candidate_sets, paste, collapse = ";", FUN.VALUE = character(1))
), file.path(outdir, "candidate_gene_sets.tsv"))

message("[INFO] Running Hallmark fgsea on the full gene-ranked universe")

# MSigDB uses different Hallmark collection codes for its two native databases:
#   human database (HS): H
#   mouse database (MM): MH
# Requesting H from db_species = "MM" raises "Unknown collection" in current
# msigdbr releases. Query the installed package first so failures remain explicit
# and reproducible if MSigDB changes its collection codes in the future.
hallmark_collection <- if (db_species == "MM") "MH" else "H"
available_collections <- msigdbr_collections(db_species = db_species)
collection_column <- intersect(c("gs_collection", "gs_cat"), names(available_collections))
if (length(collection_column) == 0) {
  stop("Could not identify the collection-code column returned by msigdbr_collections(). Columns: ",
       paste(names(available_collections), collapse = ", "))
}
collection_column <- collection_column[[1]]
available_codes <- sort(unique(as.character(available_collections[[collection_column]])))

message("[INFO] msigdbr version: ", as.character(utils::packageVersion("msigdbr")))
message("[INFO] MSigDB database: ", db_species,
        "; Hallmark collection: ", hallmark_collection)

if (!hallmark_collection %in% available_codes) {
  stop("Hallmark collection '", hallmark_collection,
       "' is unavailable for db_species='", db_species,
       "'. Available collection codes: ", paste(available_codes, collapse = ", "),
       ". Run msigdbr::msigdbr_collections(db_species = '", db_species,
       "') in this R environment to inspect the installed database.")
}

msig <- msigdbr(
  db_species = db_species,
  species = msig_species,
  collection = hallmark_collection
)
if (nrow(msig) == 0) {
  stop("msigdbr returned zero genes for database ", db_species,
       " and collection ", hallmark_collection)
}

identifier_overlap_n <- length(intersect(names(gene_ranks), unique(msig$gene_symbol)))
identifier_overlap_fraction <- identifier_overlap_n / length(gene_ranks)
if (identifier_overlap_n < 1000 || identifier_overlap_fraction < 0.10) {
  stop("Insufficient gene-symbol overlap with MSigDB: n=", identifier_overlap_n,
       ", fraction=", signif(identifier_overlap_fraction, 3),
       ". Check the GTF gene_id-to-gene_name mapping before interpreting GSEA.")
}
pathways <- split(msig$gene_symbol, msig$gs_name)
set.seed(123)
fg <- fgseaMultilevel(
  pathways = pathways,
  stats = gene_ranks,
  minSize = min_size,
  maxSize = max_size,
  eps = 0,
  nproc = 1
) %>% as_tibble() %>% arrange(padj)

fg_export <- fg %>% mutate(
  dataset = dataset_label,
  case = case_label,
  control = control_label,
  leadingEdge = vapply(leadingEdge, paste, collapse = ";", FUN.VALUE = character(1)),
  candidate_leading_edge = vapply(strsplit(leadingEdge, ";", fixed = TRUE), function(x) paste(intersect(x, all_candidate_genes), collapse = ";"), character(1)),
  candidate_leading_edge_n = vapply(strsplit(candidate_leading_edge, ";", fixed = TRUE), function(x) as.integer(sum(x != "")), integer(1))
)
write_csv(fg_export, file.path(outdir, "hallmark_fgsea_all.csv"))
write_csv(fg_export %>% filter(!is.na(padj), padj <= 0.05), file.path(outdir, "hallmark_fgsea_significant.csv"))

if (length(candidate_sets)) {
  custom <- fgseaMultilevel(
    pathways = candidate_sets,
    stats = gene_ranks,
    minSize = 5,
    maxSize = max(5000, length(gene_ranks)),
    eps = 0,
    nproc = 1
  ) %>% as_tibble() %>% arrange(padj) %>%
    mutate(dataset = dataset_label,
           leadingEdge = vapply(leadingEdge, paste, collapse = ";", FUN.VALUE = character(1)))
  write_csv(custom, file.path(outdir, "candidate_gene_set_fgsea.csv"))
}

plot_df <- fg_export %>%
  filter(!is.na(padj), is.finite(padj)) %>%
  arrange(padj, desc(abs(NES))) %>%
  slice_head(n = top_n) %>%
  mutate(pathway_label = str_replace_all(pathway, "^HALLMARK_", ""),
         pathway_label = str_replace_all(pathway_label, "_", " "),
         pathway_label = factor(pathway_label, levels = rev(pathway_label)),
         sig = -log10(pmax(padj, 1e-300)))

p <- ggplot(plot_df, aes(x = NES, y = pathway_label)) +
  geom_vline(xintercept = 0, linewidth = 0.4, linetype = 2) +
  geom_point(aes(size = sig, fill = NES), shape = 21, color = "black", stroke = 0.35) +
  geom_text(aes(label = ifelse(candidate_leading_edge_n > 0, candidate_leading_edge_n, "")),
            size = 3, color = "black") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, name = "NES") +
  scale_size_continuous(name = expression(-log[10](FDR)), range = c(3, 10)) +
  labs(
    title = paste0(dataset_label, ": Hallmark GSEA"),
    subtitle = "Numbers inside bubbles = dsRNASeeker candidate-associated leading-edge genes",
    x = paste0("Normalized enrichment score (", case_label, " vs ", control_label, ")"),
    y = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 10))

ggsave(file.path(outdir, "hallmark_GSEA_candidate_linked_bubble.svg"), p, width = 12, height = 8)
ggsave(file.path(outdir, "hallmark_GSEA_candidate_linked_bubble.png"), p, width = 12, height = 8, dpi = 600)
ggsave(file.path(outdir, "hallmark_GSEA_candidate_linked_bubble.pdf"), p, width = 12, height = 8)
ggsave(file.path(outdir, "hallmark_GSEA_candidate_linked_bubble.eps"), p, width = 7.0, height = 4.8,device = cairo_ps, units = "in")

metadata <- tibble(
  dataset = dataset_label,
  species = species_key,
  case = case_label,
  control = control_label,
  input_gene_count_rows = nrow(counts_df),
  retained_ranked_genes = length(gene_ranks),
  candidate_genes_in_summary = length(all_candidate_genes),
  candidate_genes_in_ranked_universe = length(intersect(all_candidate_genes, names(gene_ranks))),
  msigdb_identifier_overlap_n = identifier_overlap_n,
  msigdb_identifier_overlap_fraction = identifier_overlap_fraction,
  hallmark_version = unique(msig$db_version)[1]
)
write_csv(metadata, file.path(outdir, "run_metadata.csv"))
message("[OK] GSEA outputs: ", outdir)
