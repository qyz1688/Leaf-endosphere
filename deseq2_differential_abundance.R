#!/usr/bin/env Rscript

################################################################################
# DESeq2 differential abundance analysis for amplicon ASV/taxa count tables
#
# This script performs differential abundance analysis using DESeq2.
#
# Main steps:
#   - Read count table and metadata
#   - Match samples between count table and metadata
#   - Construct DESeq2 dataset
#   - Filter low-count features
#   - Run DESeq2 normalization and differential abundance testing
#   - Export full and significant results
#   - Generate bar plot and volcano plot
#
# Input:
#   1. Count table: rows = ASVs/taxa/features, columns = samples
#   2. Metadata: rows = samples; must contain sample ID and grouping column
#
# Example:
#   Rscript deseq2_differential_abundance.R count_table.csv metadata.csv Group A R bacteria
#
# Arguments:
#   1. count_table.csv   Feature count table
#   2. metadata.csv      Metadata table
#   3. Group             Grouping column in metadata
#   4. A                 Numerator group in contrast
#   5. R                 Denominator/reference group in contrast
#   6. bacteria          Output prefix
################################################################################

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(ggplot2)
  library(ggrepel)
})

# ----------------------------- arguments --------------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 6) {
  stop(
    "Usage: Rscript deseq2_differential_abundance.R <count_table.csv> <metadata.csv> <group_column> <contrast_numerator> <contrast_denominator> <prefix>\n",
    "Example: Rscript deseq2_differential_abundance.R test_otudian.csv test_design_dian.csv Group A R bacteria"
  )
}

count_file <- args[1]
meta_file <- args[2]
group_col <- args[3]
contrast_num <- args[4]
contrast_den <- args[5]
prefix <- args[6]

# ----------------------------- parameters -------------------------------------

min_count <- 1
min_sample_number <- 3
padj_cutoff <- 0.05
log2fc_cutoff <- log2(1.5)
top_label_number <- 5

outdir <- paste0("DESeq2_results_", prefix)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------- read data --------------------------------------

count_table <- read.csv(
  count_file,
  row.names = 1,
  check.names = FALSE
)

metadata <- read.csv(
  meta_file,
  check.names = FALSE
)

# If Sample column is absent, treat the first metadata column as sample ID.
if (!"Sample" %in% colnames(metadata)) {
  colnames(metadata)[1] <- "Sample"
}

if (!group_col %in% colnames(metadata)) {
  stop("Grouping column was not found in metadata: ", group_col)
}

rownames(metadata) <- metadata$Sample

# Match samples between count table and metadata.
common_samples <- intersect(colnames(count_table), rownames(metadata))

if (length(common_samples) < 3) {
  stop("Too few matched samples between count table and metadata.")
}

count_table <- count_table[, common_samples, drop = FALSE]
metadata <- metadata[common_samples, , drop = FALSE]

if (!all(colnames(count_table) == rownames(metadata))) {
  stop("Sample order between count table and metadata is not consistent.")
}

# Ensure integer count matrix.
count_table <- round(as.matrix(count_table))
storage.mode(count_table) <- "integer"

# Convert grouping column to factor and set denominator as reference.
metadata[[group_col]] <- factor(metadata[[group_col]])
metadata[[group_col]] <- relevel(metadata[[group_col]], ref = contrast_den)

# ----------------------------- DESeq2 analysis --------------------------------

dds <- DESeqDataSetFromMatrix(
  countData = count_table,
  colData = metadata,
  design = as.formula(paste0("~ ", group_col))
)

# Filter low-count features.
keep <- rowSums(counts(dds) >= min_count) >= min_sample_number
dds <- dds[keep, ]

dds <- DESeq(
  dds,
  fitType = "mean",
  minReplicatesForReplace = 7,
  parallel = FALSE
)

res <- results(
  dds,
  contrast = c(group_col, contrast_num, contrast_den)
)

res_df <- as.data.frame(res) %>%
  rownames_to_column("Feature") %>%
  arrange(padj)

write.csv(
  res_df,
  file.path(outdir, paste0(prefix, "_DESeq2_all_results.csv")),
  row.names = FALSE,
  quote = FALSE
)

sig_df <- res_df %>%
  filter(!is.na(padj)) %>%
  mutate(
    Change = case_when(
      padj < padj_cutoff & log2FoldChange > log2fc_cutoff ~ "Up",
      padj < padj_cutoff & log2FoldChange < -log2fc_cutoff ~ "Down",
      TRUE ~ "Normal"
    )
  )

write.csv(
  sig_df,
  file.path(outdir, paste0(prefix, "_DESeq2_results_with_change.csv")),
  row.names = FALSE,
  quote = FALSE
)

sig_only <- sig_df %>%
  filter(Change %in% c("Up", "Down"))

write.csv(
  sig_only,
  file.path(outdir, paste0(prefix, "_DESeq2_significant_results.csv")),
  row.names = FALSE,
  quote = FALSE
)

# ----------------------------- bar plot ---------------------------------------

bar_df <- sig_df %>%
  filter(!is.na(padj), padj < padj_cutoff) %>%
  arrange(log2FoldChange)

if (nrow(bar_df) > 0) {
  bar_df$Direction <- ifelse(bar_df$log2FoldChange >= 0, "Enriched", "Depleted")
  bar_df$Feature <- factor(bar_df$Feature, levels = bar_df$Feature)

  p_bar <- ggplot(bar_df, aes(x = Feature, y = log2FoldChange, fill = Direction)) +
    geom_col(width = 0.85, color = "grey30") +
    coord_flip() +
    scale_fill_manual(values = c("Enriched" = "#ED0000", "Depleted" = "#3F88C5")) +
    labs(
      x = NULL,
      y = "log2 fold change",
      fill = NULL
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_text(size = 7)
    )

  ggsave(
    file.path(outdir, paste0(prefix, "_DESeq2_log2FC_barplot.pdf")),
    p_bar,
    width = 7,
    height = max(4, 0.18 * nrow(bar_df))
  )
}

# ----------------------------- volcano plot -----------------------------------

plot_df <- sig_df %>%
  drop_na(padj, log2FoldChange) %>%
  mutate(
    neg_log10_padj = -log10(padj),
    Change = case_when(
      padj < padj_cutoff & log2FoldChange > log2fc_cutoff ~ "Up",
      padj < padj_cutoff & log2FoldChange < -log2fc_cutoff ~ "Down",
      TRUE ~ "Normal"
    )
  )

top_up <- plot_df %>%
  filter(Change == "Up") %>%
  arrange(padj) %>%
  slice_head(n = top_label_number)

top_down <- plot_df %>%
  filter(Change == "Down") %>%
  arrange(padj) %>%
  slice_head(n = top_label_number)

label_df <- bind_rows(top_up, top_down)

p_volcano <- ggplot(plot_df, aes(x = log2FoldChange, y = neg_log10_padj)) +
  geom_point(aes(color = Change), size = 2.5, alpha = 0.7) +
  scale_color_manual(
    values = c("Up" = "#ED0000", "Down" = "#3F88C5", "Normal" = "grey70")
  ) +
  geom_text_repel(
    data = label_df,
    aes(label = Feature),
    box.padding = 0.5,
    max.overlaps = Inf,
    segment.color = "grey50",
    size = 3
  ) +
  geom_vline(
    xintercept = c(-log2fc_cutoff, log2fc_cutoff),
    linetype = "dashed",
    color = "grey30"
  ) +
  geom_hline(
    yintercept = -log10(padj_cutoff),
    linetype = "dashed",
    color = "grey30"
  ) +
  labs(
    x = "log2 fold change",
    y = "-log10 adjusted P value",
    color = NULL
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "right"
  )

ggsave(
  file.path(outdir, paste0(prefix, "_DESeq2_volcano.pdf")),
  p_volcano,
  width = 6,
  height = 5
)

# ----------------------------- normalized counts ------------------------------

norm_counts <- counts(dds, normalized = TRUE) %>%
  as.data.frame() %>%
  rownames_to_column("Feature")

write.csv(
  norm_counts,
  file.path(outdir, paste0(prefix, "_DESeq2_normalized_counts.csv")),
  row.names = FALSE,
  quote = FALSE
)

# ----------------------------- session info -----------------------------------

sink(file.path(outdir, paste0(prefix, "_sessionInfo.txt")))
sessionInfo()
sink()

message("DESeq2 differential abundance analysis finished. Results are saved in: ", outdir)
