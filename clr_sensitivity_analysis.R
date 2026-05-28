#!/usr/bin/env Rscript

################################################################################
# CLR-based sensitivity analysis for amplicon ASV data
#
# This script performs a compositional-data sensitivity analysis using:
#   - ASV filtering
#   - pseudocount addition
#   - centered log-ratio (CLR) transformation
#   - Euclidean distance on CLR-transformed data, equivalent to Aitchison distance
#   - PCoA visualization
#   - PERMANOVA for overall compartment effects
#   - PERMANOVA within each compartment to test health effects
#
# Input:
#   1. ASV table: rows = ASVs/features, columns = samples
#   2. Metadata: must contain sample ID, Compartment, and Health columns
#
# Example:
#   Rscript clr_sensitivity_analysis.R ASV_table_bacteria.csv metadata.csv bacteria
#   Rscript clr_sensitivity_analysis.R ASV_table_fungi.csv metadata.csv fungi
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(ape)
  library(compositions)
  library(ggplot2)
})

# ----------------------------- arguments --------------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop(
    "Usage: Rscript clr_sensitivity_analysis.R <ASV_table.csv> <metadata.csv> <prefix>\n",
    "Example: Rscript clr_sensitivity_analysis.R ASV_table_bacteria.csv metadata.csv bacteria"
  )
}

asv_file <- args[1]
meta_file <- args[2]
prefix <- args[3]

# ----------------------------- parameters -------------------------------------

min_total_abundance <- 10
min_prevalence <- 2
pseudocount <- 1
permutations <- 999

compartment_levels <- c("Leaf endosphere", "Rhizosphere", "Root endosphere")
health_levels <- c("Diseased", "Healthy")

outdir <- paste0("CLR_sensitivity_", prefix)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------- read data --------------------------------------

asv <- read.csv(
  asv_file,
  row.names = 1,
  check.names = FALSE
)

meta <- read.csv(
  meta_file,
  check.names = FALSE
)

# The first column of metadata is assumed to be the sample ID if "ID" is absent.
if (!"ID" %in% colnames(meta)) {
  colnames(meta)[1] <- "ID"
}

meta <- meta %>%
  column_to_rownames("ID")

# Match samples between ASV table and metadata.
common_samples <- intersect(colnames(asv), rownames(meta))

if (length(common_samples) < 3) {
  stop("Too few matched samples between ASV table and metadata.")
}

asv <- asv[, common_samples, drop = FALSE]
meta <- meta[common_samples, , drop = FALSE]

if (!all(colnames(asv) == rownames(meta))) {
  stop("Sample order between ASV table and metadata is not consistent.")
}

if (!all(c("Compartment", "Health") %in% colnames(meta))) {
  stop("Metadata must contain 'Compartment' and 'Health' columns.")
}

# ----------------------------- ASV filtering ----------------------------------

asv_filt <- asv[
  rowSums(asv) > min_total_abundance &
    rowSums(asv > 0) >= min_prevalence,
  ,
  drop = FALSE
]

if (nrow(asv_filt) < 3) {
  stop("Too few ASVs remained after filtering.")
}

# Samples x ASVs
asv_t <- t(asv_filt)

# ----------------------------- CLR transformation -----------------------------

asv_pc <- asv_t + pseudocount

asv_clr <- compositions::clr(asv_pc)
asv_clr <- as.data.frame(asv_clr)

# Euclidean distance calculated from CLR-transformed data is Aitchison distance.
dist_clr <- dist(asv_clr, method = "euclidean")

# ----------------------------- PCoA -------------------------------------------

pcoa_res <- ape::pcoa(dist_clr)

pcoa_df <- data.frame(
  ID = rownames(asv_clr),
  PCoA1 = pcoa_res$vectors[, 1],
  PCoA2 = pcoa_res$vectors[, 2]
) %>%
  left_join(
    meta %>% rownames_to_column("ID"),
    by = "ID"
  )

var_exp <- round(100 * pcoa_res$values$Relative_eig[1:2], 2)

pcoa_df$Compartment <- factor(
  pcoa_df$Compartment,
  levels = compartment_levels
)

pcoa_df$Health <- factor(
  pcoa_df$Health,
  levels = health_levels
)

p <- ggplot(
  pcoa_df,
  aes(
    x = PCoA1,
    y = PCoA2,
    color = Compartment,
    shape = Health
  )
) +
  geom_point(size = 3, alpha = 0.85) +
  stat_ellipse(
    aes(group = Compartment, color = Compartment),
    type = "t",
    linewidth = 0.7,
    linetype = 2
  ) +
  labs(
    x = paste0("PCoA 1 (", var_exp[1], "%)"),
    y = paste0("PCoA 2 (", var_exp[2], "%)"),
    color = "Compartment",
    shape = "Health"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "right"
  )

ggsave(
  filename = file.path(outdir, paste0(prefix, "_CLR_PCoA.pdf")),
  plot = p,
  width = 5.67,
  height = 4.32
)

ggsave(
  filename = file.path(outdir, paste0(prefix, "_CLR_PCoA.png")),
  plot = p,
  width = 6,
  height = 5,
  dpi = 300
)

write.csv(
  pcoa_df,
  file.path(outdir, paste0(prefix, "_CLR_PCoA_coordinates.csv")),
  row.names = FALSE
)

# ----------------------------- PERMANOVA --------------------------------------

adonis_compartment <- adonis2(
  dist_clr ~ Compartment,
  data = meta,
  permutations = permutations
)

write.csv(
  as.data.frame(adonis_compartment),
  file.path(outdir, paste0(prefix, "_CLR_PERMANOVA_compartment.csv"))
)

permanova_by_compartment <- lapply(unique(meta$Compartment), function(comp) {

  samples_use <- rownames(meta)[meta$Compartment == comp]

  dist_sub <- as.dist(as.matrix(dist_clr)[samples_use, samples_use])
  meta_sub <- meta[samples_use, , drop = FALSE]

  res <- adonis2(
    dist_sub ~ Health,
    data = meta_sub,
    permutations = permutations
  )

  out <- as.data.frame(res)
  out$Compartment <- comp
  out$Term <- rownames(out)

  return(out)
}) %>%
  bind_rows()

write.csv(
  permanova_by_compartment,
  file.path(outdir, paste0(prefix, "_CLR_PERMANOVA_health_by_compartment.csv")),
  row.names = FALSE
)

# ----------------------------- session info -----------------------------------

sink(file.path(outdir, paste0(prefix, "_sessionInfo.txt")))
sessionInfo()
sink()

message("CLR sensitivity analysis finished. Results are saved in: ", outdir)
