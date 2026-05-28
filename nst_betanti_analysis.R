#!/usr/bin/env Rscript

################################################################################
# NST / betaMNTD / betaNTI analysis pipeline
#
# This script calculates phylogenetic community assembly metrics using NST::pNST.
# It outputs:
#   - pairwise betaMNTD and betaNTI values
#   - within-group stochastic contribution: proportion of |betaNTI| < 2
#   - within-group deterministic contribution: proportion of |betaNTI| > 2
#
# Input:
#   1. OTU/ASV table: rows = ASVs/OTUs, columns = samples
#   2. Group file: rows = samples; must contain sample ID and grouping column
#   3. Phylogenetic tree in Newick format
#
# Example:
#   Rscript nst_betanti_analysis.R eyeotu.csv groupeye.csv Eye_bac_all.nwk treatment bacteria
################################################################################

suppressPackageStartupMessages({
  library(NST)
  library(picante)
  library(ape)
  library(dplyr)
  library(readr)
})

# ----------------------------- arguments --------------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  stop(
    "Usage: Rscript nst_betanti_analysis.R <otu_table.csv> <group_file.csv> <tree.nwk> <group_column> <prefix>\n",
    "Example: Rscript nst_betanti_analysis.R eyeotu.csv groupeye.csv Eye_bac_all.nwk treatment bacteria"
  )
}

otu_file <- args[1]
group_file <- args[2]
tree_file <- args[3]
group_col <- args[4]
prefix <- args[5]

# ----------------------------- parameters -------------------------------------

rand_number <- 1000
nworker <- 1
set_seed <- 123
abundance_weighted <- TRUE

outdir <- paste0("NST_betaNTI_", prefix)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------- read data --------------------------------------

otu_raw <- read.csv(
  otu_file,
  row.names = 1,
  check.names = FALSE
)

# Input OTU table is assumed to be features x samples.
# pNST requires communities as samples x species.
otu <- data.frame(t(otu_raw), check.names = FALSE)

group <- read.csv(
  group_file,
  check.names = FALSE
)

# If the first column is not named SampleID, rename it.
if (!"SampleID" %in% colnames(group)) {
  colnames(group)[1] <- "SampleID"
}

if (!group_col %in% colnames(group)) {
  stop("The grouping column was not found in the group file: ", group_col)
}

rownames(group) <- group$SampleID

# Match samples between OTU table and group file.
common_samples <- intersect(rownames(otu), rownames(group))

if (length(common_samples) < 3) {
  stop("Too few matched samples between OTU table and group file.")
}

otu <- otu[common_samples, , drop = FALSE]
group <- group[common_samples, , drop = FALSE]

# Remove ASVs/OTUs absent from all matched samples.
otu <- otu[, colSums(otu) > 0, drop = FALSE]

# Read and prune tree.
tree <- read.tree(tree_file)
tree <- prune.sample(otu, tree)

# Keep only taxa present in both OTU table and tree.
common_taxa <- intersect(colnames(otu), tree$tip.label)

if (length(common_taxa) < 3) {
  stop("Too few matched taxa between OTU table and phylogenetic tree.")
}

otu <- otu[, common_taxa, drop = FALSE]
tree <- keep.tip(tree, common_taxa)

# Ensure group file is ordered exactly as OTU table.
group <- group[rownames(otu), , drop = FALSE]

# ----------------------------- pNST analysis ----------------------------------

set.seed(set_seed)

pnst <- pNST(
  comm = otu,
  tree = tree,
  group = group[, group_col, drop = FALSE],
  phylo.shuffle = TRUE,
  taxo.null.model = NULL,
  pd.wd = tempdir(),
  abundance.weighted = abundance_weighted,
  rand = rand_number,
  nworker = nworker,
  SES = TRUE,
  RC = FALSE
)

saveRDS(
  pnst,
  file = file.path(outdir, paste0(prefix, "_pNST_result.rds"))
)

# Pairwise betaMNTD and betaNTI results.
betaMNTD <- pnst$index.pair

write.csv(
  betaMNTD,
  file.path(outdir, paste0(prefix, "_pairwise_betaMNTD_betaNTI.csv")),
  row.names = FALSE,
  quote = FALSE
)

# ----------------------------- within-group summary ---------------------------

group_levels <- unique(group[[group_col]])
summary_list <- list()

for (grp in group_levels) {

  samples_use <- rownames(group)[group[[group_col]] == grp]

  beta_sub <- betaMNTD %>%
    filter(name1 %in% samples_use & name2 %in% samples_use)

  if (nrow(beta_sub) == 0) {
    warning("No within-group pairwise comparison found for group: ", grp)
    next
  }

  # Choose abundance-weighted betaNTI if available.
  if ("bNTI.wt" %in% colnames(beta_sub)) {
    beta_nti <- beta_sub$bNTI.wt
    beta_nti_col <- "bNTI.wt"
  } else if ("bNTI" %in% colnames(beta_sub)) {
    beta_nti <- beta_sub$bNTI
    beta_nti_col <- "bNTI"
  } else {
    stop("No betaNTI column found in pNST index.pair output.")
  }

  total_pairs <- length(beta_nti)
  stochastic_pairs <- sum(abs(beta_nti) < 2, na.rm = TRUE)
  deterministic_pairs <- sum(abs(beta_nti) > 2, na.rm = TRUE)

  summary_list[[as.character(grp)]] <- data.frame(
    Group = grp,
    BetaNTI_column = beta_nti_col,
    Pair_number = total_pairs,
    Stochastic_pair_number = stochastic_pairs,
    Deterministic_pair_number = deterministic_pairs,
    Stochastic_contribution = stochastic_pairs / total_pairs,
    Deterministic_contribution = deterministic_pairs / total_pairs,
    stringsAsFactors = FALSE
  )
}

summary_df <- bind_rows(summary_list)

write.csv(
  summary_df,
  file.path(outdir, paste0(prefix, "_betaNTI_group_summary.csv")),
  row.names = FALSE,
  quote = FALSE
)

# ----------------------------- session info -----------------------------------

sink(file.path(outdir, paste0(prefix, "_sessionInfo.txt")))
sessionInfo()
sink()

message("NST / betaNTI analysis finished. Results are saved in: ", outdir)
