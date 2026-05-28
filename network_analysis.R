#!/usr/bin/env Rscript

################################################################################
# Co-occurrence network analysis for amplicon ASV tables
#
# Input:
#   1. ASV table: rows = ASVs/features, columns = samples
#   2. Metadata: first column = SampleID, with a grouping column
#
# Main steps:
#   - Split ASV table by group
#   - Convert counts to relative abundance within each sample
#   - Calculate Spearman correlations among ASVs
#   - Adjust P values using Benjamini-Hochberg FDR correction
#   - Keep strong and significant correlations
#   - Build undirected igraph networks
#   - Remove isolated nodes
#   - Calculate network modules and basic topological properties
#   - Export edge tables, node tables, network property table, and PDF network plots
################################################################################

suppressPackageStartupMessages({
  library(igraph)
  library(Hmisc)
  library(dplyr)
  library(readr)
})

# ----------------------------- user parameters --------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop(
    "Usage: Rscript network_analysis.R <asv_table.txt> <metadata.txt> <group_column>\n",
    "Example: Rscript network_analysis.R asv_table.txt metadata.txt Group"
  )
}

asv_file <- args[1]
meta_file <- args[2]
group_col <- args[3]

outdir <- "network_results"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "edges"), showWarnings = FALSE)
dir.create(file.path(outdir, "nodes"), showWarnings = FALSE)
dir.create(file.path(outdir, "plots"), showWarnings = FALSE)

rho_cutoff <- 0.8
p_cutoff <- 0.05
min_prevalence <- 2
min_total_abundance <- 1
top_module_number <- 20
layout_seed <- 7

# Colors used for modules. Edges are colored by correlation sign in the final plot.
module_cols <- c(
  "#FFB300", "#803E75", "#FF6800", "#A6BDD7", "#117F9D",
  "#C10020", "#CEA262", "#817066", "#007D34", "#F6768E",
  "#00538A", "#FF7A5C", "#53377A", "#FF8E00", "#B32851",
  "#F4C800", "#7F180D", "#93AA00", "#593315", "#F13A13",
  "#232C16", "#DCD300", "#8CA0A0", "#A4E300", "#006E00",
  "#A1CAF1", "#D783FF", "#92C5DE", "#FFC857", "#119DA4"
)
grey_col <- "#C1C1C1"

# ----------------------------- helper functions -------------------------------

read_asv_table <- function(file) {
  x <- read.delim(file, header = TRUE, row.names = 1, check.names = FALSE, sep = "\t")
  x <- as.data.frame(x)
  x[] <- lapply(x, as.numeric)
  x[is.na(x)] <- 0
  return(x)
}

read_metadata <- function(file) {
  meta <- read.delim(file, header = TRUE, check.names = FALSE, sep = "\t")
  if (!"SampleID" %in% colnames(meta)) {
    colnames(meta)[1] <- "SampleID"
  }
  return(meta)
}

make_network <- function(otu, rho_cutoff = 0.8, p_cutoff = 0.05) {
  # Remove rare features before correlation
  otu <- otu[rowSums(otu) > min_total_abundance & rowSums(otu > 0) >= min_prevalence, , drop = FALSE]

  if (nrow(otu) < 3 || ncol(otu) < 3) {
    return(NULL)
  }

  # Convert count table to relative abundance by sample
  otu_rel <- sweep(otu, 2, colSums(otu), FUN = "/")
  otu_rel[is.na(otu_rel)] <- 0

  # Hmisc::rcorr requires variables in columns, samples in rows
  cor_res <- Hmisc::rcorr(t(as.matrix(otu_rel)), type = "spearman")
  r_mat <- cor_res$r
  p_mat <- cor_res$P

  diag(r_mat) <- 0
  diag(p_mat) <- 1

  # BH correction across all pairwise tests
  p_adj <- matrix(
    p.adjust(as.vector(p_mat), method = "BH"),
    nrow = nrow(p_mat),
    ncol = ncol(p_mat),
    dimnames = dimnames(p_mat)
  )

  r_mat[p_adj > p_cutoff | abs(r_mat) < rho_cutoff] <- 0
  r_mat[is.na(r_mat)] <- 0

  g <- graph_from_adjacency_matrix(r_mat, mode = "undirected", weighted = TRUE, diag = FALSE)
  g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
  g <- delete_vertices(g, which(degree(g) == 0))

  if (vcount(g) < 3 || ecount(g) < 2) {
    return(NULL)
  }

  E(g)$correlation <- E(g)$weight
  E(g)$sign <- ifelse(E(g)$correlation > 0, "positive", "negative")
  E(g)$weight <- abs(E(g)$correlation)

  set.seed(layout_seed)
  V(g)$module <- membership(cluster_fast_greedy(g, weights = E(g)$weight))

  module_size <- sort(table(V(g)$module), decreasing = TRUE)
  top_modules <- names(module_size)[seq_len(min(top_module_number, length(module_size)))]

  module_color_map <- module_cols[seq_along(top_modules)]
  names(module_color_map) <- top_modules

  V(g)$color <- grey_col
  V(g)$color[as.character(V(g)$module) %in% top_modules] <-
    module_color_map[as.character(V(g)$module[as.character(V(g)$module) %in% top_modules])]
  V(g)$frame.color <- V(g)$color
  V(g)$label <- NA

  E(g)$color <- ifelse(E(g)$correlation > 0, "#D73027", "#4575B4")

  return(g)
}

calc_properties <- function(g, group_name) {
  if (is.null(g)) return(NULL)

  positive_edges <- sum(E(g)$correlation > 0)
  negative_edges <- sum(E(g)$correlation < 0)

  comp <- components(g)
  largest_nodes <- names(comp$membership[comp$membership == which.max(comp$csize)])
  g_largest <- induced_subgraph(g, vids = largest_nodes)

  data.frame(
    Group = group_name,
    Nodes = vcount(g),
    Edges = ecount(g),
    Positive_edges = positive_edges,
    Negative_edges = negative_edges,
    Average_degree = mean(degree(g)),
    Density = edge_density(g, loops = FALSE),
    Average_clustering_coefficient = transitivity(g, type = "average", isolates = "zero"),
    Modularity = modularity(cluster_fast_greedy(g, weights = E(g)$weight)),
    Average_path_length_largest_component = average.path.length(g_largest, directed = FALSE),
    Diameter_largest_component = diameter(g_largest, directed = FALSE),
    stringsAsFactors = FALSE
  )
}

export_network <- function(g, group_name) {
  if (is.null(g)) return(NULL)

  edge_df <- as_data_frame(g, what = "edges") %>%
    rename(Source = from, Target = to) %>%
    mutate(Group = group_name)

  node_df <- data.frame(
    ASV = V(g)$name,
    Group = group_name,
    Degree = degree(g),
    Module = V(g)$module,
    Color = V(g)$color,
    stringsAsFactors = FALSE
  )

  write.table(edge_df, file.path(outdir, "edges", paste0(group_name, "_edges.txt")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(node_df, file.path(outdir, "nodes", paste0(group_name, "_nodes.txt")),
              sep = "\t", row.names = FALSE, quote = FALSE)

  pdf(file.path(outdir, "plots", paste0(group_name, "_network.pdf")),
      width = 6, height = 6, encoding = "MacRoman")
  set.seed(layout_seed)
  lay <- layout_with_fr(g, niter = 1000, grid = "nogrid")
  par(mar = c(0, 0, 2, 0), font.main = 4)
  plot(
    g,
    layout = lay,
    vertex.size = 2,
    vertex.label = NA,
    vertex.color = V(g)$color,
    vertex.frame.color = V(g)$frame.color,
    edge.color = E(g)$color,
    edge.width = 0.5
  )
  title(main = paste0(group_name, ": Nodes=", vcount(g), ", Edges=", ecount(g)))
  dev.off()
}

# ----------------------------- main analysis ----------------------------------

asv <- read_asv_table(asv_file)
meta <- read_metadata(meta_file)

if (!group_col %in% colnames(meta)) {
  stop("The group column was not found in metadata: ", group_col)
}

common_samples <- intersect(colnames(asv), meta$SampleID)
if (length(common_samples) < 3) {
  stop("Too few matched samples between ASV table and metadata.")
}

asv <- asv[, common_samples, drop = FALSE]
meta <- meta[match(common_samples, meta$SampleID), , drop = FALSE]

groups <- unique(meta[[group_col]])
groups <- groups[!is.na(groups)]

all_properties <- list()

for (grp in groups) {
  message("Processing group: ", grp)

  sample_use <- meta$SampleID[meta[[group_col]] == grp]
  otu_sub <- asv[, sample_use, drop = FALSE]

  g <- make_network(otu_sub, rho_cutoff = rho_cutoff, p_cutoff = p_cutoff)

  if (is.null(g)) {
    warning("No valid network was generated for group: ", grp)
    next
  }

  safe_grp <- gsub("[^A-Za-z0-9_\\-]", "_", grp)

  export_network(g, safe_grp)
  all_properties[[safe_grp]] <- calc_properties(g, safe_grp)
}

properties_df <- bind_rows(all_properties)
write.table(properties_df, file.path(outdir, "network_properties.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

message("Network analysis finished. Results are saved in: ", outdir)
