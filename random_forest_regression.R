#!/usr/bin/env Rscript

################################################################################
# Random Forest regression and core-feature association analysis
#
# This script performs Random Forest regression using microbial features or
# core ASVs as explanatory variables and a plant phenotype as the response.
#
# Main steps:
#   - Read feature/phenotype table
#   - Optionally subset samples by a metadata/grouping column
#   - Train Random Forest regression model
#   - Predict phenotype values using the fitted model
#   - Plot observed vs predicted values
#   - Export variable importance table
#   - Plot top important variables
#   - Perform repeated cross-validation with rfcv
#   - Test correlations between selected variables and phenotype
#
# Input:
#   A tab-delimited table with samples as rows and variables as columns.
#   The response variable and optional grouping column must be included.
#
# Example:
#   Rscript random_forest_regression.R data.txt SW crop wheat 10
#
# Arguments:
#   1. data.txt          Input table
#   2. SW                Response variable
#   3. crop              Optional grouping column; use NA if not needed
#   4. wheat             Optional group level; use NA if not needed
#   5. 10                Number of top variables to export/plot
################################################################################

suppressPackageStartupMessages({
  library(randomForest)
  library(tidyverse)
  library(ggplot2)
  library(ggpubr)
  library(ggpmisc)
  library(ggExtra)
  library(reshape2)
})

# ----------------------------- arguments --------------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  stop(
    "Usage: Rscript random_forest_regression.R <data.txt> <response_variable> <group_column_or_NA> <group_level_or_NA> <top_n>\n",
    "Example: Rscript random_forest_regression.R data.txt SW crop wheat 10\n",
    "Example without subsetting: Rscript random_forest_regression.R data.txt Leaf_area NA NA 10"
  )
}

data_file <- args[1]
response_var <- args[2]
group_col <- args[3]
group_level <- args[4]
top_n <- as.numeric(args[5])

# ----------------------------- parameters -------------------------------------

set.seed(123)
cv_seed <- 111
cv_repeats <- 5
cv_fold <- 5
cv_step <- 0.8

outdir <- paste0("RF_results_", response_var)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------- read and prepare data ---------------------------

dat <- read.table(
  data_file,
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE
)

if (!response_var %in% colnames(dat)) {
  stop("Response variable was not found in the input table: ", response_var)
}

# Optional sample subsetting by a grouping column.
if (group_col != "NA" && group_level != "NA") {
  if (!group_col %in% colnames(dat)) {
    stop("Grouping column was not found in the input table: ", group_col)
  }
  dat <- dat[dat[[group_col]] == group_level, , drop = FALSE]
  dat <- dat[, setdiff(colnames(dat), group_col), drop = FALSE]
}

# Keep only complete cases.
dat <- dat[complete.cases(dat), , drop = FALSE]

# Convert all predictors to numeric when possible.
predictor_cols <- setdiff(colnames(dat), response_var)
dat[predictor_cols] <- lapply(dat[predictor_cols], function(x) as.numeric(as.character(x)))

# Remove predictors with all NA or zero variance.
valid_predictors <- predictor_cols[
  sapply(dat[predictor_cols], function(x) sum(!is.na(x)) > 2 && sd(x, na.rm = TRUE) > 0)
]

dat <- dat[, c(valid_predictors, response_var), drop = FALSE]
dat <- dat[complete.cases(dat), , drop = FALSE]

if (nrow(dat) < 5) {
  stop("Too few samples remained after filtering.")
}

if (length(valid_predictors) < 2) {
  stop("Too few valid predictors remained after filtering.")
}

# ----------------------------- random forest model -----------------------------

rf_formula <- as.formula(paste(response_var, "~ ."))

set.seed(123)
rf_model <- randomForest(
  rf_formula,
  data = dat,
  importance = TRUE
)

saveRDS(
  rf_model,
  file = file.path(outdir, paste0(response_var, "_random_forest_model.rds"))
)

capture.output(
  rf_model,
  file = file.path(outdir, paste0(response_var, "_random_forest_model_summary.txt"))
)

# ----------------------------- observed vs predicted ---------------------------

rf_preds <- predict(rf_model, newdata = dat)

pred_df <- data.frame(
  Observed = dat[[response_var]],
  Predicted = rf_preds
)

write.csv(
  pred_df,
  file.path(outdir, paste0(response_var, "_observed_predicted.csv")),
  row.names = TRUE,
  quote = FALSE
)

p_pred <- ggplot(pred_df, aes(x = Observed, y = Predicted)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  stat_poly_eq(
    aes(label = paste(..adj.rr.label.., sep = "~~")),
    formula = y ~ x,
    parse = TRUE,
    size = 5,
    color = "black",
    label.x = 0.05,
    label.y = 0.95
  ) +
  labs(
    x = paste0("Observed ", response_var),
    y = paste0("Predicted ", response_var)
  ) +
  theme_bw() +
  theme(panel.grid = element_blank())

p_pred_marginal <- ggExtra::ggMarginal(p_pred, type = "histogram", fill = "transparent")

ggsave(
  file.path(outdir, paste0(response_var, "_observed_vs_predicted.pdf")),
  p_pred_marginal,
  width = 5.2,
  height = 5.0
)

# ----------------------------- variable importance -----------------------------

importance_mat <- rf_model$importance

importance_df <- tibble(
  Variable = rownames(importance_mat),
  IncMSE = importance_mat[, "%IncMSE"],
  IncNodePurity = importance_mat[, "IncNodePurity"]
) %>%
  arrange(desc(IncMSE))

write.csv(
  importance_df,
  file.path(outdir, paste0(response_var, "_variable_importance.csv")),
  row.names = FALSE,
  quote = FALSE
)

top_importance <- importance_df %>%
  slice_head(n = top_n)

write.csv(
  top_importance,
  file.path(outdir, paste0(response_var, "_top", top_n, "_variables.csv")),
  row.names = FALSE,
  quote = FALSE
)

p_imp_mse <- ggplot(top_importance, aes(x = reorder(Variable, IncMSE), y = IncMSE)) +
  geom_segment(aes(xend = Variable, y = 0, yend = IncMSE)) +
  geom_point(size = 3, alpha = 0.8) +
  coord_flip() +
  labs(x = NULL, y = "%IncMSE") +
  theme_bw() +
  theme(panel.grid.major.y = element_blank())

p_imp_purity <- ggplot(top_importance, aes(x = reorder(Variable, IncNodePurity), y = IncNodePurity)) +
  geom_segment(aes(xend = Variable, y = 0, yend = IncNodePurity)) +
  geom_point(size = 3, alpha = 0.8) +
  coord_flip() +
  labs(x = NULL, y = "IncNodePurity") +
  theme_bw() +
  theme(panel.grid.major.y = element_blank())

p_imp <- ggpubr::ggarrange(p_imp_mse, p_imp_purity, ncol = 2)

ggsave(
  file.path(outdir, paste0(response_var, "_top", top_n, "_importance.pdf")),
  p_imp,
  width = 10,
  height = 6
)

# ----------------------------- repeated rfcv -----------------------------------

set.seed(cv_seed)

x <- dat[, setdiff(colnames(dat), response_var), drop = FALSE]
y <- dat[[response_var]]

rf_cv <- replicate(
  cv_repeats,
  rfcv(x, y, cv.fold = cv_fold, step = cv_step),
  simplify = FALSE
)

rf_cv_df <- data.frame(sapply(rf_cv, "[[", "error.cv"))
rf_cv_df$Variable_number <- as.numeric(rownames(rf_cv_df))

rf_cv_long <- reshape2::melt(
  rf_cv_df,
  id.vars = "Variable_number",
  variable.name = "Repeat",
  value.name = "CV_error"
)

rf_cv_mean <- rf_cv_long %>%
  group_by(Variable_number) %>%
  summarise(
    Mean_CV_error = mean(CV_error, na.rm = TRUE),
    SD_CV_error = sd(CV_error, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  rf_cv_mean,
  file.path(outdir, paste0(response_var, "_rfcv_mean_error.csv")),
  row.names = FALSE,
  quote = FALSE
)

p_cv <- ggplot(rf_cv_mean, aes(x = Variable_number, y = Mean_CV_error)) +
  geom_line() +
  geom_point() +
  labs(x = "Number of variables", y = "Cross-validation error") +
  theme_bw() +
  theme(panel.grid = element_blank())

ggsave(
  file.path(outdir, paste0(response_var, "_rfcv_error.pdf")),
  p_cv,
  width = 5,
  height = 4
)

# ----------------------------- correlations for top variables ------------------

cor_results <- data.frame(
  Variable = character(),
  Correlation = numeric(),
  P_value = numeric(),
  stringsAsFactors = FALSE
)

for (v in top_importance$Variable) {
  cor_test <- cor.test(dat[[v]], dat[[response_var]], method = "pearson")
  cor_results <- rbind(
    cor_results,
    data.frame(
      Variable = v,
      Correlation = unname(cor_test$estimate),
      P_value = cor_test$p.value
    )
  )
}

cor_results$Significance <- cut(
  cor_results$P_value,
  breaks = c(0, 0.001, 0.01, 0.05, 1),
  labels = c("***", "**", "*", "NS"),
  include.lowest = TRUE
)

write.csv(
  cor_results,
  file.path(outdir, paste0(response_var, "_top", top_n, "_correlations.csv")),
  row.names = FALSE,
  quote = FALSE
)

p_cor_heat <- ggplot(
  cor_results,
  aes(x = factor(Variable, levels = top_importance$Variable), y = response_var)
) +
  geom_tile(aes(fill = Correlation), color = "white", linewidth = 0) +
  geom_text(
    aes(label = paste0(round(Correlation, 3), "\n", Significance)),
    color = "black",
    size = 4
  ) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Pearson\ncorrelation"
  ) +
  labs(x = "Variables", y = response_var) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

ggsave(
  file.path(outdir, paste0(response_var, "_top", top_n, "_correlation_heatmap.pdf")),
  p_cor_heat,
  width = 7,
  height = 3
)

# Scatterplots for top variables.
plot_df <- dat[, c(top_importance$Variable, response_var), drop = FALSE] %>%
  pivot_longer(
    cols = all_of(top_importance$Variable),
    names_to = "Variable",
    values_to = "Value"
  )

p_scatter <- ggplot(plot_df, aes(x = Value, y = .data[[response_var]])) +
  geom_point(alpha = 0.8, size = 2) +
  geom_smooth(method = "lm", se = TRUE) +
  facet_wrap(~ Variable, scales = "free_x") +
  labs(x = "Feature value", y = response_var) +
  theme_bw() +
  theme(panel.grid = element_blank())

ggsave(
  file.path(outdir, paste0(response_var, "_top", top_n, "_scatterplots.pdf")),
  p_scatter,
  width = 8,
  height = 6
)

# ----------------------------- session info -----------------------------------

sink(file.path(outdir, paste0(response_var, "_sessionInfo.txt")))
sessionInfo()
sink()

message("Random Forest analysis finished. Results are saved in: ", outdir)
