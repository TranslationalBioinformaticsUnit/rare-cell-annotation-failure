# ==============================================================================
# Mixed-effects logistic regression that predicts poor annotation of rare cell
# types from (i) proximity to the nearest common cell type and (ii) class
# imbalance, with SHAP attribution and leave-one-dataset-out cross-validation.
#
# Manuscript: "Rare cell type annotation exposes predictable failure modes of
#              representation learning in single-cell foundation models"
#
# INPUT   data/LR_data.csv  (one row per rare cell-type instance)
#   Case                     analysis case, e.g. "CellPLM_embeddings"
#   Dataset                  dataset identifier (random intercept)
#   low_f1                   TRUE if the F1-score of the rare cell type is < 0.5
#                            ("poorly annotated"), otherwise FALSE
#   Euclidean_min            proximity: Euclidean distance between the centroid of
#                            the rare cell type and that of its nearest common
#                            cell type
#   Imbalance_Ratio_Log_max  class imbalance: largest log-transformed imbalance
#                            ratio (common / rare cells) over all common cell
#                            types.
#
#
# USAGE   Run from the repository root:   Rscript analysis.R
#         To analyze another case, change CASE below.
#
# TESTED  R 4.4.0 | dplyr 1.1.4 | pROC 1.18.5 | PRROC 1.4 | lme4 1.1-36 |
#         fastshap 0.1.1
# ==============================================================================


# ---- 0. Settings -------------------------------------------------------------
DATA_FILE  <- file.path("data", "LR_data.csv")
CASE       <- "CellPLM_embeddings"  # value of the `Case` column to analyze
N_BOOT     <- 1000                  # bootstrap resamples for the LODO AUROC CI
SEED       <- 123                    # one seed for the whole script (SHAP, bootstrap)


# ---- 1. Packages -------------------------------------------------------------
# install.packages(c("dplyr", "pROC", "PRROC", "lme4", "fastshap"))
required <- c("dplyr", "pROC", "PRROC", "lme4", "fastshap")
missing  <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop("Please install the missing packages: ", paste(missing, collapse = ", "))
}

library(dplyr)
library(pROC)    # roc(), auc()
library(PRROC)   # pr.curve()
library(lme4)    # glmer(), glmerControl()
# fastshap is called as fastshap::explain()

set.seed(SEED)


# ---- 2. Load data ------------------------------------------------------------
if (!file.exists(DATA_FILE)) {
  stop("Input file not found: ", DATA_FILE,
       "\nRun the script from the repository root and place LR_data.csv in data/.")
}
dat_all <- read.csv(DATA_FILE)

required_cols <- c("Case", "Dataset", "low_f1", "Euclidean_min",
                   "Imbalance_Ratio_Log_max")
stopifnot(all(required_cols %in% names(dat_all)))

cat("Available cases:", paste(unique(dat_all$Case), collapse = ", "), "\n")
stopifnot(CASE %in% dat_all$Case)

dat <- dat_all[dat_all$Case == CASE, ]
dat$low_f1 <- as.logical(dat$low_f1)     # TRUE = poorly annotated (F1 < 0.5)

cat(sprintf("Case: %s | rare cell-type instances: %d | poorly annotated: %d | datasets: %d\n",
            CASE, nrow(dat), sum(dat$low_f1), length(unique(dat$Dataset))))


# ---- 3. Final model: Euclidean_min + Imbalance_Ratio_Log_max -----------------
formula_final <- low_f1 ~ Euclidean_min + Imbalance_Ratio_Log_max + (1 | Dataset)

glmm_final <- glmer(formula_final, data = dat, family = binomial)
print(summary(glmm_final))

# In-sample predictions (include the dataset random intercept)
dat$pred_prob <- predict(glmm_final, newdata = dat, type = "response")

# AUROC
roc_in   <- roc(dat$low_f1, dat$pred_prob,
                levels = c(FALSE, TRUE), direction = "<", quiet = TRUE)
auroc_in <- as.numeric(auc(roc_in))

# AUPRC (preferred over AUROC when the poor class is rare)
pr_in    <- pr.curve(scores.class0 = dat$pred_prob[dat$low_f1],
                     scores.class1 = dat$pred_prob[!dat$low_f1])
auprc_in <- pr_in$auc.integral

cat(sprintf("In-sample AUROC = %.3f | AUPRC = %.3f\n", auroc_in, auprc_in))


# ---- 4. Feature contributions (SHAP) -----------------------------------------
# Monte Carlo Shapley values (fastshap) for the fixed-effects predictions,
# on the probability scale (random effects set to zero via re.form = NA).
pred_wrapper <- function(object, newdata) {
  predict(object, newdata = newdata, type = "response", re.form = NA)
}

X <- dat %>% dplyr::select(Euclidean_min, Imbalance_Ratio_Log_max)

shap_values <- fastshap::explain(
  object        = glmm_final,
  feature_names = colnames(X),
  X             = X,
  pred_wrapper  = pred_wrapper,
  nsim          = 1000
)

shap_importance <- data.frame(
  Feature     = colnames(shap_values),
  MeanAbsSHAP = colMeans(abs(shap_values))
)
shap_importance$PercentContribution <-
  100 * shap_importance$MeanAbsSHAP / sum(shap_importance$MeanAbsSHAP)
shap_importance <- shap_importance[order(-shap_importance$PercentContribution), ]

print(shap_importance)


# ---- 5. Leave-one-dataset-out cross-validation (LODO-CV) ---------------------
# Each dataset is held out in turn; the model is refitted on the remaining
# datasets and predicts the held-out dataset (new random-effect level allowed).
# Out-of-fold predictions are pooled across datasets.
lodo_control    <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
datasets        <- unique(dat$Dataset)
oof_predictions <- list()

for (test_ds in datasets) {
  
  train_data <- dat %>% filter(Dataset != test_ds)
  test_data  <- dat %>% filter(Dataset == test_ds)
  
  if (length(unique(train_data$low_f1)) < 2) {
    cat("Skipping", test_ds, "- training set has only one class\n")
    next
  }
  
  model_cv <- tryCatch(
    glmer(formula_final, data = train_data, family = binomial, control = lodo_control),
    error = function(e) NULL
  )
  if (is.null(model_cv)) {
    cat("Model failed for", test_ds, "\n")
    next
  }
  
  oof_predictions[[length(oof_predictions) + 1]] <- data.frame(
    Dataset    = test_ds,
    true_label = test_data$low_f1,
    pred_prob  = predict(model_cv, newdata = test_data, type = "response",
                         allow.new.levels = TRUE)
  )
}

oof_df <- bind_rows(oof_predictions)

cat("\nOut-of-fold instances:", nrow(oof_df), "\n")
cat("Poorly annotated:", sum(oof_df$true_label), "\n")
cat("Accurately annotated:", sum(!oof_df$true_label), "\n")
cat("Datasets included:", length(unique(oof_df$Dataset)), "\n\n")

# Pooled AUROC and AUPRC
roc_oof   <- roc(oof_df$true_label, oof_df$pred_prob,
                 levels = c(FALSE, TRUE), direction = "<", quiet = TRUE)
auroc_oof <- as.numeric(auc(roc_oof))

cat(sprintf("LODO-CV AUROC = %.3f ", auroc_oof))


# ---- 6. Bootstrap 95% CI of the pooled LODO AUROC ----------------------------
# Instances are resampled with replacement; percentile interval.
boot_auroc <- replicate(N_BOOT, {
  boot_df <- oof_df[sample(nrow(oof_df), replace = TRUE), ]
  if (length(unique(boot_df$true_label)) < 2) return(NA_real_)
  as.numeric(auc(roc(boot_df$true_label, boot_df$pred_prob,
                     levels = c(FALSE, TRUE), direction = "<", quiet = TRUE)))
})

ci_auroc <- quantile(boot_auroc, c(0.025, 0.975), na.rm = TRUE)
cat(sprintf("LODO-CV AUROC 95%% CI (bootstrap): %.3f-%.3f\n", ci_auroc[1], ci_auroc[2]))

