# visualize prediction performance of synergies 
# change location of "ss_hsa_file" for each observation 

suppressMessages(library(dplyr))
library(tibble)
library(emba)
library(usefun)
library(PRROC)
library(DT)

# Read ensemble-wise synergies file
# `ss` => models trained to steady state
ss_hsa_file = "ags_cascade_1.0/ags_cascade_1.0_many_ev_TP53_TCF_MYC/ags_cascade_1.0_ensemblewise_synergies.tab"
ss_hsa_ensemblewise_synergies = emba::get_synergy_scores(ss_hsa_file)


# Read observed synergies file
observed_synergies_file = 'ags_cascade_1.0/observed_synergies'
observed_synergies = emba::get_observed_synergies(observed_synergies_file)
# 1 (positive/observed synergy) or 0 (negative/not observed) for all tested drug combinations
observed = sapply(ss_hsa_ensemblewise_synergies$perturbation %in% observed_synergies, as.integer)

# Make a data table
pred_hsa = dplyr::bind_cols(ss_hsa_ensemblewise_synergies %>% rename(ss_score = score),
                            tibble::as_tibble_col(observed, column_name = "observed"))

# Get ROC statistics (`roc_res$AUC` holds the ROC AUC)
roc_res = usefun::get_roc_stats(df = pred_hsa, pred_col = "ss_score", label_col = "observed")


# Plot ROC
my_palette = RColorBrewer::brewer.pal(n = 9, name = "Set1")

plot(x = roc_res$roc_stats$FPR, y = roc_res$roc_stats$TPR,
     type = 'l', lwd = 3, col = my_palette[1], main = 'ROC curve, Ensemble-wise synergies (HSA)',
     xlab = 'False Positive Rate (FPR)', ylab = 'True Positive Rate (TPR)')
legend('bottomright', title = 'AUC', col = my_palette[1], pch = 19,
       legend = paste(round(roc_res$AUC, digits = 2), "Calibrated"), cex = 1.3)
grid(lwd = 0.5)
abline(a = 0, b = 1, col = 'lightgrey', lty = 'dotdash', lwd = 1.2)

# Get PR statistics (`pr_res$auc.davis.goadrich` holds the PR AUC)
# NOTE: PRROC considers by default that larger prediction values indicate the
# positive class labeling. For us, the synergy scores belonging to the positive
# or synergy class (observed = 1) are the lower ones, so we need to
# reverse the scores to correctly calculate the PR curve
pr_res = PRROC::pr.curve(scores.class0 = pred_hsa %>% pull(ss_score) %>% (function(x) {-x}),
                         weights.class0 = pred_hsa %>% pull(observed), curve = TRUE, rand.compute = TRUE)

plot(pr_res, main = 'PR curve, Ensemble-wise synergies (HSA)',
     auc.main = FALSE, color = my_palette[1], rand.plot = TRUE)
legend('topright', title = 'AUC', col = my_palette[1], pch = 19,
       legend = paste(round(pr_res$auc.davis.goadrich, digits = 2), "Calibrated"))
grid(lwd = 0.5)

