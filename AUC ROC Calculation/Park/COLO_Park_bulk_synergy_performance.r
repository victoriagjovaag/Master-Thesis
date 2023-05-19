suppressMessages(library(dplyr))
library(tibble)
library(emba)
library(usefun)
library(PRROC)
library(DT)
library(ggplot2)

# Specify the path to the main directory containing the subdirectories with ensemble-wise synergy files
main_directory <- "COLO_205_Park/Synergy_results/Bootstrap/COLO_205_Park_58_nodes"

# Identify all subdirectories within the main directory
subdirectories <- list.dirs(path = main_directory, full.names = TRUE, recursive = FALSE)

# Identify all ensemble-wise synergy files within the subdirectories
synergy_files <- list.files(path = subdirectories, pattern = "ensemblewise_synergies.tab", full.names = TRUE, recursive = TRUE)

# Initialize a list to store the AUC results for each ensemble-wise synergy file
auc_results <- list()

# Loop over each ensemble-wise synergy file and calculate the AUC
for (synergy_file in synergy_files) {
  # Get the name of the directory containing the file
  directory <- dirname(synergy_file)
  # Extract the name of the ensemble-wise synergy file
  synergy_name <- gsub("ensemblewise_synergies.tab", "", basename(synergy_file))
  # Read the ensemble-wise synergies file
  ss_ensemblewise_synergies <- emba::get_synergy_scores(synergy_file)
  # Read observed synergies file
  observed_synergies_file <- file.path("COLO_205_Park/observed_synergies_Park_COLO_205")
  observed_synergies <- emba::get_observed_synergies(observed_synergies_file)
  # 1 (positive/observed synergy) or 0 (negative/not observed) for all tested drug combinations
  observed <- sapply(ss_ensemblewise_synergies$perturbation %in% observed_synergies, as.integer)
  # Make a data table
  pred_hsa <- dplyr::bind_cols(ss_ensemblewise_synergies %>% rename(ss_score = score),
                               tibble::as_tibble_col(observed, column_name = "observed"))
  # Get ROC statistics
  roc_res <- usefun::get_roc_stats(df = pred_hsa, pred_col = "ss_score", label_col = "observed")
  # Get PR statistics
  pr_res <- PRROC::pr.curve(scores.class0 = pred_hsa %>% pull(ss_score) %>% (function(x) {-x}),
                            weights.class0 = pred_hsa %>% pull(observed), curve = TRUE, rand.compute = TRUE)
  # Store the AUC results for this ensemble-wise synergy file
  auc_results[[synergy_name]] <- c(roc_res$AUC, pr_res$auc.davis.goadrich)
}

# Combine the AUC results into a data frame
auc_df <- data.frame(synergy_name = names(auc_results), 
                     auc_roc = sapply(auc_results, `[`, 1),
                     auc_pr = sapply(auc_results, `[`, 2))

# Print the AUC results
print(auc_df)

# Plot a histogram of the AUC ROC values
ggplot(auc_df, aes(x = auc_roc)) +
  geom_histogram(bins = 20, fill = "skyblue", color = "black") +
  xlab("AUC ROC") +
  ylab("Frequency")

  # Plot a histogram of the AUC PR values
ggplot(auc_df, aes(x = auc_pr)) +
  geom_histogram(bins = 20, fill = "skyblue", color = "black") +
  xlab("AUC PR") +
  ylab("Frequency")

# Save the AUC results to a CSV file with columns separated by commas
#write.table(auc_df, file = "Bulk_synergy_results/COLO-205_Park_PROFILE_inverted.csv", sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Calculate the mean, median and standard deviation of  AUC results
mean_auc_roc <- mean(auc_df$auc_roc)
median_auc_roc <- median(auc_df$auc_roc)
sd_auc_roc <- sd(auc_df$auc_roc)
mean_auc_pr <- mean(auc_df$auc_pr)
median_auc_pr <- median(auc_df$auc_pr)


# Print the mean, median and standard deviation of AUC results
cat("Mean AUC ROC:", mean_auc_roc, "\n")
cat("Median AUC ROC:", median_auc_roc, "\n")
cat("Standard deviation AUC ROC:", sd_auc_roc, "\n")
cat("Mean AUC PR:", mean_auc_pr, "\n")
cat("Median AUC PR:", median_auc_pr, "\n")

# Calculate the 95% confidence interval for the mean AUC ROC using a t-test
ci_auc_roc<- t.test(auc_df$auc_roc, conf.level = 0.95)$conf.int
print(ci_auc_roc)

# Calculate the 95% confidence interval for the mean AUC PR using a t-test
ci_auc_pr <- t.test(auc_df$auc_pr, conf.level = 0.95)$conf.int

# Print the confidence interval for the AUC results
cat(sprintf("95%% Confidence Interval for the Mean AUC ROC: [%.3f, %.3f]\n", ci_auc_roc[1], ci_auc_roc[2]))
cat(sprintf("95%% Confidence Interval for the Mean AUC PR: [%.3f, %.3f]\n", ci_auc_pr[1], ci_auc_pr[2]))

# Example AUC ROC value
auc_value <- 0.40

# Calculate the z-score for the AUC value
z_score <- (auc_value - mean_auc_roc) / (sd_auc_roc)

# Calculate the p-value
p_value <- 2* (1 - pnorm(abs(z_score)))


# Print the p-value
cat("The p-value for an AUC of", auc_value, "given a mean AUC of", mean_auc_roc, "and a standard deviation of", sd_auc_roc, "is", p_value, "\n")