
# importing necessary packages
library(tidyverse)
library(here)
library(OmnipathR)
library(dorothea)
library(decoupleR)
library(workflowr)
library(rmarkdown)
library(org.Hs.eg.db)
library(arsenal)
library(dbplyr)
library(tidyselect)
source("/node_conversion/node_name_conversion_Park.r") #loading file for node name conversion from HGNC

# importing the network (dorothea regulons)
dorothea_df <- dorothea_hs %>%
  dplyr::filter(confidence %in% c("A", "B", "C")) %>% #using confidence levels A, B and C
  dplyr::select(target, tf, mor) %>%
  as.data.frame()
dorothea_df$likelihood <- 1

# importing the matrix with molecular features
# changes done in csv file from CMP: replaced spaces in headers to underscores and removed all ""
COLO_differential_analysis <- read_csv("data/CMP_COLO-205_rnaseq.csv") %>%
  tibble::column_to_rownames("Symbol") %>%
  dplyr::select(TPM_value) %>%
  as.matrix()

# finding NAs AND Infs
 COLO_differential_analysis[is.na(COLO_differential_analysis),] # 0 - no need for removal

# excluding NAs from the df (if necessary)
# COLO_differential_analysis_without_Nas <- na.omit(COLO_differential_analysis)

# running decopleR: consensus
decoupleR <- decouple(mat= COLO_differential_analysis, network = dorothea_df, .source='tf', 
.target='target', minsize = 5)                                   
as.data.frame

results_consensus <- run_consensus(decoupleR, include_time=FALSE, seed = NULL)
results_consensus

#filtering for significant p-values (using threshould of 0.05)
COLO_gitsbe_df <- subset(results_consensus, p_value < 0.2)
COLO_gitsbe_df

#scaling the data 
min <- min(COLO_gitsbe_df$score)
max <- max(COLO_gitsbe_df$score)

scaled_data <- COLO_gitsbe_df

for (i in 1:ncol(scaled_data)) {
    if(grepl("score", colnames(scaled_data)[i])) {
      scaled_data[ , i] <- (scaled_data[ , i] - min) / (max - min)
        }
    }

#binarizing the scaled data (using threshold of 0.5)
scaled_data$score <- ifelse(scaled_data$score > 0.5, "1", "0")
scaled_data

#selecting relevant columns
COLO_gitsbe_relevant <- dplyr::select(scaled_data, source, score)

# converting node names to match node names in Park
converted_node_names <- node_name_conversion(df = COLO_gitsbe_relevant)
converted_node_names

# removing nodes that are not in Park network
Park_nodes <- read_csv("node_conversion/Park_nodes.csv") %>% 
as.data.frame()

filtered_nodes_list <- intersect(converted_node_names$source, Park_nodes$Park_nodes) %>%
as.list()

nodes_for_training_data <- converted_node_names[converted_node_names$source %in% filtered_nodes_list, ]
nodes_for_training_data

# format df to training data format
  cat("Condition",file="training_files/COLO-205/Park/training_file.txt",sep="\n")
  cat("-",file="training_files/COLO-205/Park/training_file.txt",sep="\n", append = TRUE)
  cat("Response",file="training_files/COLO-205/Park/training_file.txt",sep="\n", append = TRUE)
  for (i in 1:nrow(nodes_for_training_data)) {
  node_value_string = paste0(nodes_for_training_data$source[i], ":", nodes_for_training_data$score[i])
  cat(node_value_string, file="training_files/COLO-205/Park/training_file.txt", append=TRUE)
  cat("\t", file="training_files/COLO-205/Park/training_file.txt", append=TRUE)
  }
  cat("\nWeight:1",file="training_files/COLO-205/Park/training_file.txt",sep="\n", append = TRUE)