# importing necessary packages
library(tidyverse)
library(here)
library(workflowr)
library(rmarkdown)
library(arsenal)
library(dbplyr) 
library(tidyselect)
source("path_to_node_name_conversion_file/node_name_conversion.r") # select approperiate model file


# load PROFILE output 
PROFILE_nodes <- read_csv("/path_to_output/PROFILE_output.csv") %>% # importing the nodes from PROFILE output
as.data.frame()

#binarising the normalised data 
#PROFILE_nodes$ex.value <- ifelse(PROFILE_nodes$ex.value > 0.5, "1", "0")
#PROFILE_nodes

#selecting relevant columns
PROFILE_nodes <- dplyr::select(PROFILE_nodes, nodes, ex.value)

# converting node names to match node names in the model network
converted_node_names <- node_name_conversion(df = PROFILE_nodes)
converted_node_names

# removing nodes that are not in the model
model_nodes <- read_csv("path_to_file_with_model_nodes/model_nodes.csv") %>% # importing the nodes in the model
as.data.frame()

filtered_nodes_list <- intersect(PROFILE_nodes$nodes, model_nodes$model_nodes) %>%
as.list()

nodes_for_training_data <- PROFILE_nodes[PROFILE_nodes$nodes %in% filtered_nodes_list, ]
nodes_for_training_data

 # format df to training data format
  cat("Condition",file="training_PROFILE.txt",sep="\n")
  cat("-",file="training_PROFILE.txt",sep="\n", append = TRUE)
  cat("Response",file="training_PROFILE.txt",sep="\n", append = TRUE)
  for (i in 1:nrow(PROFILE_nodes)) {
  node_value_string = paste0(PROFILE_nodes$nodes[i], ":", PROFILE_nodes$ex.value[i])
  cat(node_value_string, file="training_PROFILE.txt", append=TRUE)
  cat("\t", file="/training_PROFILE.txt", append=TRUE)
  }
  cat("\nWeight:1",file="training_PROFILE.txt",sep="\n", append = TRUE)
  
  