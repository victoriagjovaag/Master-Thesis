library(dplyr)
library(tidyr)

# Load the data
data <- read.csv("Data/Bootstrap/Lu_boot_dataset.csv", header=TRUE, stringsAsFactors=FALSE)

# Set the number of datasets to generate
num_datasets <- 100

# Set the directory to save the output files
output_dir <- "training_files/Bootstrap/Lu/42_nodes"

# Create the output directory if it doesn't already exist
if (!file.exists(output_dir)) {
  dir.create(output_dir)
}

# Loop through each dataset
for (i in 1:num_datasets) {
  
  # Randomly select five biological entities
  entities <- sample(data$Lu_nodes, 42)
  
  # Randomly select a value for each entity
  values <- sample(c(0, 1), 42, replace = TRUE)
  
  # Create a new dataset with only the selected entities and values
  new_data <- data %>%
    filter(Lu_nodes %in% entities) %>%
    mutate(value = ifelse(Lu_nodes %in% entities, values, NA)) %>%
    drop_na(value)
  
  # Set the output file name
  output_file <- paste0(output_dir, "/", "42_nodes_", i, ".txt")
  
  # Write the output file
  cat("Condition\n-\nResponse\n", file=output_file)
cat(paste(entities, ":", new_data$value, sep="", collapse="\t"), "\nWeight:1", file=output_file, append=TRUE)

}
