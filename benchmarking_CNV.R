
library(dplyr)


######## CONTROLS ########

# Import each coverage file and store it in the list with basename as the name
cov_dir_path <- "all_amplicon_coverage_files"
files <- list.files(path = cov_dir_path, pattern = "\\.txt$", full.names = TRUE)

list_cov_files <- list()

for (file in files) {

  file_name <- basename(file)
  list_cov_files[[file_name]] <- read.table(file, header = TRUE)
}


# Extract controls
extracted_samples <- data.frame(SampleID = character(), FromFile = character(), stringsAsFactors = FALSE)

for (file_name in names(list_cov_files)) {

  df <- list_cov_files[[file_name]]
  
  # Filter rows based on the SampleID patterns found in CNV controls DD2 and PM4
  extracted <- df[grepl("(?i)Dd2|PM-4|PM2_4", df$SampleID), ] # this for ALL controls: df[grepl("(?i)3D7|Dd2|PM|HB3", df$SampleID), ]
  
  if (nrow(extracted) > 0) {
    extracted$run <- file_name
    extracted_samples <- rbind(extracted_samples, extracted)
  }
}

controls <- extracted_samples %>%
  distinct(SampleID, run)

dim(controls)

controls$expected_CNV <- ifelse(grepl("(?i)Dd2", controls$SampleID), "MDR1", "")

# Add PM2 and PM3 to PM controls
pm_indices <- grepl("(?i)PM", controls$SampleID)
pm_rows <- controls[pm_indices, ]
pm_rows$expected_CNV <- "PM2"  # Assign "PM2" to one duplicated row
pm_rows2 <- pm_rows  # Create a copy of the duplicated rows
pm_rows2$expected_CNV <- "PM3"  # Assign "PM3" to the other duplicated row

# Combine the duplicated rows with original controls
controls <- rbind(controls, pm_rows, pm_rows2)
controls <- controls[controls$expected_CNV != "",]
controls <- controls[order(controls$run), ]


# exclude shit controls (pending...)



######## RESULTS ########

# Import each results file and store it in the list with basename as the name
files <- list.files(pattern = "\\_CNV_results.csv$", full.names = TRUE)

list_results_files <- list()

for (file in files) {
  
  file_name <- basename(file)
  list_results_files[[file_name]] <- read.csv(file, header = TRUE)
  
}


# Extract results
extracted_results <- data.frame(SampleID = character(), FromFile = character(), stringsAsFactors = FALSE)

for (file_name in names(list_results_files)) {
  
  df <- list_results_files[[file_name]]
  
  # Filter rows based on the SampleID patterns found in CNV controls DD2 and PM4
  extracted <- df[grepl("(?i)Dd2|PM-4|PM2_4", df$SampleID), ] # this for ALL controls: df[grepl("(?i)3D7|Dd2|PM|HB3", df$SampleID), ]
  
  if (nrow(extracted) > 0) {
    extracted$run <- file_name
    extracted_results <- rbind(extracted_results, extracted)
  }
}

extracted_results <- extracted_results %>%
  rename(observed_CNV = loi)

results <- extracted_results %>%
  distinct(SampleID, run, observed_CNV)

dim(results)



######## BENCHMARKING ########

# Merge expected and observed CNV for the CNV controls
controls$run <- gsub("_RESULTS_v0\\.1\\.8_amplicon_coverage\\.txt", "", controls$run)
results$run <- gsub("_RESULTS_v0.1.8_amplicon_coverage_CNV_results.csv", "", results$run)

merged_data <- merge(controls, results, by = c("SampleID", "run"), all =T)

#deal with controls that have multiple CNV like PM
merged_data <- merged_data[merged_data$expected_CNV == merged_data$observed_CNV | is.na(merged_data$observed_CNV),]

# da numbas demselves
total <- length(merged_data$expected_CNV)
TP <- sum(!is.na(merged_data$expected_CNV == merged_data$observed_CNV)) # Calculate True Positives (TP)
FP <- sum(!is.na(merged_data$observed_CNV) & is.na(merged_data$expected_CNV)) # Calculate False Positives (FP)
FN <- sum(is.na(merged_data$observed_CNV) & !is.na(merged_data$expected_CNV)) # Calculate False Negatives (FN)

# evaluation metrics
precision <- TP / (TP + FP)
recall <- TP / (TP + FN)
f1_score <- 2 * precision * recall / (precision + recall)

print(paste("Precision:", precision))
print(paste("Recall:", recall))
print(paste("F1-score:", f1_score))
