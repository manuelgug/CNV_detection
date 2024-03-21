library(ggplot2)
library(ggpubr)
library(dplyr)


data <- read.table("CNV_runs_sample_coverage/REACT2_NextSeq01_amplicon_coverage_DD2.txt", header = T)
data <- data[,-3:-4]

#remove neg controls and undetermined
data <- data[!grepl("(?i)neg", data$SampleID), ]
data <- data[!grepl("(?i)Undetermined", data$SampleID), ]

#samples to be flagged as low quality (below a median of 100 read count)
co <- data %>%
  group_by(SampleID)%>%
  summarize(median_read_count = median(OutputPostprocessing))

samples_below_median_100 <- co[co$median_read_count < 100,]$SampleID 


# calculate proportions of amplicons for each sample
data_norm <- data %>% 
  group_by(SampleID) %>%
  summarize(NORM_OutputPostprocessing = OutputPostprocessing / sum(OutputPostprocessing))

data_norm <- as.data.frame(cbind(data_norm, Locus = data$Locus))


# add loci of interest column
loci_of_interest <- readRDS("loci_of_interest.RDS")

loi <- character(nrow(data))

# Loop through each row of 'data'
for (i in seq_len(nrow(data_norm))) {
  # Check if 'Locus' matches any element in 'loci_of_interest'
  matching_loi <- sapply(loci_of_interest, function(x) any(grepl(data_norm$Locus[i], x)))
  # If there's a match, assign the name of the matching element to 'loi'
  if (any(matching_loi)) {
    loi[i] <- names(loci_of_interest)[which(matching_loi)]
  } else {
    loi[i] <- NA  # If no match is found, assign "Not Found"
  }
}

data_norm$loi <- loi


#extract single-copy controls
controls <- data_norm[!grepl("(?i)Dd2|PM|HB3", data_norm$SampleID) & grepl("(?i)3D7", data_norm$SampleID), ]
controls_nozero <- controls[controls$NORM_OutputPostprocessing != 0, ]

#calculate thresholds from single-copy controls
thresholds_controls<- controls_nozero %>% 
  group_by(SampleID) %>%
  summarize(q_up = quantile(NORM_OutputPostprocessing,0.99),
            q_down = quantile(NORM_OutputPostprocessing,0.025),
            max = max(NORM_OutputPostprocessing),
            min = min(NORM_OutputPostprocessing))


# ##############################################################
#find best control (most reads total)
# controls_raw <- data[!grepl("(?i)Dd2|PM|HB3", data$SampleID) & grepl("(?i)3D7", data$SampleID), ]
# 
# total_read_counts <- controls_raw %>%
#   group_by(SampleID) %>%
#   summarize(read_count_total = sum(OutputPostprocessing))
# 
# best_control <- total_read_counts[which.max(total_read_counts$read_count_total),][1]
# 
# 
# #normalize proportions using best control
# best_control_propotions <- controls_nozero[controls_nozero$SampleID ==best_control$SampleID,]
# 
# best_control_propotions
# data_norm
# 
# merged_table<- merge(data_norm, best_control_propotions[c("Locus", "NORM_OutputPostprocessing")], by ="Locus")
# 
# merged_table$fold_change_probs<- merged_table$NORM_OutputPostprocessing.x / merged_table$NORM_OutputPostprocessing.y
##############################################################

# VISUALIZATION #

#grab a sample
unique_samples <- sort(unique(merged_table$SampleID))

# Calculate the median of OutputPostprocessing
up_threshold_q <- mean(thresholds_controls$q_up)
down_threshold_q <- mean(thresholds_controls$q_down)
up_threshold <- max(thresholds_controls$max)
down_threshold <- min(thresholds_controls$min)


# Initialize an empty list to store plots
plot_list <- list()

# Loop through each sample in unique_samples
for (xsample in unique_samples) {
  # Subset the data for the current sample
  sample <- data_norm[data_norm$SampleID == xsample, ]
  
  # Create scatterplot
  p <- ggplot(sample, aes(x = reorder(Locus, -NORM_OutputPostprocessing), y = NORM_OutputPostprocessing, color = loi)) +
    geom_point(size=4, alpha =0.7) +
    geom_hline(yintercept = up_threshold, linetype = "dashed", color = "red") +
    #geom_hline(yintercept = 2, linetype = "dashed", color = "red") +
    labs(x = "Locus", y = "Read Proportions", title = xsample) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6.5))
  
  # Add the plot to the list
  plot_list[[length(plot_list) + 1]] <- p
}

# Arrange plots in a grid
grid <- ggarrange(plotlist = plot_list[[1:9]], nrow = 3, ncol =3 )



# ACTUAL RESULTS #

# quantitative
data_norm$CNV_result <- ifelse(data_norm$NORM_OutputPostprocessing > max(thresholds_controls$max), "CNV", "single")
data_norm$QC<- ifelse(data_norm$SampleID %in% samples_below_median_100, "median_below_100", "good")

#check results: good QC samples that have CNV for loci of interest
results <- data_norm[!is.na(data_norm$loi)& data_norm$CNV_result == "CNV" & data_norm$QC == "good",]

length(unique(results$SampleID))

# # quantitative
# merged_table$CNV_result <- ifelse(merged_table$fold_change_probs > 2, "CNV", "single")
# merged_table$QC<- ifelse(merged_table$SampleID %in% samples_below_median_100, "median_below_100", "good")
# 
# #check results: good QC samples that have CNV for loci of interest
# results <- merged_table[!is.na(merged_table$loi)& merged_table$CNV_result == "CNV" & merged_table$QC == "good",]



# NOTES
# TEST THIS METHOD WITH ALL CONTROLS FROM ALL RUNS
# BENCHMARK AGAINST estCNV
# ASSESS FALSE POSITIVES USING quantile99 vs max
