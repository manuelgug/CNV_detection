library(ggplot2)
library(ggpubr)
library(dplyr)
library(cowplot)


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

data_filtered <- data[!data$SampleID %in% samples_below_median_100, ]


#exclude one HRP3 amplicon wat too outlierish:
#data_filtered <- data_filtered[!data_filtered$Locus =="Pf3D7_13_v3-2816310-2816554-1B",]


# calculate proportions of amplicons for each sample
data_norm <- data_filtered %>% 
  group_by(SampleID) %>%
  summarize(NORM_OutputPostprocessing = OutputPostprocessing / sum(OutputPostprocessing))

data_norm <- as.data.frame(cbind(data_norm, Locus = data_filtered$Locus))

unique_samples <- sort(unique(data_norm$SampleID))



# samples to be flagged as low quality by slope
sample_slopes <- data.frame(SampleID = character(), Slope = numeric(), stringsAsFactors = FALSE)

# Loop through each sample in unique_samples
for (sample_id in unique_samples) {
  # Subset the data for the current sample
  sample <- data_norm[data_norm$SampleID == sample_id, ]
  
  # Sort and preprocess the data as needed
  sample_slope <- as.data.frame(cbind(x = 1:length(sample$NORM_OutputPostprocessing), y = sort(sample$NORM_OutputPostprocessing, decreasing = TRUE)))
  sample_slope <- sample_slope[-c((length(sample_slope$y) - 100):length(sample_slope$y)), ]  # Remove 100 least present amplicons
  sample_slope <- sample_slope[-c(1:100), ]  # Remove 50 most and least present amplicons
  
  # Fit a linear regression model
  model <- lm(y ~ x, data = sample_slope)
  
  # Extract the slope coefficient
  slope <- coef(model)[[2]]
  
  # Add the sample ID and slope to the data frame
  sample_slopes <- rbind(sample_slopes, data.frame(SampleID = sample_id, Slope = slope))
}

hist(sample_slopes$Slope, breaks = 150)

# Calculate quantiles
q_up <- quantile(sample_slopes$Slope, probs = 0.99) #samples with many amplicons at 0 (MAYBE HERE'S NO NEED TO REMOVE THESE?)
q_down <- quantile(sample_slopes$Slope, probs = 0.01) #samples with overall low yield or 

# Add lines for quantiles 99 and 1
abline(v = q_down, col = "red", lty = 2) 
abline(v = q_up, col = "blue", lty = 2) 

sample_slopes <- sample_slopes[order(sample_slopes$Slope), ]

bad_slope_samples <- sample_slopes[sample_slopes$Slope > q_up | sample_slopes$Slope < q_down,]$SampleID



# add loci of interest column
loci_of_interest <- readRDS("loci_of_interest.RDS")

loi <- character(nrow(data_filtered))

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
  summarize(max = max(NORM_OutputPostprocessing),
            min = min(NORM_OutputPostprocessing))

# top 5 most abundant amplicons in controls
ranking_amps <- controls_nozero %>%
  group_by(Locus) %>%
  summarize(mean_NORM_OutputPostprocessing = mean(NORM_OutputPostprocessing))

ranking_amps <- ranking_amps%>%
  arrange(desc(mean_NORM_OutputPostprocessing))

# remove top 5 mos abundant amplicons from data_norm (TO AVOID FALSE POSITIVES)
data_norm <- data_norm[!data_norm$Locus %in% ranking_amps$Locus[1:5],]


# ##############################################################
# #find best control (most reads total)
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

# Set thresholds
up_threshold <- max(thresholds_controls$max)
up_threshold_mean <- mean(thresholds_controls$max)


plot_list <- list()

# Loop through each sample in unique_samples
for (xsample in unique_samples) {
  # Subset the data for the current sample
  sample <- data_norm[data_norm$SampleID == xsample, ]
  
  # Create scatterplot
  p <- ggplot(sample, aes(x = reorder(Locus, -NORM_OutputPostprocessing), y = NORM_OutputPostprocessing, color = loi)) +
    geom_point(size=4, alpha =0.7) +
    geom_hline(yintercept = up_threshold_mean, linetype = "dashed", color = "red") +
    #geom_hline(yintercept = 2, linetype = "dashed", color = "red") +
    labs(x = "Locus", y = "Read Proportions", title = xsample) +
    #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6.5))+
    theme(axis.text.x = element_blank(),  # Remove amplicon names
          axis.title.x = element_blank())
  
  # Add the plot to the list
  plot_list[[length(plot_list) + 1]] <- p
}

# Assuming 'plot_list' contains a list of ggplot objects
combined_plot <- plot_grid(plotlist = plot_list, nrow = round(sqrt(length(unique_samples))), ncol = round(sqrt(length(unique_samples))))
ggsave("grid_of_plots.pdf", plot = combined_plot, width = 60, height = 50, dpi = 300, limitsize = FALSE)

# ACTUAL RESULTS #

# quantitative
data_norm$CNV_result <- ifelse(data_norm$NORM_OutputPostprocessing > up_threshold_mean, "CNV", "single")
data_norm$QC<- ifelse(data_norm$SampleID %in% samples_below_median_100, "bad", "good")
data_norm$slope<- ifelse(data_norm$SampleID %in% bad_slope_samples, "bad", "good")

#check results: good QC samples that have CNV for loci of interest
results <- data_norm[!is.na(data_norm$loi) & 
                       data_norm$CNV_result == "CNV" & 
                       data_norm$QC == "good" &
                       data_norm$slope == "good",]

#remove single-copy controls (since i used the mean of the max value of single-copy controls, it's expected that some end up as false positives. it's ok)
results <- results[!(grepl("(?i)3D7", results$SampleID) & !grepl("(?i)(Dd2|PM|HB3)", results$SampleID)), ]
results

print(paste0("There are ", length(unique(results$SampleID)), " samples with good QC that have CNV for loci of interest out of ", length(unique_samples)))

write.csv(results, "CNV_results.csv", row.names = F)



# NOTES
#   one hrp2 amplicon often comes up as cnv whe it's probably not. remove it?
# esclude consistently high amplicons across samples from the run to avoid false positives? (does it affect true positives? how much?)
# TEST THIS METHOD WITH ALL CONTROLS FROM ALL RUNS
# BENCHMARK AGAINST estCNV
# ASSESS FALSE POSITIVES USING quantile99 vs max
