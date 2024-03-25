library(ggplot2)
library(ggpubr)
library(dplyr)
library(cowplot)


###
# parameters to tune:

# 1) amount of amplicons to keep for the slope analysis
# sample_slope[(length(sample_slope$y) - 29):length(sample_slope$y), ]

# 2) quantiles to use as thresholds
# q_up <- quantile(sample_slopes$Slope, probs = 0.95) 
# q_down <- quantile(sample_slopes$Slope, probs = 0.05)

# 3) min or mean max proportion as threshold from single-copy controls
# up_threshold <- min(thresholds_controls$max)

# 4) amplicons to exclude

####

data <- read.table("CNV_runs_sample_coverage/SMC2_NextSeq01_amplicon_coverage_DD2.txt", header = T)
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
  
  #exclude amps that were not amplified
  sample <- sample[sample$NORM_OutputPostprocessing != 0,]
  
  # Sort and preprocess the data as needed
  sample_slope <- as.data.frame(cbind(x = 1:length(sample$NORM_OutputPostprocessing), y = sort(sample$NORM_OutputPostprocessing, decreasing = TRUE)))
  #sample_slope <- sample_slope[-c((length(sample_slope$y) - 20):length(sample_slope$y)), ]  # Remove 100 least present amplicons
  #sample_slope <- sample_slope[-c(1:50), ]  # Remove 50 most present amplicons
  
  # keep 50 least abundant amplicons (most informative)
  sample_slope <- sample_slope[(length(sample_slope$y) - 29):length(sample_slope$y), ]
  
  # Fit a linear regression model
  model <- lm(y ~ x, data = sample_slope)
  
  # Extract the slope coefficient
  slope <- coef(model)[[2]]
  
  # Add the sample ID and slope to the data frame
  sample_slopes <- rbind(sample_slopes, data.frame(SampleID = sample_id, Slope = slope))
}

#test for normality
shapiro_test <- shapiro.test(sample_slopes$Slope)
p_value <- shapiro_test$p.value #not significant means normal means solid run

# # Calculate quantiles
# q_up <- quantile(sample_slopes$Slope, probs = 0.99) #samples with many amplicons at 0 (MAYBE HERE'S NO NEED TO REMOVE THESE?)
# q_down <- quantile(sample_slopes$Slope, probs = 0.05) #samples with overall low yield or 

# # calculate sds
 sd_down <- mean(sample_slopes$Slope) - sd(sample_slopes$Slope)*2
 sd_up <- mean(sample_slopes$Slope) + sd(sample_slopes$Slope)*2

histogram_slopes <- ggplot(sample_slopes, aes(x = Slope)) +
  geom_histogram(bins = 150, fill = "skyblue", color = "skyblue", alpha = 0.7) +
  labs(x = "Slope", y = "Frequency", title = "") +
  theme_minimal()+
  annotate("text", x = Inf, y = Inf, label = paste("Shapiro-Wilk p-value:", round(p_value, 4)), 
           hjust = 1, vjust = 1, size = 4)+
  geom_vline(xintercept = sd_up, col = "blue", linetype = "dashed") +
  geom_vline(xintercept = sd_down, col = "red", linetype = "dashed")

histogram_slopes

sample_slopes <- sample_slopes[order(sample_slopes$Slope), ]

bad_slope_samples <- sample_slopes[sample_slopes$Slope > sd_up | sample_slopes$Slope < sd_down,]$SampleID



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

# top 5 most abundant amplicons in controls
ranking_amps <- controls_nozero %>%
  group_by(Locus) %>%
  summarize(mean_NORM_OutputPostprocessing = mean(NORM_OutputPostprocessing))

ranking_amps <- ranking_amps%>%
  arrange(desc(mean_NORM_OutputPostprocessing))

# remove top 3 most abundant amplicons from controls_nozero (TO AVOID FALSE POSITIVES)
controls_nozero <- controls_nozero[!controls_nozero$Locus %in% ranking_amps$Locus[1:5],]

#exclude some HRP3 and PM2 ampliconS that too outlierish:
controls_nozero <- controls_nozero[!controls_nozero$Locus %in% c("Pf3D7_13_v3-2816310-2816554-1B","Pf3D7_13_v3-2814583-2814832-2", "Pf3D7_14_v3-294506-294753-1B", "Pf3D7_13_v3-2814606-2814855-1B"),]

#calculate thresholds from single-copy controls
thresholds_controls<- controls_nozero %>% 
  group_by(SampleID) %>%
  summarize(max = max(NORM_OutputPostprocessing),
            min = min(NORM_OutputPostprocessing))

# remove outliers by IQR method
Q1 <- quantile(thresholds_controls$max, 0.25)
Q3 <- quantile(thresholds_controls$max, 0.75)

IQR <- Q3 - Q1

lower_bound <- Q1 - 1.5 * IQR
upper_bound <- Q3 + 1.5 * IQR

thresholds_controls <- thresholds_controls[thresholds_controls$max >= lower_bound & thresholds_controls$max <= upper_bound,]

#keep top controls with max > 0.01 when possible
if (any(thresholds_controls$max >= 0.01)){
  thresholds_controls <- thresholds_controls[thresholds_controls$max >= 0.01,]
}

 
  
# remove top 5 most abundant amplicons from data_norm (TO AVOID FALSE POSITIVES)
data_norm <- data_norm[!data_norm$Locus %in% ranking_amps$Locus[1:5],]

#exclude some HRP3 and PM2 ampliconS that too outlierish:
data_norm <- data_norm[!data_norm$Locus %in% c("Pf3D7_13_v3-2816310-2816554-1B","Pf3D7_13_v3-2814583-2814832-2", "Pf3D7_14_v3-294506-294753-1B", "Pf3D7_13_v3-2814606-2814855-1B"),]


#data_norm[data_norm$SampleID == "N1979742_7_S187" & !is.na(data_norm$loi) ,]


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

# RESULTS #

# Set threshold for CNV
if (any(thresholds_controls$max >= 0.01)){
  up_threshold <- min(thresholds_controls$max)
}else{
  up_threshold <- max(thresholds_controls$max)
  }

data_norm$CNV_result <- ifelse(data_norm$NORM_OutputPostprocessing > up_threshold, "CNV", "single")
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



# VISUALIZATION #

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
    #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6.5))+
    theme(axis.text.x = element_blank(),  # Remove amplicon names
          axis.title.x = element_blank())
  
  # Add the plot to the list
  plot_list[[length(plot_list) + 1]] <- p
}

# Assuming 'plot_list' contains a list of ggplot objects
combined_plot <- plot_grid(plotlist = plot_list, nrow = ceiling(sqrt(length(unique_samples))), ncol = ceiling(sqrt(length(unique_samples))))
ggsave("grid_of_plots.pdf", plot = combined_plot, width = 60, height = 50, dpi = 300, limitsize = FALSE)


# NOTES
#   one hrp2 amplicon often comes up as cnv whe it's probably not. remove it?
# esclude consistently high amplicons across samples from the run to avoid false positives? (does it affect true positives? how much?)
# TEST THIS METHOD WITH ALL CONTROLS FROM ALL RUNS
# BENCHMARK AGAINST estCNV
# ASSESS FALSE POSITIVES USING quantile99 vs max
