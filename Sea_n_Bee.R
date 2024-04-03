library(ggplot2)
library(ggpubr)
library(dplyr)
library(cowplot)
library(tools)

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
# Define the path to the directory containing all files
dir_path <- "all_amplicon_coverage_files"

# Get a list of all files in the directory
files <- list.files(path = dir_path, pattern = "\\.txt$", full.names = TRUE)

process_file <- function(coverage_file) {
  
  cat("\n", "Processing", basename(coverage_file), "\n")
  
  loci_of_interest <- readRDS("loci_of_interest.RDS")
  
  data <- read.table(coverage_file, header = T)
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
  
  if (nrow(data[!grepl("(?i)Dd2|PM|HB3", data_filtered$SampleID) & grepl("(?i)3D7", data_filtered$SampleID), ]) == 0) {
    
    print(paste0("No single-copy controls for ", basename(coverage_file), ". Unable to look for CNV"))
    return()  # Exit the function when there are no single-copy controls
    
  }
  
  # calculate proportions of amplicons for each sample
  data_norm <- data_filtered %>% 
    group_by(SampleID) %>%
    summarize(NORM_OutputPostprocessing = OutputPostprocessing / sum(OutputPostprocessing))
  
  data_norm <- as.data.frame(cbind(data_norm, Locus = data_filtered$Locus))
  
  unique_samples <- sort(unique(data_norm$SampleID))
  
  
  ##### PCA ######
  
  library(tidyr)
  library(ggrepel)
  
  # Pivot the data to have SampleID as rows and Locus as columns
  data_wide <- pivot_wider(data_norm, names_from = Locus, values_from = NORM_OutputPostprocessing)
  
  data_for_pca <- data_wide[, -1]
  
  # Remove constant or zero variance columns
  constant_cols <- sapply(data_for_pca, function(x) length(unique(x)) == 1)
  zero_variance_cols <- apply(data_for_pca, 2, function(x) var(x) == 0)
  cols_to_remove <- constant_cols | zero_variance_cols
  data_for_pca_filtered <- data_for_pca[, !cols_to_remove]
  
  # remove loci of interest from this analysis
  loci_names <- unlist(loci_of_interest)
  cols_to_keep <- setdiff(names(data_for_pca_filtered), loci_names)
  data_for_pca_filtered <- data_for_pca_filtered[, cols_to_keep]
  
  
  # #boxplot of amplicon proportions
  data_for_pca_filtered_long <- gather(data_for_pca_filtered, key = "Locus", value = "Value")
  
  # Calculate the sd of each Locus
  sd_ranking <- data_for_pca_filtered_long %>%
    group_by(Locus) %>%
    summarize(median_value = sd(Value, na.rm = TRUE))
  
  # Sort the levels of Locus based on the median values
  data_for_pca_filtered_long$Locus <- factor(data_for_pca_filtered_long$Locus, levels = sd_ranking$Locus[order(sd_ranking$median_value)])
  
  # # Plot boxplots with sorted Locus levels
  # ggplot(data_for_pca_filtered_long, aes(x = Locus, y = log(Value))) +
  #   geom_boxplot() +
  #   labs(title = "Boxplots of Locus Values",
  #        x = "Locus",
  #        y = "Value") +
  #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  # Rotate x-axis labels for better readability
  # 
  
  
  # extract n amplicons with lowest sd:
  low_sd_amps <- levels(data_for_pca_filtered_long$Locus)[1:125]
  
  data_for_pca_filtered <- data_for_pca_filtered[low_sd_amps]
  
  
  # Perform PCA
  pca_result <- prcomp(data_for_pca_filtered, scale. = TRUE)  # Remove SampleID column before PCA, then scale data if necessary
  
  # Extract PC scores
  pc_scores <- pca_result$x
  
  # Convert PC scores to dataframe and add SampleID column
  pc_scores_df <- as.data.frame(pc_scores)
  pc_scores_df$SampleID <- data_wide$SampleID
  
  # Calculate Mahalanobis distance
  mah_dist <- mahalanobis(pc_scores_df[, c("PC1", "PC2")], colMeans(pc_scores_df[, c("PC1", "PC2")]), cov(pc_scores_df[, c("PC1", "PC2")]))
  alpha <- 0.05
  mah_threshold <- qchisq(1 - alpha, df = 2)  # Chi-square threshold
  
  # Identify outliers
  pc_scores_df <- pc_scores_df %>%
    mutate(is_outlier = ifelse(mah_dist > mah_threshold, TRUE, FALSE))
  
  #isolate outliers for later
  outlier_samples <- pc_scores_df[pc_scores_df$is_outlier,]$SampleID
  
  # Calculate percentage variance explained by each PC
  variance_explained <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 2)
  
  # Plot PCA results with outliers colored in red and percentage variance explained
  pca_outlier <- ggplot(pc_scores_df, aes(x = PC1, y = PC2, color = is_outlier, label = SampleID)) +
    geom_point() +
    geom_text_repel(data = subset(pc_scores_df, is_outlier), aes(label = SampleID), size = 3, color = "red", segment.color = "red", segment.size = 0.5) + 
    scale_color_manual(values = c("black", "red")) + 
    labs(title = "PCA Plot of Amplicon Proportions",
         x = paste0("PC1 (", variance_explained[1], "%)"),
         y = paste0("PC2 (", variance_explained[2], "%)")) +
    theme_minimal()
  
  
  
  
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
    
    # keep 30 least abundant amplicons (most informative)
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
  data_norm$outliers<- ifelse(data_norm$SampleID %in% outlier_samples, "bad", "good")
  
  #check results: good QC samples that have CNV for loci of interest
  results <- data_norm[!is.na(data_norm$loi) & 
                         data_norm$CNV_result == "CNV" & 
                         data_norm$QC == "good" &
                         data_norm$slope == "good" &
                         data_norm$outliers == "good",]
  
  #remove single-copy controls (since i used the mean of the max value of single-copy controls, it's expected that some end up as false positives. it's ok)
  results <- results[!(grepl("(?i)3D7", results$SampleID) & !grepl("(?i)(Dd2|PM|HB3)", results$SampleID)), ]
  results
  
  print(paste0("There are ", length(unique(results$SampleID)), " samples with good QC that have CNV for loci of interest out of ", length(unique_samples)))
  
  write.csv(results, paste0(file_path_sans_ext(basename(coverage_file)), "_CNV_results.csv"), row.names = F)
  
  
  
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
  
  #output images
  ggsave(paste0(file_path_sans_ext(basename(coverage_file)), "_histogram_slopes.pdf"), plot = histogram_slopes, width = 6, height = 4, dpi = 300, limitsize = FALSE)
  ggsave(paste0(file_path_sans_ext(basename(coverage_file)), "_pca_outliers.pdf"), plot = pca_outlier, width = 10, height = 8, dpi = 300, limitsize = FALSE)
  ggsave(paste0(file_path_sans_ext(basename(coverage_file)),"_grid_of_plots.pdf"), plot = combined_plot, width = 60, height = 50, dpi = 300, limitsize = FALSE)

  
  cat("COMPLETE!", "\n")
}


for (file in files) {
  process_file(file)
}
