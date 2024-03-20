library(ggplot2)
library(dplyr)


data <- read.table("230623_M07977_0012_000000000-KJGHP_amplicon_coverage_DD2.txt", header = T)
data <- data[,-3:-4]

#normalize read count for each sample
data_norm <- data %>% 
  group_by(SampleID) %>%
  summarize(NORM_OutputPostprocessing = OutputPostprocessing / sum(OutputPostprocessing))

data_norm <- as.data.frame(cbind(data_norm, Locus = data$Locus))


# add loci of interest
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

# Add the 'loi' column to the dataframe
data_norm$loi <- loi

#extract single-copy controls
controls <- data_norm[!grepl("(?i)Dd2|PM|HB3", data_norm$SampleID) & grepl("(?i)3D7", data_norm$SampleID), ]
controls_nozero <- controls[controls$NORM_OutputPostprocessing != 0, ]

#calculate thresholds deom single-copy controls
thresholds_controls<- controls_nozero %>% 
  group_by(SampleID) %>%
  summarize(q_up = quantile(NORM_OutputPostprocessing,0.975),
            q_down = quantile(NORM_OutputPostprocessing,0.025),
            max = max(NORM_OutputPostprocessing),
            min = min(NORM_OutputPostprocessing))



#grab a sample
unique_samples <- unique(data_norm$SampleID)

xsample <- unique_samples[92]

sample <- data_norm[data_norm$SampleID == xsample , ] #& data_norm$NORM_OutputPostprocessing != 0

# Calculate the median of OutputPostprocessing
up_threshold_q <- mean(thresholds_controls$q_up)
down_threshold_q <- mean(thresholds_controls$q_down)
up_threshold <- max(thresholds_controls$max)
down_threshold <- min(thresholds_controls$min)

# Scatterplot with sorted samples and median line
ggplot(sample, aes(x = reorder(Locus, -NORM_OutputPostprocessing), y = log(NORM_OutputPostprocessing), color = loi)) +
  geom_point(size=4, alpha =0.7) +
  geom_hline(yintercept = log(up_threshold), linetype = "dashed", color = "red") +
  #geom_hline(yintercept = log(down_threshold_q), linetype = "dashed", color = "blue") +
  labs(x = "Locus", y = "Read Count", title = xsample) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6.5))


