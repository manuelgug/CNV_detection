# Sea'n'Bee 🌊🐝

![logo](https://github.com/manuelgug/CNV_detection/blob/main/seanbee_Logo_.png)

__Sea'n'Bee__ flags potential CNV (Copy Number Variation) amplicons from the [mad4hatter](https://github.com/EPPIcenter/mad4hatter) panel for subsequent qPCR confirmation.

## Inputs
1. Output amplicon coverage file from mad4hatter (`v0.1.8` was used, although coverage files should be pretty similar across versions).
2. Loci of interest provided in an RDS file.

## Outputs
1. CSV file containing amplicons of loci of interest flagged as CNV for each sample.
2. PCA of with labels for outlier samples 
3. Histogram of slopes with 2σ upper and lower thresholds.
4. Plots of amplicon proportions for each sample.

## Sea'n'Bee Workflow

### QC1: Sample filtering
- Remove samples with a median read count < 100

### Normalization
- Calculate proportions of amplicons from read counts

### QC2: Outlier identification
- Subset n least variable amplicons across all samples
- Perform a PCA of amplicon proportions from samples
- Calculate Mahalanobis distance of PC1 and PC2
- Flag outliers given alpha = 0.05

### QC3: Slope analysis
- Sort amplicons by proportion and subset the 30 least abundant ones
- Fit a linear regression model and calculate slope
- Perform Shapiro-Wilk test on slope values to assess run quality
- Remove samples with abnormal slopes (> and < 2σ from the mean slope)

### QC4: Outlier amplicons removal
- Remove 5 most consistently abundant amplicons from the data
- Exclude outlierish HRP3 and PM2 amplicons

### CNV cutoff calculation from single-copy controls
- Exclude amplicons with no reads
- Calculate max and min amplicon proportions for each control
- When there are many single-copy controls, remove outlier controls based on max amplicon proportions using IQR method
- If possible, keep only controls with max proportion > 0.01
- If any control has max proportion >= 0.01, set lowest max proportion as CNV threshold, else use highest max proportion

### CNV Flagging
- Flag loci of interest that passed all QC steps with proportion above cutoff as CNV.

## Example Outputs

![logo](https://github.com/manuelgug/CNV_detection/blob/main/dd2_gradient_.png)

*CNV detection across a gradient of 3D7-DD2 mixes. Red line is the CNV cutoff. Correct detection of CNV for MDR1 was possible down to 25% of DD2 prevalence in the mix. Quantitative results can be seen in the *example_output.csv* file.*
