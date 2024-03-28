# Sea'n'Bee 🌊🐝

![logo](https://github.com/manuelgug/CNV_detection/blob/main/seanbee_Logo_.png)

__Sea'n'Bee__ identifies and flags potential CNV (Copy Number Variation) amplicons from the [mad4hatter](https://github.com/EPPIcenter/mad4hatter) panel for subsequent qPCR confirmation.

## Inputs
1. Output amplicon coverage file from mad4hatter (*v0.1.8* was used, although coverage files should be pretty similar across versions).
2. Loci of interest provided in an RDS file.

## Outputs
1. CSV file containing amplicons of loci of interest flagged as CNV for each sample.
2. Histogram of slopes with 2 sigma upper and lower thresholds.
3. Plots of amplicon proportions for each sample.

## Steps

### QC1: Sample filtering
- Remove samples with a median read count below 100.

### Normalization
- Calculate proportions of amplicons from read counts.

### QC2: Slope analysis
- Sort amplicons by proportion and subset the 30 least abundant ones.
- Fit a linear regression model and calculate slope.
- Perform Shapiro-Wilk test on slope values to assess run quality.
- Remove samples with abnormal slopes (> and < 2 sigma from the mean slope value).

### QC3: Outlier amplicons removal
- Remove 5 most consistently abundant amplicons from the data.
- Exclude outlierish HRP3 and PM2 amplicons.

### CNV cutoff calculation from single-copy controls
- exclude amplicons with no reads
- calculate max and min amplicon proportions for each control
- when there are many single-copy controls, remove outlier controls based on max amplicon proportions using IQR method
- if possible, keep only controls with max proportion > 0.01
- if any control has max proportion >= 0.01, set lowest max proportion as CNV threshold, else use highest max proportion

### CNV Flagging
- Flag loci of interest that passed all QC steps with proportion above cutoff as CNV.

## Example Outputs

![logo](https://github.com/manuelgug/CNV_detection/blob/main/dd2_gradient.png)
Plots depicting CNV detection across a gradient of 3D7-DD2 mixes. Red line is the CNV cutoff. Correct detection of CNV for MDR1 was possible down to 25% of DD2 prevalence in the mix. Quantitative results can be seen in the *example_output.csv* file.
