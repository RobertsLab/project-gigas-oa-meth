#!/gscratch/srlab/programs/R-3.6.2/bin/R
## Job Name
#SBATCH --job-name=manchester-roslin
## Allocation Definition
#SBATCH --account=srlab
#SBATCH --partition=srlab
## Resources
## Nodes
#SBATCH --nodes=1
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=10-00:00:00
## Memory per node
#SBATCH --mem=500G
##turn on e-mail notification
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yaaminiv@uw.edu
## Specify the working directory for this job
#SBATCH --chdir=/gscratch/scrubbed/yaaminiv/Manchester/analyses/methylKit

---
title: "Identification of DML"
author: "Yaamini Venkataraman"
date: "03/10/2021"
output: github_document
---

In this script, I'll get preliminary methylation information from oyster samples. Additionally, identify differentially methylated loci (DML) between pH treatments. I have any covariate information for sex and stage, but not tank, shell height, width, or wet weight because the some samples were mixed up during histology processing.

# Prepare R Markdown file

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# Install packages

```{r}
install.packages("devtools") #Install the devtools package
require(devtools) #Load devtools
#source("http://bioconductor.org/biocLite.R") #Source package from bioconductor
#biocLite("methylKit") #Install methylkit
#install_github("al2na/methylKit", build_vignettes = FALSE, repos = BiocManager::repositories(), dependencies = TRUE) #Install more methylKit options
require(methylKit) #Load methylkit
```

```{r}
require(dplyr) #Load dplyr
```

# Obtain session information

```{r}
sessionInfo()
```

# Download files from `gannet`

```{bash}
wget -r -l1 --no-parent --no-directories -A.merged_CpG_evidence.cov https://gannet.fish.washington.edu/spartina/project-gigas-oa-meth/output/bismark-roslin/ #Download files from gannet to the current directory (where script is saved). Exclude directory structure
```

```{bash}
mkdir ../output/06-methylKit #Make methylkit output directory
```

```{bash}
mv *.merged_CpG_evidence.cov ../output/06-methylKit #Move files from current directory to analyses folder
```

## Process methylation data

```{r}
analysisFiles <- list("zr3616_1_R1_val_1_val_1_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov",
                      "zr3616_2_R1_val_1_val_1_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov",
                      "zr3616_3_R1_val_1_val_1_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov",
                      "zr3616_4_R1_val_1_val_1_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov",
                      "zr3616_5_R1_val_1_val_1_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov",
                      "zr3616_6_R1_val_1_val_1_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov",
                      "zr3616_7_R1_val_1_val_1_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov",
                      "zr3616_8_R1_val_1_val_1_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov") #Put all .cov files into a list for analysis.
```

```{r}
sampleMetadata <- data.frame("sampleID" = c("1", "2", "3", "4",
                                            "5", "6", "7", "8"),
                             "slide-label" = c("02-T1", "04-T3", "06-T1", "08-T2",
                                               "11-T4", "UK-02", "UK-04", "UK-06"),
                             "pHTreatment" = c(rep(1, times = 4),
                                               rep(0, times = 4))) #Create dataframe with metadata. 1 = low, 0 = ambient. Slide label information from histology processing
head(sampleMetadata) #Confirm dataframe creation
```

I'll use `methRead` to create a methylation object from the coverage files, and include sample ID and treatment information.

```{r}
setwd(dir = "../output/06-methylKit/") #Set working directory inside code chunk
processedFiles <- methylKit::methRead(analysisFiles,
                                      sample.id = list("1", "2", "3", "4",
                                                       "5", "6", "7", "8"),
                                      assembly = "Roslin",
                                      treatment = sampleMetadata$pHTreatment,
                                      pipeline = "bismarkCoverage",
                                      mincov = 2) #Process files. Treatment specified based on ploidy status. Use mincov = 2 to quickly process reads.
```

```{r}
processedFilteredFilesCov5 <- methylKit::filterByCoverage(processedFiles,
                                                          lo.count = 5, lo.perc = NULL,
                                                          high.count = NULL, high.perc = 99.9) %>%
  methylKit::normalizeCoverage(.) #Filter coverage information for minimum 5x coverage, and remove PCR duplicates by excluding data in the 99.9th percentile of coverage with hi.perc = 99.9. Normalize coverage between samples to avoid over-sampling reads from one sample during statistical testing
```

```{r}
save.image("../output/06-methylKit/methylKit.RData") #Save R Data in case R crashes
#load("../output/06-methylKit/methylKit.RData") #Load R Data
```

# Characterize general methylation

## Sample-specific descriptive statistics

### Specify working directory for output

```{bash}
mkdir ../output/06-methylKit/general-stats
```

```{r}
nFiles <- 8 #Count number of samples
fileName <- data.frame("nameBase" = rep("../output/06-methylKit/general-stats/percent-CpG-methylation", times = nFiles),
                       "nameBase2" = rep("../output/06-methylKit/general-stats/percent-CpG-coverage", times = nFiles),
                       "sample.ID" = c(1:8)) #Create new dataframe for filenames.
head(fileName) #Confirm dataframe creation
```

```{r}
fileName$actualFileName1 <- paste(fileName$nameBase, "-Filtered", "-5xCoverage", "-Sample", fileName$sample.ID, ".jpeg", sep = "") #Create a new column for the full filename for filtered + 5x coverage + specific sample's percent CpG methylation plot
fileName$actualFileName2 <- paste(fileName$nameBase2, "-Filtered", "-5xCoverage", "-Sample", fileName$sample.ID, ".jpeg", sep = "") #Create a new column for the full filename for filtered + 5x coverage + specific sample's percent CpG coverage plot
head(fileName) #Confirm column creation
```

### Create plots

```{r}
for(i in 1:nFiles) { #For each data file
  jpeg(filename = fileName$actualFileName1[i], height = 1000, width = 1000) #Save file with designated name
  methylKit::getMethylationStats(processedFilteredFilesCov5[[i]], plot = TRUE, both.strands = FALSE) #Get %CpG methylation information
  dev.off() #Turn off plotting device
} #Plot and save %CpG methylation information
```

```{r}
for(i in 1:nFiles) { #For each data file
  jpeg(filename = fileName$actualFileName2[i], height = 1000, width = 1000) #Save file with designated name
  methylKit::getCoverageStats(processedFilteredFilesCov5[[i]], plot = TRUE, both.strands = FALSE) #Get CpG coverage information
  dev.off() #Turn off plotting device
} #Plot and save CpG coverage information
```

## Comparative analysis

```{r}
methylationInformationFilteredCov5 <- methylKit::unite(processedFilteredFilesCov5,
                                                       destrand = FALSE,
                                                       mc.cores = 2) #Combine all processed files into a single table. Use destrand = FALSE to not destrand. By default only bases with data in all samples will be kept
head(methylationInformationFilteredCov5) #Confirm unite
```

```{r}
save.image("../output/06-methylKit/methylKit.RData") #Save R Data in case R crashes
#load("../output/06-methylKit/methylKit.RData") #Load R Data
```

```{r}
clusteringInformationFilteredCov5 <- methylKit::clusterSamples(methylationInformationFilteredCov5, dist = "correlation", method = "ward", plot = FALSE) #Save cluster information as a new object
```

```{r}
jpeg(filename = "../output/06-methylKit/general-stats/Full-Sample-Pearson-Correlation-Plot-FilteredCov5Destrand.jpeg", height = 1000, width = 1000) #Save file with designated name
methylKit::getCorrelation(methylationInformationFilteredCov5, plot = TRUE) #Understand correlation between methylation patterns in different samples
dev.off()
```

```{r}
jpeg(filename = "../output/06-methylKit/general-stats/Full-Sample-CpG-Methylation-Clustering-FilteredCov5Destrand.jpeg", height = 1000, width = 1000) #Save file with designated name
methylKit::clusterSamples(methylationInformationFilteredCov5, dist = "correlation", method = "ward", plot = TRUE) #Cluster samples based on correlation coefficients
dev.off()
```

```{r}
jpeg(filename = "../output/06-methylKit/general-stats/Full-Sample-Methylation-PCA-FilteredCov5Destrand.jpeg", height = 1000, width = 1000) #Save file with designated name
methylKit::PCASamples(methylationInformationFilteredCov5) #Run a PCA analysis on percent methylation for all samples
dev.off() #Turn off plotting device
```

```{r}
jpeg(filename = "../output/06-methylKit/general-stats/Full-Sample-Methylation-Screeplot-FilteredCov5Destrand.jpeg", height = 1000, width = 1000) #Save file with designated name
methylKit::PCASamples(methylationInformationFilteredCov5, screeplot = TRUE) #Run the PCA analysis and plot variances against PC number in a screeplot
dev.off()
```

```{r}
save.image("../output/06-methylKit/methylKit.RData") #Save R Data in case R crashes
#load("../output/06-methylKit/methylKit.RData") #Load R Data
```

# Differentially methylated loci

I have three female samples and one indeterminate sample per treatment. Although I produced comparative plots for I will identify DML for these datasets separately.

## Female samples

I'll first separate the female samples from the full dataset, then conduct comparative analyses within this subset.

```{r}
methylationInformationFilteredCov5Fem <- methylKit::reorganize(methylationInformationFilteredCov5,
                                                               sample.ids = c("1", "3", "4",
                                                                              "6", "7", "8"),
                                                               treatment = c(rep(1, times = 3),
                                                                             rep(0, times = 3))) #Get female sample information from methylBase object created by unite by specifying sample IDs. Include treatment information as well
```

```{r}
jpeg(filename = "../output/06-methylKit/general-stats/Full-Sample-CpG-Methylation-Clustering-FilteredCov5Destrand-Fem.jpeg", height = 1000, width = 1000) #Save file with designated name
methylKit::clusterSamples(methylationInformationFilteredCov5Fem, dist = "correlation", method = "ward", plot = TRUE) #Cluster female samples based on correlation coefficients
dev.off()
```

```{r}
jpeg(filename = "../output/06-methylKit/general-stats/Full-Sample-Methylation-PCA-FilteredCov5Destrand-Fem.jpeg", height = 1000, width = 1000) #Save file with designated name
methylKit::PCASamples(methylationInformationFilteredCov5Fem) #Run a PCA analysis on percent methylation for female samples
dev.off() #Turn off plotting device
```

```{r}
jpeg(filename = "../output/06-methylKit/general-stats/Full-Sample-Methylation-Screeplot-FilteredCov5Destrand-Fem.jpeg", height = 1000, width = 1000) #Save file with designated name
methylKit::PCASamples(methylationInformationFilteredCov5Fem, screeplot = TRUE) #Run the PCA analysis and plot variances against PC number in a screeplot for female samples
dev.off()
```

### Create covariate matrix

I don't have covariate information (tank, wet weight, shell length) for all samples, since some tissues were mixed up during histology processing. However, I do have sex and stage information for each sample that will likely contribute to methylation patterns. Information is from [this spreadsheet](https://github.com/RobertsLab/project-oyster-oa/blob/master/data/Manchester/2018-02-27-Gigas-Histology-Classification.csv), and slide label information is consistent with `sampleMetadata`.

```{r}
covariateMetadataFem <- data.frame("stage" = c("2", "1", "1",
                                               "1", "2", "3")) #Create dataframe with stage information. 1/2/3 = early female maturation, with larger numbers indicating more maturity
head(covariateMetadataFem) #Check dataframe format
```

### Identify DML

```{r}
#Code that was used to test calculateDiffMeth parameters

#differentialMethylationStatsTreatment <- methylKit::calculateDiffMeth(methylationInformationFilteredCov5)
#differentialMethylationStatsTreatment2 <- methylKit::calculateDiffMeth(methylationInformationFilteredCov5, overdispersion = "MN", test = "Chisq")
#head(differentialMethylationStatsTreatment) #Look at differential methylation statistics
#head(differentialMethylationStatsTreatment2) #Look at differential methylation statistics
```

```{r}
differentialMethylationStatsTreatmentFem <- methylKit::calculateDiffMeth(methylationInformationFilteredCov5Fem, covariates = covariateMetadata, overdispersion = "MN", test = "Chisq") #Calculate differential methylation statistics and include covariate information.
head(differentialMethylationStatsTreatment) #Look at differential methylation statistics
```

```{r}
diffMethStatsTreatment25Fem <- methylKit::getMethylDiff(differentialMethylationStatsTreatmentFem, difference = 25, qvalue = 0.01) #Identify loci that are at least 25% different
length(diffMethStatsTreatment25Fem$chr) #Count the number of DML
head(diffMethStatsTreatment25Fem) #Confirm creation
```

```{r}
diffMethStatsTreatment50Fem <- methylKit::getMethylDiff(differentialMethylationStatsTreatmentFem, difference = 50, qvalue = 0.01) #Identify loci that are at least 50% different
length(diffMethStatsTreatment50Fem$chr) #Count the number of DML
head(diffMethStatsTreatment50Fem) #Confirm creation
```
```{r}
save.image("../output/06-methylKit/methylKit.RData") #Save R Data in case R crashes
#load("../output/06-methylKit/methylKit.RData") #Load R Data
```

```{r}
write.csv(diffMethStatsTreatment25, "../output/06-methylKit/DML-pH-25-Cov5-Fem.csv") #Save table as .csv
write.csv(diffMethStatsTreatment50, "../output/06-methylKit/DML-pH-50-Cov5-Fem.csv") #Save table as .csv
```

### Randomization trial

```{r}
# Code used to test initial randomization trials
# for (i in 1:1000) { #For each iteration
#   pHRandTreat <- sample(sampleMetadata$pHTreatment,
#                         size = length(sampleMetadata$pHTreatment),
#                         replace = FALSE) #Randomize treatment information
#   pHRandDML <- methylKit::reorganize(processedFiles,
#                                      sample.ids = sampleMetadata$sampleID,
#                                      treatment = pHRandTreat) %>%
#     methylKit::filterByCoverage(., lo.count = 5, lo.perc = NULL,
#                                 hi.count = NULL, hi.perc = 99.9) %>%
#     methylKit::normalizeCoverage(.) %>%
#     methylKit::unite(., destrand = FALSE) %>%
#     methylKit::calculateDiffMeth(.) %>%
#     methylKit::getMethylDiff(., difference = 50, qvalue = 0.01) %>%
#     nrow() -> pHRandTest25[i] #Reorganize, filter, normalize, and unite samples. Calculate diffMeth and obtain DML that are 50% and have a q-value of 0.01. Save the number of DML identified
# }
```

```{r}
sampleMetadataFem <- data.frame("sampleID" = c("1", "3", "4",
                                               "6", "7", "8"),
                             "slide-label" = c("02-T1", "06-T1", "08-T2",
                                               "UK-02", "UK-04", "UK-06"),
                             "pHTreatment" = c(rep(1, times = 3),
                                               rep(0, times = 3))) #Create dataframe with metadata for female samples. 1 = low, 0 = ambient. Slide label information from histology processing
head(sampleMetadataFem) #Confirm dataframe creation
```

#### 25% difference

```{r}
pHRandTest25Fem <- NULL #Create an empty matrix to store randomization trial results for 25% difference
```

```{r}
for (i in 1:1000) { #For each iteration
  pHRandTreat <- sample(sampleMetadataFem$pHTreatment,
                        size = length(sampleMetadataFem$pHTreatment),
                        replace = FALSE) #Randomize female treatment information
  pHRandDML <- methylKit::reorganize(methylationInformationFilteredCov5Fem,
                                     sample.ids = sampleMetadataFem$sampleID,
                                     treatment = pHRandTreat) %>%
    methylKit::calculateDiffMeth(., covariates = covariateMetadataFem,
                                 overdispersion = "MN",
                                 test = "Chisq") %>%
    methylKit::getMethylDiff(., difference = 25, qvalue = 0.01) %>%
    nrow() -> pHRandTest25Fem[i] #Shuffle treatments from methylBase object created by unite. Calculate diffMeth and obtain DML that are 25% and have a q-value of 0.01. Save the number of DML identified
}
```

```{r}
pHRandTestResults25Fem <- t.test(pHRandTest25Fem, alternative = "less", mu = nrow(diffMethStatsTreatment25Fem)) #Conduct a one-sample t-test to determine the probability of identifying more DML by chance
pHRandTestResults25Fem$statistic #t
pHRandTestResults25Fem$parameter #df
pHRandTestResults25Fem$p.value #P-value
```

```{r}
jpeg(filename = "../output/06-methylKit/rand-test/Random-Distribution-diff25-Fem.jpeg", height = 1000, width = 1000) #Save file with designated name
hist(pHRandTest25Fem, plot = TRUE, main = "", xlab = "Number of Female-DML (25% difference)") #Plot a histogram with randomization trial results
abline(v = nrow(diffMethStatsTreatment25Fem)) #Add a vertical line for the number of DML identified with getMethylDiff
dev.off()
```

```{r}
save.image("../output/06-methylKit/methylKit.RData") #Save R Data in case R crashes
#load("../output/06-methylKit/methylKit.RData") #Load R Data
```

#### 50% difference

```{r}
pHRandTest50Fem <- NULL #Create an empty matrix to store randomization trial results for 50% difference
```

```{r}
for (i in 1:1000) { #For each iteration
  pHRandTreat <- sample(sampleMetadataFem$pHTreatment,
                        size = length(sampleMetadataFem$pHTreatment),
                        replace = FALSE) #Randomize female treatment information
  pHRandDML <- methylKit::reorganize(methylationInformationFilteredCov5Fem,
                                     sample.ids = sampleMetadataFem$sampleID,
                                     treatment = pHRandTreat) %>%
    methylKit::calculateDiffMeth(., covariates = covariateMetadataFem,
                                 overdispersion = "MN",
                                 test = "Chisq") %>%
    methylKit::getMethylDiff(., difference = 50, qvalue = 0.01) %>%
    nrow() -> pHRandTest50Fem[i] #Shuffle treatments from methylBase object created by unite. Calculate diffMeth and obtain DML that are 50% and have a q-value of 0.01. Save the number of DML identified
}
```

```{r}
pHRandTestResults50Fem <- t.test(pHRandTest50Fem, alternative = "less", mu = nrow(diffMethStatsTreatment50Fem)) #Conduct a one-sample t-test to determine the probability of identifying more DML by chance
pHRandTestResults50Fem$statistic #t
pHRandTestResults50Fem$parameter #df
pHRandTestResults50Fem$p.value #P-value
```

```{r}
jpeg(filename = "../output/06-methylKit/rand-test/Random-Distribution-diff50-Fem.jpeg", height = 1000, width = 1000) #Save file with designated name
hist(pHRandTest50Fem, plot = TRUE, main = "", xlab = "Number of Female-DML (50% difference)") #Plot a histogram with randomization trial results
abline(v = nrow(diffMethStatsTreatment50Fem)) #Add a vertical line for the number of DML identified with getMethylDiff
dev.off()
```

```{r}
save.image("../output/06-methylKit/methylKit.RData") #Save R Data in case R crashes
#load("../output/06-methylKit/methylKit.RData") #Load R Data
```

## Indeterminate samples

```{r}
methylationInformationFilteredCov5Ind <- methylKit::reorganize(methylationInformationFilteredCov5,
                                                               sample.ids = c("2",
                                                                              "5"),
                                                               treatment = c(1, 0)) #Get indeterminate sample information from methylBase object created by unite by specifying sample IDs. Include treatment information as well
```

```{r}
jpeg(filename = "../output/06-methylKit/general-stats/Full-Sample-CpG-Methylation-Clustering-FilteredCov5Destrand-Ind.jpeg", height = 1000, width = 1000) #Save file with designated name
methylKit::clusterSamples(methylationInformationFilteredCov5Ind, dist = "correlation", method = "ward", plot = TRUE) #Cluster indeterminate samples based on correlation coefficients
dev.off()
```

```{r}
jpeg(filename = "../output/06-methylKit/general-stats/Full-Sample-Methylation-PCA-FilteredCov5Destrand-Ind.jpeg", height = 1000, width = 1000) #Save file with designated name
methylKit::PCASamples(methylationInformationFilteredCov5Ind) #Run a PCA analysis on percent methylation for indeterminate samples
dev.off() #Turn off plotting device
```

```{r}
jpeg(filename = "../output/06-methylKit/general-stats/Full-Sample-Methylation-Screeplot-FilteredCov5Destrand-Ind.jpeg", height = 1000, width = 1000) #Save file with designated name
methylKit::PCASamples(methylationInformationFilteredCov5Ind, screeplot = TRUE) #Run the PCA analysis and plot variances against PC number in a screeplot for indeterminate samples
dev.off()
```

### Identify DML

```{r}
differentialMethylationStatsTreatmentInd <- methylKit::calculateDiffMeth(methylationInformationFilteredCov5Ind) #Calculate differential methylation statistics.
head(differentialMethylationStatsTreatment) #Look at differential methylation statistics
```

```{r}
diffMethStatsTreatment25Ind <- methylKit::getMethylDiff(differentialMethylationStatsTreatmentInd, difference = 25, qvalue = 0.01) #Identify loci that are at least 25% different
length(diffMethStatsTreatment25Ind$chr) #Count the number of DML
head(diffMethStatsTreatment25Ind) #Confirm creation
```

```{r}
diffMethStatsTreatment50Ind <- methylKit::getMethylDiff(differentialMethylationStatsTreatmentInd, difference = 50, qvalue = 0.01) #Identify loci that are at least 50% different
length(diffMethStatsTreatment50Ind$chr) #Count the number of DML
head(diffMethStatsTreatment50Ind) #Confirm creation
```
```{r}
save.image("../output/06-methylKit/methylKit.RData") #Save R Data in case R crashes
#load("../output/06-methylKit/methylKit.RData") #Load R Data
```

```{r}
write.csv(diffMethStatsTreatment25, "../output/06-methylKit/DML-pH-25-Cov5-Ind.csv") #Save table as .csv
write.csv(diffMethStatsTreatment50, "../output/06-methylKit/DML-pH-50-Cov5-Ind.csv") #Save table as .csv
```

### Randomization trial

```{r}
sampleMetadataInd <- data.frame("sampleID" = c("2",
                                               "5"),
                             "slide-label" = c("04-T3",
                                               "11-T4"),
                             "pHTreatment" = c(1,
                                               0)) #Create dataframe with metadata for indeterminate samples. 1 = low, 0 = ambient. Slide label information from histology processing
head(sampleMetadataInd) #Confirm dataframe creation
```

#### 25% difference

```{r}
pHRandTest25Ind <- NULL #Create an empty matrix to store randomization trial results for 25% difference
```

```{r}
for (i in 1:1000) { #For each iteration
  pHRandTreat <- sample(sampleMetadataInd$pHTreatment,
                        size = length(sampleMetadataInd$pHTreatment),
                        replace = FALSE) #Randomize indeterminate treatment information
  pHRandDML <- methylKit::reorganize(methylationInformationFilteredCov5Ind,
                                     sample.ids = sampleMetadataInd$sampleID,
                                     treatment = pHRandTreat) %>%
    methylKit::calculateDiffMeth(.) %>%
    methylKit::getMethylDiff(., difference = 25, qvalue = 0.01) %>%
    nrow() -> pHRandTest25Ind[i] #Shuffle treatments from methylBase object created by unite. Calculate diffMeth and obtain DML that are 25% and have a q-value of 0.01. Save the number of DML identified
}
```

```{r}
pHRandTestResults25Ind <- t.test(pHRandTest25Ind, alternative = "less", mu = nrow(diffMethStatsTreatment25Ind)) #Conduct a one-sample t-test to determine the probability of identifying more DML by chance
pHRandTestResults25Ind$statistic #t
pHRandTestResults25Ind$parameter #df
pHRandTestResults25Ind$p.value #P-value
```

```{r}
jpeg(filename = "../output/06-methylKit/rand-test/Random-Distribution-diff25-Ind.jpeg", height = 1000, width = 1000) #Save file with designated name
hist(pHRandTest25Ind, plot = TRUE, main = "", xlab = "Number of Indeterminate-DML (25% difference)") #Plot a histogram with randomization trial results
abline(v = nrow(diffMethStatsTreatment25Ind)) #Add a vertical line for the number of DML identified with getMethylDiff
dev.off()
```

```{r}
save.image("../output/06-methylKit/methylKit.RData") #Save R Data in case R crashes
#load("../output/06-methylKit/methylKit.RData") #Load R Data
```

#### 50% difference

```{r}
pHRandTest50Ind <- NULL #Create an empty matrix to store randomization trial results for 50% difference
```

```{r}
for (i in 1:1000) { #For each iteration
  pHRandTreat <- sample(sampleMetadataInd$pHTreatment,
                        size = length(sampleMetadataInd$pHTreatment),
                        replace = FALSE) #Randomize female treatment information
  pHRandDML <- methylKit::reorganize(methylationInformationFilteredCov5Ind,
                                     sample.ids = sampleMetadataInd$sampleID,
                                     treatment = pHRandTreat) %>%
    methylKit::calculateDiffMeth(.) %>%
    methylKit::getMethylDiff(., difference = 50, qvalue = 0.01) %>%
    nrow() -> pHRandTest50Ind[i] #Shuffle treatments from methylBase object created by unite. Calculate diffMeth and obtain DML that are 50% and have a q-value of 0.01. Save the number of DML identified
}
```

```{r}
pHRandTestResults50Ind <- t.test(pHRandTest50Ind, alternative = "less", mu = nrow(diffMethStatsTreatment50Ind)) #Conduct a one-sample t-test to determine the probability of identifying more DML by chance
pHRandTestResults50Ind$statistic #t
pHRandTestResults50Ind$parameter #df
pHRandTestResults50Ind$p.value #P-value
```

```{r}
jpeg(filename = "../output/06-methylKit/rand-test/Random-Distribution-diff50-Ind.jpeg", height = 1000, width = 1000) #Save file with designated name
hist(pHRandTest50Ind, plot = TRUE, main = "", xlab = "Number of Indeterminate-DML (50% difference)") #Plot a histogram with randomization trial results
abline(v = nrow(diffMethStatsTreatment50Ind)) #Add a vertical line for the number of DML identified with getMethylDiff
dev.off()
```

```{r}
save.image("../output/06-methylKit/methylKit.RData") #Save R Data in case R crashes
#load("../output/06-methylKit/methylKit.RData") #Load R Data
```
