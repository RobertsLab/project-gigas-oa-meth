# Load packages

require(devtools)
require(methylKit, lib.loc = "/gscratch/home/yaaminiv/R/x86_64-pc-linux-gnu-library/3.6/")
require(dplyr)
sessionInfo()

## Process methylation data

analysisFiles <- list("zr3616_1_R1_val_1_val_1_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov",
                      "zr3616_2_R1_val_1_val_1_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov",
                      "zr3616_3_R1_val_1_val_1_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov",
                      "zr3616_4_R1_val_1_val_1_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov",
                      "zr3616_5_R1_val_1_val_1_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov",
                      "zr3616_6_R1_val_1_val_1_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov",
                      "zr3616_7_R1_val_1_val_1_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov",
                      "zr3616_8_R1_val_1_val_1_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov")

sampleMetadata <- data.frame("sampleID" = c("1", "2", "3", "4", "5", "6", "7", "8"),
                             "slide-label" = c("02-T1", "04-T3", "06-T1", "08-T2", "11-T4", "UK-02", "UK-04", "UK-06"),
                             "pHTreatment" = c(rep(1, times = 4), rep(0, times = 4)))
head(sampleMetadata) #Confirm dataframe creation

processedFiles <- methylKit::methRead(analysisFiles,
                                      sample.id = list("1", "2", "3", "4",
                                                       "5", "6", "7", "8"),
                                      assembly = "Roslin",
                                      treatment = sampleMetadata$pHTreatment,
                                      pipeline = "bismarkCoverage",
                                      mincov = 2)

processedFilteredFilesCov5 <- methylKit::filterByCoverage(processedFiles,
                                                          lo.count = 5, lo.perc = NULL,
                                                          high.count = NULL, high.perc = 99.9) %>%
  methylKit::normalizeCoverage(.)

save.image("methylKit.RData") #Save R Data in case R crashes

# Characterize general methylation

## Sample-specific descriptive statistics

### Specify working directory for output

nFiles <- 8
fileName <- data.frame("nameBase" = rep("general-stats/percent-CpG-methylation", times = nFiles),
                       "nameBase2" = rep("general-stats/percent-CpG-coverage", times = nFiles),
                       "sample.ID" = c(1:8))
head(fileName)

fileName$actualFileName1 <- paste(fileName$nameBase, "-Filtered", "-5xCoverage", "-Sample", fileName$sample.ID, ".jpeg", sep = "")
fileName$actualFileName2 <- paste(fileName$nameBase2, "-Filtered", "-5xCoverage", "-Sample", fileName$sample.ID, ".jpeg", sep = "")
head(fileName)

### Create plots

for(i in 1:nFiles) {
  jpeg(filename = fileName$actualFileName1[i], height = 1000, width = 1000)
  methylKit::getMethylationStats(processedFilteredFilesCov5[[i]], plot = TRUE, both.strands = FALSE)
  dev.off()
}

for(i in 1:nFiles) {
  jpeg(filename = fileName$actualFileName2[i], height = 1000, width = 1000)
  methylKit::getCoverageStats(processedFilteredFilesCov5[[i]], plot = TRUE, both.strands = FALSE)
  dev.off()
}

## Comparative analysis

methylationInformationFilteredCov5 <- methylKit::unite(processedFilteredFilesCov5,
                                                       destrand = FALSE,
                                                       mc.cores = 2)
head(methylationInformationFilteredCov5)

save.image("methylKit.RData")

clusteringInformationFilteredCov5 <- methylKit::clusterSamples(methylationInformationFilteredCov5, dist = "correlation", method = "ward", plot = FALSE)

jpeg(filename = "general-stats/Full-Sample-Pearson-Correlation-Plot-FilteredCov5Destrand.jpeg", height = 1000, width = 1000)
methylKit::getCorrelation(methylationInformationFilteredCov5, plot = TRUE)
dev.off()

jpeg(filename = "general-stats/Full-Sample-CpG-Methylation-Clustering-FilteredCov5Destrand.jpeg", height = 1000, width = 1000)
methylKit::clusterSamples(methylationInformationFilteredCov5, dist = "correlation", method = "ward", plot = TRUE)
dev.off()

jpeg(filename = "general-stats/Full-Sample-Methylation-PCA-FilteredCov5Destrand.jpeg", height = 1000, width = 1000)
methylKit::PCASamples(methylationInformationFilteredCov5)
dev.off()

jpeg(filename = "general-stats/Full-Sample-Methylation-Screeplot-FilteredCov5Destrand.jpeg", height = 1000, width = 1000)
methylKit::PCASamples(methylationInformationFilteredCov5, screeplot = TRUE)
dev.off()

save.image("methylKit.RData")

# Differentially methylated loci

## Female samples

methylationInformationFilteredCov5Fem <- methylKit::reorganize(methylationInformationFilteredCov5,
                                                               sample.ids = c("1", "3", "4",
                                                                              "6", "7", "8"),
                                                               treatment = c(rep(1, times = 3),
                                                                             rep(0, times = 3)))

jpeg(filename = "general-stats/Full-Sample-CpG-Methylation-Clustering-FilteredCov5Destrand-Fem.jpeg", height = 1000, width = 1000)
methylKit::clusterSamples(methylationInformationFilteredCov5Fem, dist = "correlation", method = "ward", plot = TRUE)
dev.off()

jpeg(filename = "general-stats/Full-Sample-Methylation-PCA-FilteredCov5Destrand-Fem.jpeg", height = 1000, width = 1000)
methylKit::PCASamples(methylationInformationFilteredCov5Fem)
dev.off() #Turn off plotting device

jpeg(filename = "general-stats/Full-Sample-Methylation-Screeplot-FilteredCov5Destrand-Fem.jpeg", height = 1000, width = 1000)
methylKit::PCASamples(methylationInformationFilteredCov5Fem, screeplot = TRUE)
dev.off()

### Create covariate matrix

covariateMetadataFem <- data.frame("stage" = c("2", "1", "1",
                                               "1", "2", "3"))
head(covariateMetadataFem)

### Identify DML

differentialMethylationStatsTreatmentFem <- methylKit::calculateDiffMeth(methylationInformationFilteredCov5Fem, covariates = covariateMetadataFem, overdispersion = "MN", test = "Chisq")
head(differentialMethylationStatsTreatment)

diffMethStatsTreatment25Fem <- methylKit::getMethylDiff(differentialMethylationStatsTreatmentFem, difference = 25, qvalue = 0.01)
length(diffMethStatsTreatment25Fem$chr)
head(diffMethStatsTreatment25Fem)

diffMethStatsTreatment50Fem <- methylKit::getMethylDiff(differentialMethylationStatsTreatmentFem, difference = 50, qvalue = 0.01)
length(diffMethStatsTreatment50Fem$chr)
head(diffMethStatsTreatment50Fem)

save.image("methylKit.RData")


write.csv(diffMethStatsTreatment25, "DML/DML-pH-25-Cov5-Fem.csv")
write.csv(diffMethStatsTreatment50, "DML/DML-pH-50-Cov5-Fem.csv")

### Randomization trial

sampleMetadataFem <- data.frame("sampleID" = c("1", "3", "4",
                                               "6", "7", "8"),
                             "slide-label" = c("02-T1", "06-T1", "08-T2",
                                               "UK-02", "UK-04", "UK-06"),
                             "pHTreatment" = c(rep(1, times = 3),
                                               rep(0, times = 3)))
head(sampleMetadataFem)

#### 25% difference

pHRandTest25Fem <- NULL

for (i in 1:1000) {
  pHRandTreat <- sample(sampleMetadataFem$pHTreatment,
                        size = length(sampleMetadataFem$pHTreatment),
                        replace = FALSE)
  pHRandDML <- methylKit::reorganize(methylationInformationFilteredCov5Fem,
                                     sample.ids = sampleMetadataFem$sampleID,
                                     treatment = pHRandTreat) %>%
    methylKit::calculateDiffMeth(., covariates = covariateMetadataFem,
                                 overdispersion = "MN",
                                 test = "Chisq") %>%
    methylKit::getMethylDiff(., difference = 25, qvalue = 0.01) %>%
    nrow() -> pHRandTest25Fem[i]
}

pHRandTestResults25Fem <- t.test(pHRandTest25Fem, alternative = "less", mu = nrow(diffMethStatsTreatment25Fem))
pHRandTestResults25Fem$statistic
pHRandTestResults25Fem$parameter
pHRandTestResults25Fem$p.value

jpeg(filename = "rand-test/Random-Distribution-diff25-Fem.jpeg", height = 1000, width = 1000)
hist(pHRandTest25Fem, plot = TRUE, main = paste("t =", pHRandTestResults25Fem$statistic, "df = ", pHRandTestResults25Fem$parameter, "P-value =", pHRandTestResults25Fem$p.value), xlab = "Number of Female-DML (25% difference)")
abline(v = nrow(diffMethStatsTreatment25Fem))
dev.off()

save.image("methylKit.RData")

#### 50% difference

pHRandTest50Fem <- NULL

for (i in 1:1000) {
  pHRandTreat <- sample(sampleMetadataFem$pHTreatment,
                        size = length(sampleMetadataFem$pHTreatment),
                        replace = FALSE)
  pHRandDML <- methylKit::reorganize(methylationInformationFilteredCov5Fem,
                                     sample.ids = sampleMetadataFem$sampleID,
                                     treatment = pHRandTreat) %>%
    methylKit::calculateDiffMeth(., covariates = covariateMetadataFem,
                                 overdispersion = "MN",
                                 test = "Chisq") %>%
    methylKit::getMethylDiff(., difference = 50, qvalue = 0.01) %>%
    nrow() -> pHRandTest50Fem[i]
}

pHRandTestResults50Fem <- t.test(pHRandTest50Fem, alternative = "less", mu = nrow(diffMethStatsTreatment50Fem))
pHRandTestResults50Fem$statistic
pHRandTestResults50Fem$parameter
pHRandTestResults50Fem$p.value

jpeg(filename = "rand-test/Random-Distribution-diff50-Fem.jpeg", height = 1000, width = 1000)
hist(pHRandTest50Fem, plot = TRUE, main = paste("t =", pHRandTestResults50Fem$statistic, "df = ", pHRandTestResults50Fem$parameter, "P-value =", pHRandTestResults50Fem$p.value), xlab = "Number of Female-DML (50% difference)")
abline(v = nrow(diffMethStatsTreatment50Fem))
dev.off()

save.image("methylKit.RData")

## Indeterminate samples

methylationInformationFilteredCov5Ind <- methylKit::reorganize(methylationInformationFilteredCov5,
                                                               sample.ids = c("2",
                                                                              "5"),
                                                               treatment = c(1, 0))

jpeg(filename = "general-stats/Full-Sample-CpG-Methylation-Clustering-FilteredCov5Destrand-Ind.jpeg", height = 1000, width = 1000)
methylKit::clusterSamples(methylationInformationFilteredCov5Ind, dist = "correlation", method = "ward", plot = TRUE)
dev.off()

jpeg(filename = "general-stats/Full-Sample-Methylation-PCA-FilteredCov5Destrand-Ind.jpeg", height = 1000, width = 1000)
methylKit::PCASamples(methylationInformationFilteredCov5Ind)
dev.off()

jpeg(filename = "general-stats/Full-Sample-Methylation-Screeplot-FilteredCov5Destrand-Ind.jpeg", height = 1000, width = 1000)
methylKit::PCASamples(methylationInformationFilteredCov5Ind, screeplot = TRUE)
dev.off()

### Identify DML

differentialMethylationStatsTreatmentInd <- methylKit::calculateDiffMeth(methylationInformationFilteredCov5Ind)
head(differentialMethylationStatsTreatment)

diffMethStatsTreatment25Ind <- methylKit::getMethylDiff(differentialMethylationStatsTreatmentInd, difference = 25, qvalue = 0.01)
head(diffMethStatsTreatment25Ind)

diffMethStatsTreatment50Ind <- methylKit::getMethylDiff(differentialMethylationStatsTreatmentInd, difference = 50, qvalue = 0.01)
head(diffMethStatsTreatment50Ind)

save.image("methylKit.RData")

write.csv(diffMethStatsTreatment25, "DML/DML-pH-25-Cov5-Ind.csv")
write.csv(diffMethStatsTreatment50, "DML/06-methylKit/DML-pH-50-Cov5-Ind.csv")

### Randomization trial

sampleMetadataInd <- data.frame("sampleID" = c("2",
                                               "5"),
                             "slide-label" = c("04-T3",
                                               "11-T4"),
                             "pHTreatment" = c(1,
                                               0))
head(sampleMetadataInd)

#### 25% difference

pHRandTest25Ind <- NULL

for (i in 1:1000) {
  pHRandTreat <- sample(sampleMetadataInd$pHTreatment,
                        size = length(sampleMetadataInd$pHTreatment),
                        replace = FALSE)
  pHRandDML <- methylKit::reorganize(methylationInformationFilteredCov5Ind,
                                     sample.ids = sampleMetadataInd$sampleID,
                                     treatment = pHRandTreat) %>%
    methylKit::calculateDiffMeth(.) %>%
    methylKit::getMethylDiff(., difference = 25, qvalue = 0.01) %>%
    nrow() -> pHRandTest25Ind[i]
}

pHRandTestResults25Ind <- t.test(pHRandTest25Ind, alternative = "less", mu = nrow(diffMethStatsTreatment25Ind))
pHRandTestResults25Ind$statistic
pHRandTestResults25Ind$parameter
pHRandTestResults25Ind$p.value

jpeg(filename = "rand-test/Random-Distribution-diff25-Ind.jpeg", height = 1000, width = 1000)
hist(pHRandTest25Ind, plot = TRUE, main = paste("t =", pHRandTestResults25Ind$statistic, "df = ", pHRandTestResults25Ind$parameter, "P-value =", pHRandTestResults25Ind$p.value), xlab = "Number of Indeterminate-DML (25% difference)")
abline(v = nrow(diffMethStatsTreatment25Ind))
dev.off()

save.image("methylKit.RData")

#### 50% difference

for (i in 1:1000) {
  pHRandTreat <- sample(sampleMetadataInd$pHTreatment,
                        size = length(sampleMetadataInd$pHTreatment),
                        replace = FALSE)
  pHRandDML <- methylKit::reorganize(methylationInformationFilteredCov5Ind,
                                     sample.ids = sampleMetadataInd$sampleID,
                                     treatment = pHRandTreat) %>%
    methylKit::calculateDiffMeth(.) %>%
    methylKit::getMethylDiff(., difference = 50, qvalue = 0.01) %>%
    nrow() -> pHRandTest50Ind[i]
}

pHRandTestResults50Ind <- t.test(pHRandTest50Ind, alternative = "less", mu = nrow(diffMethStatsTreatment50Ind))
pHRandTestResults50Ind$statistic
pHRandTestResults50Ind$parameter
pHRandTestResults50Ind$p.value

jpeg(filename = "rand-test/Random-Distribution-diff50-Ind.jpeg", height = 1000, width = 1000)
hist(pHRandTest50Ind, plot = TRUE, main = paste("t =", pHRandTestResults50Ind$statistic, "df = ", pHRandTestResults50Ind$parameter, "P-value =", pHRandTestResults50Ind$p.value), xlab = "Number of Indeterminate-DML (50% difference)")
abline(v = nrow(diffMethStatsTreatment50Ind))
dev.off()

save.image("methylKit.RData")
