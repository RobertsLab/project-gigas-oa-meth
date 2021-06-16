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

## All samples

### Create covariate matrix

covariateMetadataAll <- data.frame("stage" = c("2", "0", "1", "1",
                                               "0", "1", "2", "3"))

### Identify DML

differentialMethylationStatsTreatmentAll <- methylKit::calculateDiffMeth(methylationInformationFilteredCov5, covariates = covariateMetadataAll, overdispersion = "MN", test = "Chisq")
head(differentialMethylationStatsTreatmentAll)

save.image("methylKit.RData")

diffMethStatsTreatment25All <- methylKit::getMethylDiff(differentialMethylationStatsTreatmentAll, difference = 25, qvalue = 0.01)
length(diffMethStatsTreatment25All$chr) #12826 DML
head(diffMethStatsTreatment25All)

diffMethStatsTreatment50All <- methylKit::getMethylDiff(differentialMethylationStatsTreatmentAll, difference = 50, qvalue = 0.01)
length(diffMethStatsTreatment50All$chr) #1599 DML
head(diffMethStatsTreatment50All)

diffMethStatsTreatment75All <- methylKit::getMethylDiff(differentialMethylationStatsTreatmentAll, difference = 75, qvalue = 0.01)
length(diffMethStatsTreatment75All$chr) #59 DML
head(diffMethStatsTreatment75All)

save.image("methylKit.RData")

write.csv(diffMethStatsTreatment25All, "DML/DML-pH-25-Cov5-All.csv")
write.csv(diffMethStatsTreatment50All, "DML/DML-pH-50-Cov5-All.csv")
write.csv(diffMethStatsTreatment75All, "DML/DML-pH-75-Cov5-All.csv")

save.image("methylKit.RData")

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

differentialMethylationStatsTreatmentFem2 <- methylKit::calculateDiffMeth(methylationInformationFilteredCov5Fem, covariates = covariateMetadataFem, overdispersion = "MN", test = "Chisq")
head(differentialMethylationStatsTreatmentFem2)

diffMethStatsTreatment25Fem <- methylKit::getMethylDiff(differentialMethylationStatsTreatmentFem2, difference = 25, qvalue = 0.01)
length(diffMethStatsTreatment25Fem$chr) #20830 DML
head(diffMethStatsTreatment25Fem)

diffMethStatsTreatment50Fem <- methylKit::getMethylDiff(differentialMethylationStatsTreatmentFem2, difference = 50, qvalue = 0.01)
length(diffMethStatsTreatment50Fem$chr) #3709 DML
head(diffMethStatsTreatment50Fem)

diffMethStatsTreatment75Fem <- methylKit::getMethylDiff(differentialMethylationStatsTreatmentFem2, difference = 75, qvalue = 0.01)
length(diffMethStatsTreatment75Fem$chr) #315 DML
head(diffMethStatsTreatment75Fem)

save.image("methylKit.RData")

write.csv(diffMethStatsTreatment25Fem, "DML/DML-pH-25-Cov5-Fem.csv")
write.csv(diffMethStatsTreatment50Fem, "DML/DML-pH-50-Cov5-Fem.csv")
write.csv(diffMethStatsTreatment75Fem, "DML/DML-pH-75-Cov5-Fem.csv")

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
head(differentialMethylationStatsTreatmentInd)

diffMethStatsTreatment25Ind <- methylKit::getMethylDiff(differentialMethylationStatsTreatmentInd, difference = 25, qvalue = 0.01)
length(diffMethStatsTreatment25Ind$chr) #129040
head(diffMethStatsTreatment25Ind)

diffMethStatsTreatment50Ind <- methylKit::getMethylDiff(differentialMethylationStatsTreatmentInd, difference = 50, qvalue = 0.01)
length(diffMethStatsTreatment50Ind$chr) #73414
head(diffMethStatsTreatment50Ind)

diffMethStatsTreatment75Ind <- methylKit::getMethylDiff(differentialMethylationStatsTreatmentInd, difference = 75, qvalue = 0.01)
length(diffMethStatsTreatment75Ind$chr) #21891
head(diffMethStatsTreatment75Ind)

diffMethStatsTreatment100Ind <- methylKit::getMethylDiff(differentialMethylationStatsTreatmentInd, difference = 99, qvalue = 0.01)
length(diffMethStatsTreatment100Ind$chr) #2861
head(diffMethStatsTreatment100Ind)

save.image("methylKit.RData")

write.csv(diffMethStatsTreatment25Ind, "DML/DML-pH-25-Cov5-Ind.csv")
write.csv(diffMethStatsTreatment50Ind, "DML/DML-pH-50-Cov5-Ind.csv")
write.csv(diffMethStatsTreatment75Ind, "DML/DML-pH-75-Cov5-Ind.csv")
write.csv(diffMethStatsTreatment100Ind, "DML/DML-pH-100-Cov5-Ind.csv")

save.image("methylKit2.RData")
