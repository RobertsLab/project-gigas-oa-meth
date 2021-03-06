# Load packages


require(BiocGenerics, lib.loc = "/gscratch/home/yaaminiv/R/x86_64-pc-linux-gnu-library/3.6/")
require(S4Vectors, lib.loc = "/gscratch/home/yaaminiv/R/x86_64-pc-linux-gnu-library/3.6/")
require(IRanges, lib.loc = "/gscratch/home/yaaminiv/R/x86_64-pc-linux-gnu-library/3.6/")
require(GenomeInfoDb, lib.loc = "/gscratch/home/yaaminiv/R/x86_64-pc-linux-gnu-library/3.6/")
require(GenomicRanges, lib.loc = "/gscratch/home/yaaminiv/R/x86_64-pc-linux-gnu-library/3.6/")
require(methylKit, lib.loc = "/gscratch/home/yaaminiv/R/x86_64-pc-linux-gnu-library/3.6/")

require(dplyr)

sessionInfo()

# Load data

setwd("/gscratch/scrubbed/yaaminiv/Manchester/analyses/methylKit") #Set working directory
load("methylKit.RData") #Load R data from DML identification

# Female samples

sampleMetadataFem <- data.frame("sampleID" = c("1", "3", "4",
                                               "6", "7", "8"),
                             "slide-label" = c("02-T1", "06-T1", "08-T2",
                                               "UK-02", "UK-04", "UK-06"),
                             "pHTreatment" = c(rep(1, times = 3),
                                               rep(0, times = 3))) #Create metadata table
sampleMetadataFem$sampleID <- as.character(sampleMetadataFem$sampleID) #Convert to character format
head(sampleMetadataFem)

pHRandTest75Fem <- NULL

for (i in 1:1000) {
  pHRandTreat <- sample(sampleMetadataFem$pHTreatment,
                        size = length(sampleMetadataFem$pHTreatment),
                        replace = FALSE)
  pHRandDML <- methylKit::reorganize(methylationInformationFilteredCov5Fem,
                                     sample.ids = sampleMetadataFem$sampleID,
                                     treatment = pHRandTreat) %>%
    methylKit::calculateDiffMeth(., covariates = covariateMetadataFem,
                                 overdispersion = "MN",
                                 test = "Chisq", mc.cores = 40) %>%
    methylKit::getMethylDiff(., difference = 75, qvalue = 0.01) %>%
    nrow() -> pHRandTest75Fem[i]
}

save.image("methylKit.RData")

pHRandTestResults75Fem <- t.test(pHRandTest75Fem, alternative = "less", mu = nrow(diffMethStatsTreatment75Fem))
pHRandTestResults75Fem$statistic
pHRandTestResults75Fem$parameter
pHRandTestResults75Fem$p.value

jpeg(filename = "rand-test/Random-Distribution-diff75-Fem.jpeg", height = 1000, width = 1000)
hist(pHRandTest75Fem, plot = TRUE, main = paste("t =", pHRandTestResults75Fem$statistic, "df = ", pHRandTestResults75Fem$parameter, "P-value =", pHRandTestResults75Fem$p.value), xlab = "Number of Female-DML (75% difference)")
abline(v = nrow(diffMethStatsTreatment75Fem))
dev.off()

save.image("methylKit.RData")

# Indeterminate samples

sampleMetadataInd <- data.frame("sampleID" = c("2",
                                               "5"),
                             "slide-label" = c("04-T3",
                                               "11-T4"),
                             "pHTreatment" = c(1,
                                               0))
sampleMetadataInd$sampleID <- as.character(sampleMetadataInd$sampleID) #Convert to character format
head(sampleMetadataInd)

pHRandTest100Ind <- NULL

for (i in 1:1000) {
  pHRandTreat <- sample(sampleMetadataInd$pHTreatment,
                        size = length(sampleMetadataInd$pHTreatment),
                        replace = FALSE)
  pHRandDML <- methylKit::reorganize(methylationInformationFilteredCov5Ind,
                                     sample.ids = sampleMetadataInd$sampleID,
                                     treatment = pHRandTreat) %>%
    methylKit::calculateDiffMeth(.) %>%
    methylKit::getMethylDiff(., difference = 99, qvalue = 0.01) %>%
    nrow() -> pHRandTest100Ind[i]
}

save.image("methylKit.RData")

pHRandTestResults1000Ind <- t.test(pHRandTest100Ind, alternative = "less", mu = nrow(diffMethStatsTreatment100Ind))
pHRandTestResults1000Ind$statistic
pHRandTestResults1000Ind$parameter
pHRandTestResults1000Ind$p.value

jpeg(filename = "rand-test/Random-Distribution-diff100-Ind.jpeg", height = 1000, width = 1000)
hist(pHRandTest100Ind, plot = TRUE, main = paste("t =", pHRandTestResults1000Ind$statistic, "df = ", pHRandTestResults1000Ind$parameter, "P-value =", pHRandTestResults1000Ind$p.value), xlab = "Number of Indeterminate-DML (99% difference)")
abline(v = nrow(diffMethStatsTreatment100Ind))
dev.off()

save.image("methylKit.RData")
