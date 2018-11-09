## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- basic isotope pattern prediction-----------------------------------
library(masstrixR)

# generate isotopic pattern for [M+Na]+ adduct of Glucose
# default parameters are used (treshold = 0.01, resolution = 50000)
#predictIsoPattern("C6H12O6Na", "M+Na")

# use of lower resolution
#predictIsoPattern("C6H12O6Na", "M+Na", resolution = 5000)

## ----read data and annotate----------------------------------------------
# load required libraries
# library(stringr)
# 
# # read cluster
# exampleGDAClusters <- readGdaFile("..\\example\\NaAcHILICPos_Cluster.gda")
# clusterRowAnno <- exampleGDAClusters[[3]]
# clusterRowAnno$cluster <- row.names(clusterRowAnno)
# 
# # read peaks
# exampleGDAPeaks <- readGdaFile("..\\example\\NaAcHILICPos_Peaks.gda")
# peaksRowAnno <- exampleGDAPeaks[[3]]
# peaksRowAnno$Peak <- row.names(peaksRowAnno)
# 
# # get sample names
# sampleNames <- exampleGDAClusters[[2]]$SampleNames
# 
# # get peak data to reconstruct isotope pattern
# peaksData <- exampleGDAPeaks[[1]]
# 
# # read YMDB text file and annotated example data
# compoundList <- data.frame(read.table("..\\DBs\\Txt\\ymdb_20180731.txt", sep = "\t", header = T, stringsAsFactors = F, comment.char = ""))
# annotationResultsLookup <- mzLookUp(clusterRowAnno[c("cluster","m.z", "RT")], compoundList, mzTol = 0.005, tolType = "abs", c("M+H", "M+Na"))

## ------------------------------------------------------------------------
# iterate through annotation results
# for(i in 1:nrow(annotationResultsLookup)) {
#   
#   annotatedCluster <- annotationResultsLookup$cluster[i]
#   
#   # some required variables
#   intFinal <- vector()
#   tic <- 0
# 
#   # get number of cluster
#   clusterNumber <- str_extract(annotatedCluster, "\\d+")
#   
#   # isolate peaknumbers
#   peaks <- peaksRowAnno$Peak[which(peaksRowAnno$`Cluster [C]` == clusterNumber)]
# 
#   # get m/z values
#   mz <- peaksRowAnno$m.z[which(peaksRowAnno$Peak %in% peaks)]
#   
#     # get intensity values from each sample and use only the pattern with highest TIC
#   for(sample in sampleNames) {
# 
#     int <- peaksData[peaks,sample]
#     int[is.na(int)] <- 0
#     
#     if(sum(int) > tic) {
#       intFinal <- int
#     }
#     
# 
#   }
#   
#   # plot the final isotope pattern
#   plot(mz, intFinal, type = "h", ylim = c(0, max(intFinal)))
#   
#   # create Spectrum1 object from measured data
#   measured <- new("Spectrum1",
#                 mz = mz,
#                 intensity = intFinal,
#                 centroided = TRUE)
#   
#   # predict theoretical spectrum
#   theoretical <- generateIsoPattern(annotationResultsLookup$ionFormula[i], annotationResultsLookup$adductType[i], plotit = FALSE)
#   
#   annotationResultsLookup$isoMatch[i] <- 1000 - compareSpectra(measured, theoretical, fun = "dotproduct") * 1000
#   
# }



