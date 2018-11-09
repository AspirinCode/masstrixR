## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----read compound list--------------------------------------------------
library(masstrixR)

chebiExample <- system.file("extdata", "chebi_20181015.txt", package = 'masstrixR')

# read YMDB text file and create example .sqlite fil
compoundList <- data.frame(read.table(chebiExample, sep = "\t", header = T, stringsAsFactors = F, comment.char = "~"))

## ----generate sqlite-----------------------------------------------------
# list all adducts
getAdductNames()

# adducts used for DB generation
adducts <- c("M+H", "M+Na")

# create compound list for DB creation
newCompoundList <- prepareCompoundList(compoundList, adductList = adducts, rt = FALSE, ccs = FALSE, extId = FALSE)

# check if compound list is valid and create SQLite DB
if(validateCompoundList(newCompoundList, rt = FALSE, ccs = FALSE)) {
  dbFileName <- createDb(newCompoundList, "chebi_20181015_pos_MH_MNa")
}
print(dbFileName)

## ----read .gda file, message=FALSE---------------------------------------
# read example .gda file
gdaFile <- system.file("extdata", "NaAcHILICPos_Cluster.gda", package = 'masstrixR')
exampleGDA <- readGdaFile(gdaFile)

# get row annotations with m/z values
rowAnno <- exampleGDA[[3]]
rowAnno$ClusterName <- row.names(rowAnno)

#annotate
annotationResults <- mzSearch(rowAnno, dbFileName, mzTol = 0.005, mzTolType = "abs")


## ---- message=FALSE------------------------------------------------------
# perform m/z look up
annotationResultsLookup <- mzLookUp(rowAnno, newCompoundList, mzTol = 0.005, mzTolType = "abs")

