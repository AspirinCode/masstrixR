# #this is a text script for basic function testing
# source("R\\parseGda.R")
#
# #read HMDB file and create DB
# compoundList <- data.frame(read.table("DBs\\Txt\\lipid_test.txt", sep = "\t", header = T, stringsAsFactors = F, comment.char = " "))
# createDb(compoundList, "DBs\\SQLite\\wormjam_test_pos.sqlite", c("M+H", "M+Na"))
#
# #read peaklist
# #peakList <- data.frame(read.table("data\\examplePeaks.txt", sep = "\t", header = T, stringsAsFactors = F, comment.char = " "))
# peakList <- readGdaFile("data\\CelegansTandemPos.gda")
#
# peakList <- as.data.frame(read.table("clipboard", header = T, sep = "\t"))
#
# #annotate
# annotationResults <- mzSearch(peakList, "data\\HMDB_posH_20180609.sqlite", mode = "inMemory", mzTol = 0.005, tolType = "abs")
#
#
# #test lookup service
# annotationResults <- mzLookUp(peakList, compoundList, c("M+H", "M+Na", "M+NH4"))
#
#
# test <- aggregate(. ~ Row, data = annotationResults, FUN = paste, collapse = "_")
#
# test$RT <- sapply(strsplit(test$RT, "_"), function(i)
#   paste(unique(i), collapse = "_"))
#
# test$m.z <- sapply(strsplit(test$m.z, "_"), function(i)
#   paste(unique(i), collapse = "_"))
#
# test$adductType <- sapply(strsplit(test$adductType, "_"), function(i)
#   paste(unique(i), collapse = "_"))
#
# test$adductMass <- sapply(strsplit(test$adductMass, "_"), function(i)
#   paste(unique(i), collapse = "_"))
#
# test$neutralFormula <- sapply(strsplit(test$neutralFormula, "_"), function(i)
#   paste(unique(i), collapse = "_"))
#
# test$neutralMass <- sapply(strsplit(test$neutralMass, "_"), function(i)
#   paste(unique(i), collapse = "_"))
#
#
# #write.table(test, "E://clipboard.txt", row.names = F, quote = F, sep = "\t")
