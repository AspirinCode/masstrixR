#this is a text script for basic function testing

#read HMDB file and create DB
compoundList <- data.frame(read.table("data\\hmdb.txt", sep = "\t", header = T, stringsAsFactors = F))
createDb(compoundList, "HMDB_posH_20180609", c("M+H"))
