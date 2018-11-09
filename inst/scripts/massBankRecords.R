################################################################################
# read massbank record
################################################################################
test <- "E:\\Project Data\\Bio-Chemoinformatics\\R\\Packages\\masstrixR\\inst\\extdata\\MassBank"

massBankFileList <- list.files(test, full.names = TRUE)

librarySpectra <- new("Spectra")

for(massBankFile in massBankFileList) {
  print(massBankFile)

  mbRecord <- readMassBankFile(massBankFile)

  librarySpectra <- append(librarySpectra, mbRecord)

}

################################################################################
# upload to DB
################################################################################
ms2dbFileName <- createMs2Db(librarySpectra, "ms2test")

################################################################################
# test DB
################################################################################
mydb <- DBI::dbConnect(RSQLite::SQLite(), ms2dbFileName)


DBI::dbListTables(mydb)

DBI::dbGetQuery(mydb, 'SELECT * FROM spectra WHERE id = "RP000101"')
DBI::dbGetQuery(mydb, 'SELECT * FROM metaData WHERE id = "RP000101"')
DBI::dbGetQuery(mydb, 'SELECT * FROM massSpec WHERE id = "RP000101"')

#disconnect DB
DBI::dbDisconnect(mydb)

################################################################################
# example Precursor search
################################################################################
#read example spectrum

trypto <- readMgfData("E:\\Project Data\\Bio-Chemoinformatics\\R\\Packages\\masstrixR\\inst\\extdata\\trypto.mgf")
ms2 <- spectra(filterMsLevel(trypto, msLevel = 2))

results <- searchByPrecursor(ms2[[1]]@precursorMz, ms2dbFileName, mzTol = 0.005, mzTolType = "abs")

for(i in 1:length(results)) {
  plot(ms2[[1]], results[[i]])
  print(compareSpectra(ms2[[1]], results[[i]], fun = "dotproduct"))
}

