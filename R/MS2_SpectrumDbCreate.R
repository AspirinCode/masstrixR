createMs2Db <- function(librarySpectra, ms2dbFileName) {

  ##############################################################################
  # get/generate required values
  ##############################################################################
  # generate filename
  dbFileName <- paste0(ms2dbFileName, ".sqlite")

  ##############################################################################
  # data frames for upload
  ##############################################################################
  # data frames for different information
  metaData <- data.frame()
  spectra <- data.frame()
  massSpec <- data.frame()

  # iterate over library spectra and create data frames for DB upload
  for(i in 1:length(librarySpectra)){

    # create internal ID
    intID <- stringr::str_pad(i, 6, pad = "0")

    # fill data frame
    # TODO add ion formulas ??
    metaData <- rbind.data.frame(metaData, cbind(intID = intID,
                                                 id = librarySpectra@elementMetadata$id[i],
                                                 name = librarySpectra@elementMetadata$name[i],
                                                 formula = standardizeChemFormula(librarySpectra@elementMetadata$formula[i]),
                                                 precursorType = librarySpectra@elementMetadata$precursorType[i],
                                                 exactMass = librarySpectra@elementMetadata$exactMass[i],
                                                 inchi = librarySpectra@elementMetadata$inchi[i],
                                                 inchiKey = "",
                                                 smiles = librarySpectra@elementMetadata$smiles[i]))

    spectra <- rbind.data.frame(spectra, cbind(intID = intID,
                                               id = librarySpectra@elementMetadata$id[i],
                                               mz = mz(librarySpectra[[i]]),
                                               int = intensity(librarySpectra[[i]])))

    massSpec <- rbind.data.frame(massSpec, cbind(intID = intID,
                                                 id = librarySpectra@elementMetadata$id[i],
                                                 instrument = librarySpectra@elementMetadata$instrument[i],
                                                 instrumentType = librarySpectra@elementMetadata$instrumentType[i],
                                                 msType = librarySpectra@elementMetadata$msType[i],
                                                 ionMode = librarySpectra@elementMetadata$ionMode[i],
                                                 precursorMz = precursorMz(librarySpectra[[i]]),
                                                 precursorType = librarySpectra@elementMetadata$precursorType[i],
                                                 collisionEnergy = collisionEnergy(librarySpectra[[i]]),
                                                 peaksCount = peaksCount(librarySpectra[[i]])))

  }


  ##############################################################################
  # Generate SQLite DB
  ##############################################################################
  # upload to DB
  mydb <- DBI::dbConnect(RSQLite::SQLite(), dbFileName)
  DBI::dbWriteTable(mydb, "metaData", metaData, overwrite = TRUE)
  DBI::dbWriteTable(mydb, "spectra", spectra, overwrite = TRUE)
  DBI::dbWriteTable(mydb, "massSpec", massSpec, overwrite = TRUE)

  # disconnect DB
  DBI::dbDisconnect(mydb)

  # return file name
  return(dbFileName)
}
