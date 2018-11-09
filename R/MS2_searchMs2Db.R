searchByPrecursor <- function(precursorMz, ms2dbFileName, mzTol = 0.005, mzTolType = "abs", precursorType = NA) {

  # some sanity checks
  if(!is.na(precursorType) & !any(precursorType %in% getAdductNames())) {
    stop("invalid precursor type")
  }

  ######################################################
  # m/z boundaries
  ######################################################
  # calc lower and upper mz boundary
  if (mzTolType == "abs") {
    lowerMz <- precursorMz - mzTol
    upperMz <- precursorMz + mzTol
  } else if (mzTolType == "ppm") {
    lowerMz <- precursorMz - mzTol / 1e6 * precursorMz
    upperMz <- precursorMz + mzTol / 1e6 * precursorMz
  } else {
    stop("Unknown mzTolType defined!")
  }

  ################################################################################
  # connect to DB and generate query
  ################################################################################
  mydb <- DBI::dbConnect(RSQLite::SQLite(), ms2dbFileName)

  # generate query
  query <- "SELECT intID FROM massSpec WHERE (precursorMz BETWEEN :lowerMz AND :upperMz)"
  param <- list(
    lowerMz = lowerMz,
    upperMz = upperMz
  )

  # add part for precursorType
  if (!is.na(precursorType)) {
    query <- paste0(query, "AND (precursorType = :precursorType)")
    param <- c(param, list(
      precursorType = precursorType
    ))
  }

  print(query)

  print(param)

  # execute query
  resultSet <- DBI::dbSendQuery(mydb, query)
  DBI::dbBind(resultSet, param = param)

  # fetch result set into a dataframe and check if an annotation was found
  ids <- DBI::dbFetch(resultSet)

  # create new Spectra object to store results
  librarySearchResults <- new("Spectra")

  # iterate through all ids and get data
  for(id in ids$intID) {

    print(id)

    # get spectrum metadata
    metadataRs <- DBI::dbSendQuery(mydb, "SELECT * FROM metadata WHERE intID = :x")
    DBI::dbBind(metadataRs, param = list(x = id))
    metaData <- DBI::dbFetch(metadataRs)

    DBI::dbClearResult(metadataRs)

    # get spectrum
    spectrumRs <- DBI::dbSendQuery(mydb, "SELECT * FROM spectra WHERE intID = :x")
    DBI::dbBind(spectrumRs, param = list(x = id))
    spectrum <- DBI::dbFetch(spectrumRs)

    DBI::dbClearResult(spectrumRs)

    # get mass spec details
    massSpecRs <- DBI::dbSendQuery(mydb, "SELECT * FROM massSpec WHERE intID = :x")
    DBI::dbBind(massSpecRs, param = list(x = id))
    massSpec <- DBI::dbFetch(massSpecRs)

    # clear results from search
    DBI::dbClearResult(massSpecRs)

    # TODO: fix this to return values of above
    # make Spectrum2 oboject from data
    ms2spec <- new("Spectrum2",
                   merged = 0,
                   precScanNum = as.integer(1),
                   precursorMz = as.numeric(massSpec$precursorMz),
                   precursorIntensity = 100,
                   precursorCharge = as.integer(1),
                   mz = as.numeric(spectrum$mz),
                   intensity = as.numeric(spectrum$int),
                   collisionEnergy = as.numeric(massSpec$collisionEnergy),
                   centroided = TRUE)

    # make new Spectra object
    librarySearchSpectrum <- Spectra(ms2spec)

    print(metaData$id)

    # add annotations
    # from metaData table
    mcols(librarySearchSpectrum)$id <- metaData$id
    mcols(librarySearchSpectrum)$name <- metaData$name
    mcols(librarySearchSpectrum)$formula <- metaData$formula
    mcols(librarySearchSpectrum)$precursorType <- metaData$precursorType
    mcols(librarySearchSpectrum)$exactMass <- metaData$exactMass
    mcols(librarySearchSpectrum)$smiles <- metaData$smiles
    mcols(librarySearchSpectrum)$inchi <- metaData$inchi

    # from massSpec table
    mcols(librarySearchSpectrum)$instrument <- massSpec$instrument
    mcols(librarySearchSpectrum)$instrumentType <- massSpec$instrumentType
    mcols(librarySearchSpectrum)$msType <- massSpec$msType
    mcols(librarySearchSpectrum)$ionMode <- massSpec$ionMode
    mcols(librarySearchSpectrum)$precursorMz <- massSpec$precursorMz
    mcols(librarySearchSpectrum)$splash <- massSpec$splash
    mcols(librarySearchSpectrum)$numPeak <- massSpec$numPeak

    # append results to list
    librarySearchResults <- append(librarySearchResults, librarySearchSpectrum)

  }

  # disconnect DB
  DBI::dbDisconnect(mydb)

  # return search results as Spectra object
  return(librarySearchResults)
}
