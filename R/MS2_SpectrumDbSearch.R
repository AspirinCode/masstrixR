#' Function to find MS2 spectra in a MS2 library for spectral matching. The search parameters are a given precursor m/z and potentially its type of adduct.
#' @param precursorMz Mass of precursor for search
#' @param ms2dbFileName File name or path to MS2 database
#' @param mzTol Maximum allowed tolerance in Da or ppm
#' @param mzTolType Defines the error type for m/z search, "abs" is used for absolute mass error, "ppm" for relative error
#' @param precursorType String indicating the potential precursor adducts
#'
#' @import MSnbase
#' @export
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

  # execute query
  resultSet <- DBI::dbSendQuery(mydb, query)
  DBI::dbBind(resultSet, param = param)

  # fetch result set into a dataframe and check if an annotation was found
  ids <- DBI::dbFetch(resultSet)

  # create new Spectra object to store results
  librarySearchResults <- new("Spectra")

  # iterate through all ids and get data
  for(id in ids$intID) {

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

#' Function that calculates scores for a query and a list of result spectra
#' @param querySpectrum A Spectrum2 object
#' @param queryResults A Spectra object containing the results from a DB search
#'
#' @export
createResultsSet <- function(querySpectrum, queryResults, align = TRUE, mzTol = 0.005, treshold = 1,
                             plotIt = TRUE, storePlot = FALSE, prefix = "", dataPath = "", title = NA) {

  # create empty data for results
  resultSet <- data.frame()

  # loop over queryResults
  for(i in seq_along(queryResults)) {
    # perform spectral matching
    # forward matching
    forwardScore <- forwardDotProduct(querySpectrum, queryResults[[i]],
                                     align = align, mzTol = mzTol, treshold = treshold)

    # reverse matching
    reverseScore <- reverseDotProduct(querySpectrum, queryResults[[i]],
                                     align = align, mzTol = mzTol, treshold = treshold)

    # number of common peaks
    matchingPeaks <- commonPeaks(querySpectrum, queryResults[[i]],
                                 align = align, mzTol = mzTol, treshold = treshold)

    noPeaks_querySpectrum <- length(mz(querySpectrum))
    noPeaks_queryResult <- length(mz(queryResults[[i]]))

    if(is.na(title)) {
      title <- paste0(prefix, " / ",
                      queryResults[i]@elementMetadata$name,
                      ", forward: ",
                      round(forwardScore * 1000, 0),
                      "/ reverse: ",
                      round(forwardScore * 1000, 0))
    }

    # make mirror plot
    makeMirrorPlot(querySpectrum, queryResults[[i]],
                         align = align, mzTol = mzTol, treshold = treshold, title = title, plotIt = plotIt)


    # store plot if TRUE
    if(storePlot) {
      plotPath <- paste0(dataPath, "\\", prefix, ".png")
      dev.copy(png, plotName, width = 1000, height = 500)
      dev.off()
    }

    # add to result set
    resultSet <- rbind.data.frame(resultSet, cbind.data.frame(prefix = prefix,
                                                              name = queryResults[i]@elementMetadata$name,
                                                              forwardScore = forwardScore,
                                                              reverseScore = reverseScore,
                                                              matchingPeaks = matchingPeaks,
                                                              noPeaks_query = noPeaks_querySpectrum,
                                                              noPeaks_library = noPeaks_queryResult))
  }

  # return result set
  return(resultSet)
}
