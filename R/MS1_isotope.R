reconstructIsoPattern <- function(peakList, filter = c("all", "highestTic", "merge"), plotit = FALSE) {

  filter <- match.arg(filter)

  isoPatternList <- list()

  # get intensity values from each sample and use only the pattern with highest TIC
  for(sample in colnames(peakList[, !colnames(peakList) %in% c("mz")])) {

    # get intensity for specific sample
    int <- peakList[,sample]

    # check if all peaks are present
    if(!any(is.na(int))) {

      # plot the final isotope pattern
      if(plotit) plot(peakList$mz, int, type = "h", ylim = c(0, max(int)))

      # create new iso pattern
      measured <- new("isoPatternSpectrum",
                      mz = mz,
                      intensity = int,
                      centroided = TRUE,
                      monoMz = min(mz),
                      featureId = annotatedCluster,
                      sample = sample)

      isoPatternList <- c(isoPatternList, measured)
    }
  }

  # perform filtering of list
  if(filter == "highestTic") {
    isoPatternList <- unlist(isoPatternList[unlist(lapply(isoPatternList, function(x) {return(x@tic)})) == max(unlist(lapply(isoPatternList, function(x) {return(x@tic)})))])[[1]]
  } else if(filter == "merge") {

  }

  #return list
  return(isoPatternList)
}

#' This function creates an SQLite database based on a given compound list. Minimum input is a metabolite name ($name), an exact mass ($exactMass), a formula ($formula) and an InChIKey ($inchikey)
#'
#' @param compoundList List of compounds that shall be added to DB
#' @param dbName A name for the database file
#' @param adductList Vector with adducts that shall be covered in the DB.
#' @return Returns the filename of the generated database
#' @examples
#' xxx
generateIsoPattern <- function(ionFormula, adduct, plotit = TRUE, treshold = 0.01, resolution = 50000) {

  require(MSnbase)

  # get adduct calculation list
  adductCalc <- getAdductCalc()
  charge <- as.numeric(adductCalc[[adduct]][5])

  # pre-check generated ion formula
  checked<-enviPat::check_chemform(isotopes, ionFormula)

  # create isotope pattern
  centro<-enviPat::isowrap(isotopes,
                           checked,
                           resmass = FALSE,
                           resolution=resolution,
                           nknots=4,
                           spar=0.2,
                           threshold=0.01,
                           charge=charge,
                           emass=0.00054858,
                           algo=2,
                           ppm=FALSE,
                           dmz="get", # retrieve dm from R=m/dm
                           frac=1/4,
                           env="Gaussian",
                           detect="centroid",
                           plotit=plotit)

  # create Spectrum1 object from it
  centro <- as.data.frame(centro)
  colnames(centro) <- c("mz", "int")
  isotopeSpectrum <- new("Spectrum1",
                         mz = centro$mz,
                         intensity = centro$int,
                         centroided = TRUE)


  predicted <- new("isoPatternSpectrum",
                  mz = centro$mz,
                  intensity = centro$int,
                  centroided = TRUE,
                  monoMz = 0,
                  featureId = "predicted",
                  formula = ionFormula,
                  sample = sample)


  # return values
  return(isotopeSpectrum)
}
