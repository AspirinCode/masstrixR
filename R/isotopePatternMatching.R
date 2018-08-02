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

  print(charge)

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

  centro <- as.data.frame(centro)
  colnames(centro) <- c("mz", "int")

  isotopeSpectrum <- new("Spectrum1",
                         mz = centro$mz,
                         intensity = centro$int,
                         centroided = TRUE)

  # return values
  return(isotopeSpectrum)
}



