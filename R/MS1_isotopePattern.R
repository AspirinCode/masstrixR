#' Prediction of isotope pattern
#'
#' This function predicts the isotope pattern for a given ion formula and adduct definition.
#'
#' @param ionFormula The formula of the ion for which the isotopic pattern shall be calculated, e.g. C6H12O6Na for the M+Na adduct of glucose
#' @param charge Charge of the ion.
#' @param plotit Boolean value if the isotope pattern shall be plotted
#' @param treshold Intensity treshold for lowest isotope
#' @param resolution resolution of the MS
#' @return Returns a MSnbase Spectrum 1 object containing the calculated isotope pattern
#' @examples
#' predictIsoPattern()
#' @export
predictIsoPattern <- function(ionFormula, charge, plotit = FALSE, treshold = 0.01, resolution = 50000) {

  require(MSnbase)

  # get adduct calculation list
  # adductCalc <- getAdductCalc()
  # charge <- as.numeric(adductCalc[[adduct]][5])

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

  # return values
  return(isotopeSpectrum)
}
