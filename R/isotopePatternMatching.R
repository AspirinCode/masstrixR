
generateIsoPattern <- function(ionFormula, adduct, plotit = TRUE, treshold = 0.01, resolution = 50000) {

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

  # return values
  return(centro)
}



