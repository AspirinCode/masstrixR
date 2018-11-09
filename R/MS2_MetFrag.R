performMetFrag <-function(ms2spectrum, adduct, mzTol = 0.005, mzTolType = "abs") {


  #create list with settings
  settingsObject <- list()

  matrix(c(mz(isoPatternList[[1]]), mz(isoPatternList[[1]])), ncol = 2)


  scored.candidates <- run.metfrag(settingsObject)

}
