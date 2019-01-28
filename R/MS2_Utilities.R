

containsProductIon <- function(ms2spectrum, productIonMz, mzTol = 0.005, mzTolType = "abs") {

  # get mz values to search in
  mz <- mz(ms2spectrum)

  containsPi <- FALSE

  if(mzTolType == "abs") {

    containsPi <- any(abs(mz - productIonMz) < mzTol)

  } else if(mzTolType == "ppm") {

    containsPi <- any(abs((mz - productIonMz) / productIonMz * 1e6) < mzTol)

  } else {

    stop("unknown mzTolType")

  }

  #return result
  return(containsPi)

}

containsNeutralLossIon <- function(ms2spectrum, neutralLossMass, mzTol = 0.005, mzTolType = "abs") {

  productIonMz <- precursorMz(ms2spectrum) - neutralLossMass


  # get mz values to search in
  mz <- mz(ms2spectrum)

  containsNl <- FALSE

  if(mzTolType == "abs") {

    containsNl <- any(abs(mz - productIonMz) < mzTol)

  } else if(mzTolType == "ppm") {

    containsNl <- any(abs((mz - productIonMz) / productIonMz * 1e6) < mzTol)

  } else {

    stop("unknown mzTolType")

  }

  #return result
  return(containsNl)

}

containsFragmentDifference <- function() {

}

findCommonFragments <- function() {

}
