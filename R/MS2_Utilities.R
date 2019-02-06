#'
#'
#' @export
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

#'
#'
#' @export
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

#'
#'
#' @export
containsFragmentDifference <- function(ms2spectrum, fragmentMassDifference, mzTol = 0.005, mzTolType = "abs") {

  # get all mz values
  mz <- mz(ms2spectrum)

  #make all combinations
  mzCombinations <- expand.grid(mz ,mz)

  # calcluate mass differences
  mzCombinations$mzDiff <- abs(mzCombinations$Var1 - mzCombinations$Var2)

  if(mzTolType == "abs") {

    containsFragmentDiff <- any(abs(mzCombinations$mzDiff - fragmentMassDifference) < mzTol)

  } else if(mzTolType == "ppm") {

    containsFragmentDiff <- any(abs((mzCombinations$mzDiff - fragmentMassDifference) / fragmentMassDifference * 1e6) < mzTol)

  } else {

    stop("unknown mzTolType")

  }

  # return result
  return(containsFragmentDiff)

}


findCommonFragments <- function() {

}
