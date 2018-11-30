#' This function is used to match two masses from the positive and negative ionization mode to check it they are potentially derived from the same metabolite. Based on the given adducts, all combinations are tested from positive to negative ion mode and the other way. If for both cases the mass error is below the given mass error it is assumed that the masses are derived from the same metabolite.
#'
#' @param posMz m/z value of positive ionization mode
#' @param negMz m/z value of negative ionization mode
#' @param posAdducts vector with adducts that are used for the positive ionization mode
#' @param negAdducts vector with adducts that are used for the negative ionization mode
#' @param mzTol m/z error, numeric value
#' @param mzTolType type of error used, absolute (Da) or relative (ppm)
matchMz <- function(posMz, negMz, posAdducts, negAdducts, mzTol = 0.005, mzTolType = "abs") {

  # sanity checks
  # check if adduct definitions are good
  if(!all(posAdducts %in% getAllPosModeAdducts())) {
    stop("one or more adducts in the posAdducts in not fitting")
  }

  if(!all(negAdducts %in% getAllNegModeAdducts())) {
    stop("one or more adducts in the negAdducts in not fitting")
  }

  # get all adduct calculation rules
  adductCalc <- getAdductCalc()

  # make combinations
  adductCombinations <- expand.grid(posAdducts, negAdducts)

  # make df
  df <- data.frame(posAdduct = as.character(adductCombinations$Var1),
                   negAdduct = as.character(adductCombinations$Var2),
                   posMz = posMz,
                   negMz = negMz, stringsAsFactors = FALSE)

  # iterate and calculate all combinatoins
  for(i in 1:nrow(df)) {

    posAdduct <- df$posAdduct[i]
    negAdduct <- df$negAdduct[i]

    # from pos to neg
    df$theoNeutralFromPos[i] <- (df$posMz[i] - as.numeric(adductCalc[[posAdduct]][2])) / as.numeric(adductCalc[[posAdduct]][1])
    df$theoNegFromPos[i] <- df$theoNeutralFromPos[i] * as.numeric(adductCalc[[negAdduct]][1]) + as.numeric(adductCalc[[negAdduct]][2])

    # from neg to pos
    df$theoNeutralFromNeg[i] <- (df$negMz[i] - as.numeric(adductCalc[[negAdduct]][2])) / as.numeric(adductCalc[[negAdduct]][1])
    df$theoPosFromNeg[i] <- df$theoNeutralFromNeg[i] * as.numeric(adductCalc[[posAdduct]][1]) + as.numeric(adductCalc[[posAdduct]][2])
  }

  # select fitting adducts
  if(mzTolType == "abs") {
    filteredDf <- df[which(abs(df$negMz - df$theoNegFromPos) < mzTol & abs(df$posMz - df$theoPosFromNeg) < mzTol),]
  } else if(mzTolType == "ppm") {
    filteredDf <- NULL
  } else {
    stop("unknown mzTolType")
  }

  # return resulting DF
  return(filteredDf)
}

matchMassDiff <- function(mz1, mz2, mzDiff, mzTol = 0.005, mzTolType = "abs") {

  # check which mz is largerd
  if(mz1 > mz2) {

    if(abs(mz2 + mzDiff - mz1) < mzTol & abs(mz1 - mzDiff - mz2) < mzTol) {
      return(TRUE)
    } else {
      return(FALSE)
    }

  } else {

    if(abs(mz1 + mzDiff - mz2) < mzTol & abs(mz2 - mzDiff - mz1) < mzTol) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
}
