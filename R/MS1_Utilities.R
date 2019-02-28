#' Checking masses from different ion modes
#'
#' This function is used to match two masses from the positive and negative ionization mode to check it they are potentially derived from the same metabolite. Based on the given adducts, all combinations are tested from positive to negative ion mode and the other way. If for both cases the mass error is below the given mass error it is assumed that the masses are derived from the same metabolite.
#'
#' @param posMz m/z value of positive ionization mode
#' @param negMz m/z value of negative ionization mode
#' @param posAdducts vector with adducts that are used for the positive ionization mode
#' @param negAdducts vector with adducts that are used for the negative ionization mode
#' @param mzTol m/z error, numeric value
#' @param mzTolType type of error used, absolute (abs) or relative (ppm)
#'
#' @examples
#' xxx
#'
#' @author Michael Witting, \email{michael.witting@@helmholtz-muenchen.de}
#'
#' @seealso \code{\link{matchMassDiff}}
#' @export
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

  if(!is.null(filteredDf) & nrow(filteredDf) > 0) {
    matchingResult <- ""

    for(i in 1:nrow(filteredDf)) {

      if(i == 1) {
        matchingResult <- paste0(filteredDf$posAdduct[i], "<->", filteredDf$negAdduct[i])
      } else {
        matchingResult <- paste0(matchingResult, " / ", filteredDf$posAdduct[i], "<->", filteredDf$negAdduct[i])
      }
    }

    # return resulting DF
    return(matchingResult)
  } else {
    return(NA)
  }
}

#' Checking for mass diferences
#'
#' This functions is used to check of two measured masses are different by a give mass difference. This can be used to check if a neutral loss is found between two masses or if masses are connected by certain metabolic transformations, e.g. methylation, acetylation etc.
#'
#' @param mz1 First m/z value for which a mass difference shall be checked
#' @param mz2 Second m/z value for which a mass difference shall be checked
#' @param mzDiff Mass difference to check for
#' @param mzTol m/z error, numeric value
#' @param mzTolType type of error used, absolute (abs) or relative (ppm)
#'
#' @examples
#' xxx
#'
#' @author Michael Witting, \email{michael.witting@@helmholtz-muenchen.de}
#'
#' @seealso \code{\link{matchMz}}
#' @export
matchMassDiff <- function(mz1, mz2, mzDiff, mzTol = 0.005, mzTolType = "abs") {

  match <- FALSE

  # check which error type
  if(mzTolType == "abs") {

    # check which mz is largerd
    if(mz1 > mz2) {

      if(abs(mz2 + mzDiff - mz1) < mzTol & abs(mz1 - mzDiff - mz2) < mzTol) {
        match <- TRUE
      }

    } else {

      if(abs(mz1 + mzDiff - mz2) < mzTol & abs(mz2 - mzDiff - mz1) < mzTol) {
        match <- TRUE
      }
    }

  } else if(mzTolType == "ppm") {

    # check which mz is largerd
    if(mz1 > mz2) {

      if(abs(mz2 + mzDiff - mz1) < (mzTol / 1e6 * mz1) & abs(mz1 - mzDiff - mz2) < (mzTol / 1e6 * mz1)) {
        match <- TRUE
      }

    } else {

      if(abs(mz1 + mzDiff - mz2) < (mzTol / 1e6 * mz2) & abs(mz2 - mzDiff - mz1) < (mzTol / 1e6 * mz2)) {
        match <- TRUE
      }
    }

  } else {
    stop("unknown mzTolType")
  }

  # return result
  return(match)
}

#' Calculate Kendrick Mass Defect
#'
#' This function calculates the Kendrick mass defect (KMD) for a given mass.
#'
#' @param mz Mass for which the KMD shall be calculated
#'
#' @example
#' xxx
#'
#' @author Michael Witting, \email{michael.witting@@helmholtz-muenchen.de}
#'
#' @seealso \code{\link{calculateRefKendrickMassDefect}}
#' @seealso \code{\link{checkRmkd}}
#'
#' @export
calculateKendrickMassDefect <- function(mz) {

  # calculate KMD
  kendrickMass <- mz * 14 / 14.01565
  kmd <- kendrickMass %% 1

  # return Kendrick mass defect
  return(kmd)

}

#' Calculate referenced Kendrick Mass Defect
#'
#' This function calculates a referenced Kendrick Mass Defect (RKMD) for a given mass and a given reference KMD.
#'
#' @param mz Mass for which the RKMD shall be calculated
#' @param refkmd Reference Kendrick Mass Defect
#'
#' @examples
#' xxx
#'
#' @author Michael Witting, \email{michael.witting@@helmholtz-muenchen.de}
#'
#' @seealso \code{\link{calculateKendrickMassDefect}}
#' @seealso \code{\link{checkRmkd}}
#'
#' @export
calculateRefKendrickMassDefect <- function(mz, refkmd) {

  # calculate KMD
  kmd <- calculateKendrickMassDefect(mz)

  # calculate RKMD
  rmkd <- (kmd - refkmd) / 0.013399

  # return referenced Kendrick mass defect
  return(rmkd)

}

#' Check if RMKD is negative integer
#'
#' This function checks if the RKMD is a negative integer
#'
#' @param rmkd calculated RMKD
#' @param error allowed error margin for RMKD
#'
#' @author Michael Witting, \email{michael.witting@@helmholtz-muenchen.de}
#'
#' @seealso \code{\link{calculateKendrickMassDefect}}
#' @seealso \code{\link{calculateRefKendrickMassDefect}}
#'
#' @export
checkRmkd <- function(rmkd, error = 0.15) {

  # check if rmkd is in error range
  remainder <- rmkd %% -1

  if(abs(remainder) < error | abs(1 + remainder) < error) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}
