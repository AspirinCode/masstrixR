#' Calculate similarity of isotope pattern
#'
#' based on S-ratio from MS-Dial
#'
#' @export
isoPatternSimilarity <- function(x, y, mzTol = 0.005, ...) {

  # align spectra and sort according to mass
  alignedSpectra <- alignSpectra(x, y, mzTol = mzTol)
  alignedSpectra <- alignedSpectra[order(alignedSpectra$mz),]

  print(alignedSpectra)

  # prepare sum
  sum <- 0

  # iterate through all peaks and compare neighbors
  for(i in 1:(nrow(alignedSpectra) - 1)) {

    # calculate ratios for neighboring isotopes
    rLib <- alignedSpectra$intensity.bottom[i + 1] / alignedSpectra$intensity.bottom[i]
    rAct <- alignedSpectra$intensity.top[i + 1] / alignedSpectra$intensity.top[i]

    print(rLib)
    print(rAct)

    # remove infinite values
    if(rLib == Inf) {
      rLib <- 0
    }

    if(rAct == Inf) {
      rAct <- 0
    }

    # make sum of absolute differences
    sum <- sum + abs(rAct - rLib)

  }

  # calculate sRatio
  sRatio <- 1 - sum

  return(sRatio)
}

#' Calculate forward Dotproduct
#'
#'
#' @export
forwardDotProduct <- function(x, y, align = FALSE, mzTol = 0.005, treshold = 0.01, ...) {

  if(align) {

    #use aligned spectra
    alignedSpectra <- alignSpectra(x, y, mzTol = mzTol, treshold = treshold)

    # calculate dot product
    dotproduct <- dotproduct(alignedSpectra$intensity.top, alignedSpectra$intensity.bottom)

    # return result
    return(dotproduct)

  } else {

    # use binned spectra
    binnedSpectra <- bin_Spectra(x, y, ...)
    inten <- as.data.frame(lapply(binnedSpectra, intensity))
    names(inten) <- c("spec1", "spec2")

    # calculate dot product
    dotproduct <- dotproduct(inten$spec1, inten$spec2)

    # return result
    return(dotproduct)
  }
}

#' Calculate reverse Dotproduct
#'
#'
#' @export
reverseDotProduct <- function(x, y, align = FALSE, mzTol = 0.005, treshold = 0.01, ...) {

  if(align) {

    #use aligned spectra
    alignedSpectra <- alignSpectra(x, y, mzTol = mzTol, treshold = treshold)

    # use only peals with are present in spec2
    alignedSpectra <- alignedSpectra[which(alignedSpectra$intensity.bottom >0), ]

    # calculate dot product
    dotproduct <- dotproduct(alignedSpectra$intensity.top, alignedSpectra$intensity.bottom)

    # return result
    return(dotproduct)

  } else {

    # use binned spectra
    binnedSpectra <- bin_Spectra(x, y, ...)
    inten <- as.data.frame(lapply(binnedSpectra, intensity))
    names(inten) <- c("spec1", "spec2")

    # use only peals with are present in spec2
    inten <- inten[which(inten$spec2 > 0),]

    # calculate dot product
    dotproduct <- dotproduct(inten$spec1, inten$spec2)

    # return result
    return(dotproduct)
  }
}


#' Calculate reverse Dotproduct
#'
#'
#' @export
commonPeaks <- function(x, y, align = FALSE, mzTol = 0.005, treshold = 0.01, ...) {

  if(align) {

    #use aligned spectra
    alignedSpectra <- alignSpectra(x, y, mzTol = mzTol, treshold = treshold)

    commonPeaks <- nrow(alignedSpectra[which(alignedSpectra$intensity.top > 0 & alignedSpectra$intensity.bottom > 0),])

    return(commonPeaks)



  } else {

    # use binned spectra
    binnedSpectra <- bin_Spectra(x, y, ...)
    inten <- as.data.frame(lapply(binnedSpectra, intensity))
    names(inten) <- c("spec1", "spec2")

    commonPeaks <- nrow(inten[which(inten$spec1 > 0 & inten$spec2 > 0),])

    return(commonPeaks)

  }
}

#' Mass shift modified Dotproduct (GNPS)
#'
#' Inspired by GNPS
#'
#' @export
massShiftForwadDotProduct <- function() {

}

#' XRank
#'
#' Inspired b Weizmass
#'
#' @export
xRank <- function() {

}


#' Helper function to calculate dot product between two intensity vectors on same m/z scale (binned or aligned)
dotproduct <- function(x, y) {
  as.vector(x %*% y) / (sqrt(sum(x*x)) * sqrt(sum(y*y)))
}
