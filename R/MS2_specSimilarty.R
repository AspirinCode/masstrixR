#' Calculate forward Dotproduct
#'
#'
#' @export
forwardDotProduct <- function(x, y, align = FALSE, mzTol = 0.005, mzTolType = "abs", binwidth = 1L, ...) {

  if(align) {

    #use aligned spectra
    alignedSpectra <- alignSpectra(x, y, mzTol = mzTol, mzTolType = mzTolType)

    # calculate dot product
    dotproduct <- dotproduct(alignedSpectra$intensity.top, alignedSpectra$intensity.bottom)

    # return result
    return(dotproduct)

  } else {

    # use binned spectra
    binnedSpectra <- bin_Spectra(x, y, ...)

    # calculate dot product
    dotproduct <- dotproduct(binnedSpectra$intensity.top, binnedSpectra$intensity.bottom)

    # return result
    return(dotproduct)
  }
}


#' Calculate reverse Dotproduct
#'
#'
#' @export
reverseDotProduct <- function(x, y, align = FALSE, mzTol = 0.005, mzTolType = "abs", binwidth = 1L, ...) {

  if(align) {

    #use aligned spectra
    alignedSpectra <- alignSpectra(x, y, mzTol = mzTol, mzTolType = mzTolType)

    # use only peals with are present in spec2
    alignedSpectra <- alignedSpectra[which(alignedSpectra$intensity.bottom >0), ]

    # calculate dot product
    dotproduct <- dotproduct(alignedSpectra$intensity.top, alignedSpectra$intensity.bottom)

    # return result
    return(dotproduct)

  } else {

    # use binned spectra
    binnedSpectra <- bin_Spectra(x, y, ...)

    # use only peals with are present in spec2
    binnedSpectra <- binnedSpectra[which(binnedSpectra$intensity.bottom >0), ]

    # calculate dot product
    dotproduct <- dotproduct(binnedSpectra$intensity.top, binnedSpectra$intensity.bottom)

    # return result
    return(dotproduct)
  }
}


#' Calculate common peaks
#'
#'
#' @export
commonPeaks <- function(x, y, align = FALSE, mzTol = 0.005, mzTolType = "abs", binwidth = 1L, ...) {

  if(align) {

    #use aligned spectra
    alignedSpectra <- alignSpectra(x, y, mzTol = mzTol, mzTolType = mzTolType)

    commonPeaks <- nrow(alignedSpectra[which(alignedSpectra$intensity.top > 0 & alignedSpectra$intensity.bottom > 0),])

    # return number of common peaks
    return(commonPeaks)

  } else {

    # use binned spectra
    binnedSpectra <- bin_Spectra(x, y, ...)

    commonPeaks <- nrow(binnedSpectra[which(binnedSpectra$intensity.top > 0 & binnedSpectra$intensity.bottom > 0),])

    # return number of common peak
    return(commonPeaks)

  }
}

#' Mass shift modified Dotproduct (GNPS)
#'
#' Inspired by GNPS
#'
#'
massShiftForwadDotProduct <- function() {

  # TODO implement massshift alignment

}

#' XRank
#'
#' Inspired b Weizmass
#'
#'
xRank <- function() {

  # TODO implement xRank function

}


#' Helper function to calculate dot product between two intensity vectors on same m/z scale (binned or aligned)
dotproduct <- function(x, y) {
  as.vector(x %*% y) / (sqrt(sum(x*x)) * sqrt(sum(y*y)))
}
