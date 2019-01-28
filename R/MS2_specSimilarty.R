

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
