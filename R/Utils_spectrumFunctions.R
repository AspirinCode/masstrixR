bin_Spectrum <- function(object, binSize = 1L,
                         breaks = seq(floor(min(mz(object))),
                                      ceiling(max(mz(object))),
                                      by = binSize),
                         fun = sum,
                         msLevel.) {
  ## If msLevel. not missing, perform the trimming only if the msLevel
  ## of the spectrum matches (any of) the specified msLevels.
  if (!missing(msLevel.)) {
    if (!(msLevel(object) %in% msLevel.))
      return(object)
  }
  bins <- .bin_values(object@intensity, object@mz, binSize = binSize,
                      breaks = breaks, fun = fun)
  object@mz <- bins$mids
  object@intensity <- bins$x
  object@tic <- sum(object@intensity)
  object@peaksCount <- length(object@mz)
  if (validObject(object))
    return(object)
}

bin_Spectra <- function(object1, object2, binSize = 1L,
                        breaks = seq(floor(min(c(mz(object1), mz(object2)))),
                                     ceiling(max(c(mz(object1), mz(object2)))),
                                     by = binSize)) {
  breaks <- .fix_breaks(breaks, range(mz(object1), mz(object2)))
  list(bin_Spectrum(object1, breaks = breaks),
       bin_Spectrum(object2, breaks = breaks))
}

.fix_breaks <- function(brks, rng) {
  ## Assuming breaks being sorted.
  if (brks[length(brks)] <= rng[2])
    brks <- c(brks, max((rng[2] + 1e-6),
                        brks[length(brks)] + mean(diff(brks))))
  brks
}

.bin_values <- function(x, toBin, binSize = 1, breaks = seq(floor(min(toBin)),
                                                            ceiling(max(toBin)),
                                                            by = binSize),
                        fun = max) {
  if (length(x) != length(toBin))
    stop("lengths of 'x' and 'toBin' have to match.")
  fun <- match.fun(fun)
  breaks <- .fix_breaks(breaks, range(toBin))
  nbrks <- length(breaks)
  idx <- findInterval(toBin, breaks)
  ## Ensure that indices are within breaks.
  idx[which(idx < 1L)] <- 1L
  idx[which(idx >= nbrks)] <- nbrks - 1L

  ints <- double(nbrks - 1L)
  ints[unique(idx)] <- unlist(lapply(base::split(x, idx), fun),
                              use.names = FALSE)
  list(x = ints, mids = (breaks[-nbrks] + breaks[-1L]) / 2L)
}

alignSpectra <- function(x, y, mzTol = 0.005, treshold = 0.01, returnSpectra = FALSE) {

  top <- data.frame(mz = mz(x), intensity = intensity(x))
  bottom <- data.frame(mz = mz(y), intensity = intensity(y))

  top$intensity <- top$intensity / max(top$intensity) * 100
  bottom$intensity <- bottom$intensity / max(bottom$intensity) * 100

  top <- top[which(top$intensity > treshold),]
  bottom <- bottom[which(bottom$intensity > treshold),]

  ## align the m/z axis of the two spectra, the bottom spectrum is used as the reference

  for(i in 1:nrow(bottom)) {
    top[,1][bottom[,1][i] >= top[,1] - mzTol & bottom[,1][i] <= top[,1] + mzTol] <- bottom[,1][i]
  }

  alignment <- merge(top, bottom, by = 1, all = TRUE)
  if(length(unique(alignment[,1])) != length(alignment[,1])) {
    warning("the m/z tolerance is set too high")
  }

  alignment[,c(2,3)][is.na(alignment[,c(2,3)])] <- 0   # convert NAs to zero (R-Help, Sept. 15, 2004, John Fox)
  names(alignment) <- c("mz", "intensity.top", "intensity.bottom")


  if(returnSpectra) {
    # correct masses and intensites
    x@mz <- alignment$mz
    x@intensity <- alignment$intensity.top

    y@mz <- alignment$mz
    y@intensity <- alignment$intensity.bottom

    return(list(x,y))
  } else {
    return(alignment)
  }


}
