#' Prediction of isotope pattern
#'
#' This function predicts the isotope pattern for a given ion formula and adduct definition.
#'
#' @param ionFormula The formula of the ion for which the isotopic pattern shall be calculated, e.g. C6H12O6Na for the M+Na adduct of glucose
#' @param charge Charge of the ion.
#' @param plotit Boolean value if the isotope pattern shall be plotted
#' @param treshold Intensity treshold for lowest isotope
#' @param resolution resolution of the MS
#' @return Returns a MSnbase Spectrum 1 object containing the calculated isotope pattern
#' @examples
#' predictIsoPattern("C6H12O6Na", charge = 1)
#' @export
predictIsoPattern <-
  function(ionFormula,
             charge,
             plotit = FALSE,
             treshold = 0.01,
             resolution = 50000) {
    require(MSnbase)

    # pre-check generated ion formula
    checked <- enviPat::check_chemform(isotopes, ionFormula)

    # create isotope pattern
    centro <- enviPat::isowrap(
      isotopes,
      checked,
      resmass = FALSE,
      resolution = resolution,
      nknots = 4,
      spar = 0.2,
      threshold = 0.01,
      charge = charge,
      emass = 0.00054858,
      algo = 2,
      ppm = FALSE,
      dmz = "get",
      # retrieve dm from R=m/dm
      frac = 1 / 4,
      env = "Gaussian",
      detect = "centroid",
      plotit = plotit
    )

    # create Spectrum1 object from it
    centro <- as.data.frame(centro)
    colnames(centro) <- c("mz", "int")
    isotopeSpectrum <- new(
      "Spectrum1",
      mz = centro$mz,
      intensity = centro$int,
      centroided = TRUE
    )

    # return values
    return(isotopeSpectrum)
  }


#' Reconstruct isotope pattern from gda files
#'
#'
reconstructIsoPattern <- function(peaks, sampleNames) {

  clusterList <- unique(peaks$Cluster)

  isoPatternList <- new("Spectra")

  for(cluster in clusterList) {

    # get m/z values
    mz <- peaks$m.z[which(peaks$Cluster == cluster)]

    for(sample in sampleNames) {

      # get intensities for each sample
      int <- peaks[which(peaks$Cluster == cluster), sample]
      int[is.na(int)] <- 0.0

      ms1spec <- new("Spectrum1",
                     mz = mz,
                     intensity = int,
                     centroided = TRUE)

      # create new Spectra object
      isoPattern <- Spectra(ms1spec)

      # add annotation
      mcols(isoPattern)$Cluster <- cluster
      mcols(isoPattern)$sample <- sample

      isoPatternList <- append(isoPatternList, isoPattern)
    }

  }

  return(isoPatternList)
}

#' Function to score isotope pattern
#'
scoreIsoPattern <- function(filteredAnnotationResult, measuredIsotopePattern, returnPredicted = FALSE) {

  adductCalc <- getAdductCalc()

  predictedIsotopePatternList <- new("Spectra")

  # iterate over all results
  for(i in 1:nrow(filteredAnnotationResult)) {

    # get charge of annotation
    adduct <- filteredAnnotationResult$adductType[i]
    charge <- 1
    ionFormula <- filteredAnnotationResult$ionFormula[i]


    # predict isotope pattern
    predictedIsotopePattern <- predictIsoPattern(ionFormula, charge = charge, plotit = FALSE)
    predictedIsotopePattern@intensity <- predictedIsotopePattern@intensity / sum(predictedIsotopePattern@intensity)

    # make new spectra object and add predicted spectrum
    spectraClipboard <- Spectra(predictedIsotopePattern)
    mcols(spectraClipboard)$cluster <- filteredAnnotationResults$cluster[i]
    mcols(spectraClipboard)$ionFormula <- ionFormula
    mcols(spectraClipboard)$adduct <- adduct

    #append to list
    predictedIsotopePatternList <- append(predictedIsotopePatternList, spectraClipboard)

    # make new molecule from data in format for Rdisop
    molecule <- list(formula = formula,
                     score = 1,
                     exactmass = max(mz(isotopePattern)),
                     charge = charge,
                     parity = "e",
                     valid = "valid",
                     DBE = 1,
                     isotopes = list(matrix(c(mz(predictedIsotopePattern), intensity(predictedIsotopePattern)), nrow = 2, byrow = TRUE)))

    # add score to data table
    filteredAnnotationResult$score[i] <- isotopeScore(molecule, mz(measuredIsotopePattern), intensity(measuredIsotopePattern), z = charge)
  }

  # calculate final score
  filteredAnnotationResult$score <- filteredAnnotationResult$score / sum(filteredAnnotationResult$score)

  #make decision based on return option
  if(returnPredicted) {

    #return both results as list
    return(list(filteredAnnotationResult, predictedIsotopePatternList))

  } else {

  # return result
  return(filteredAnnotationResult)

  }
}


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
