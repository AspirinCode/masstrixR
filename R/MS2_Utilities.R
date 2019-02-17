#' This function allows to search for a specific product ion in MS2 spectra
#'
#' @param ms2spectrum A Spectrum2 in which the production ion shall be searched
#' @param productIonMz m/z value of product ion that shall be search for
#' @param mzTol Maximum allowed tolerance in Da or ppm
#' @param mzTolType Defines the error type for m/z search, "abs" is used for absolute mass error, "ppm" for relative error
#'
#' @import MSnbase
#' @export
containsProductIon <- function(ms2spectrum, productIonMz, mzTol = 0.005, mzTolType = "abs", multiplePi = "any") {

  # get mz values to search in
  mz <- mz(ms2spectrum)

  containsPi <- FALSE

  # check if multiple masses are defined
  if(length(productIonMz) > 1) {

    # generate list of checks
    booleanVec <- unlist(lapply(productIonMz, function(x, y) {
      if(mzTolType == "abs") {

        any(abs(x - y) < mzTol)

      } else if(mzTolType == "ppm") {

        any(abs((x - y) / y * 1e6) < mzTol)

      }
    }, y = mz))

    if(multiplePi == "any") {

      containsPi <- any(booleanVec)

    } else if(multiplePi == "all") {

      containsPi <- all(booleanVec)

    }

  } else {

    # search for single mass
    if(mzTolType == "abs") {

      containsPi <- any(abs(mz - productIonMz) < mzTol)

    } else if(mzTolType == "ppm") {

      containsPi <- any(abs((mz - productIonMz) / productIonMz * 1e6) < mzTol)

    } else {

      stop("unknown mzTolType")

    }
  }

  #return result
  return(containsPi)

}

#' This function allows to search for a fragment m/z in MS2 spectra that match to a specified neutral from the precusor
#'
#' @param ms2spectrum A Spectrum2 in which the production ion shall be searched
#' @param neutralLossMass Mass of neutral loss that shal be searched
#' @param mzTol Maximum allowed tolerance in Da or ppm
#' @param mzTolType Defines the error type for m/z search, "abs" is used for absolute mass error, "ppm" for relative error
#'
#' @import MSnbase
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

#' This function allows to search for a specific mass difference between two fragments
#'
#' @param ms2spectrum A Spectrum2 in which the production ion shall be searched
#' @param fragmentMassDifference Mass of neutral loss that shal be searched
#' @param mzTol Maximum allowed tolerance in Da or ppm
#' @param mzTolType Defines the error type for m/z search, "abs" is used for absolute mass error, "ppm" for relative error
#' @import MSnbase
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

# remove???
findCommonFragments <- function() {

}

#'
#' @import MSnbase
#' @export
makeDiffSpectra <- function(x, y, align = TRUE, mzTol = 0.005, mzTolType = "abs", binwidth = 1L, ...) {

  if(align) {

    #use aligned spectra
    alignedSpectra <- alignSpectra(x, y, mzTol = mzTol, mzTolType = mzTolType)

    # select only unique peaks
    uniquePeaks <- alignedSpectra[-which(alignedSpectra$intensity.top > 0 & alignedSpectra$intensity.bottom > 0),]

    # make Spectrum2 oboject from data
    diffSpec <- new("Spectrum2",
                   merged = 0,
                   precScanNum = as.integer(1),
                   precursorMz = precursorMz(x),
                   precursorIntensity = 100,
                   precursorCharge = precursorCharge(x),
                   mz = uniquePeaks$mz,
                   intensity = uniquePeaks$intensity.top,
                   collisionEnergy = collisionEnergy(x),
                   centroided = TRUE)

    # return updated spectrum
    return(diffSpec)

  } else {

    # use binned spectra
    binnedSpectra <- bin_Spectra(x, y, ...)

    # select only unique peaks
    uniquePeaks <- binnedSpectra[-which(binnedSpectra$intensity.top > 0 & binnedSpectra$intensity.bottom > 0),]

    # make Spectrum2 oboject from data
    diffSpec <- new("Spectrum2",
                    merged = 0,
                    precScanNum = as.integer(1),
                    precursorMz = precursorMz(x),
                    precursorIntensity = 100,
                    precursorCharge = precursorCharge(x),
                    mz = uniquePeaks$mz,
                    intensity = uniquePeaks$intensity.top,
                    collisionEnergy = collisionEnergy(x),
                    centroided = TRUE)

    # return updated spectrum
    return(diffSpec)

  }

}

#'
#' @import MSnbase
#' @import xcms
#' @export
filterIsotopePeaks <- function(x, mzTol = 0.005, mzTolType = "abs") {

  #create matrix for xcmsSet
  mat<-matrix(0,ncol=12,nrow=length(mz(x)))

  colnames(mat)<-c("mz","mzmin","mzmax","rt","rtmin","rtmax","into","intf","maxo","maxf","i","sn")
  mat[,"mz"]<-mz(x)
  mat[,"into"]<-intensity(x)
  mat[,"rt"]<-1
  mat[,"rtmin"]<-mat[,"rt"]-0.1
  mat[,"intf"]<-intensity(x)
  mat[,"rtmax"]<-mat[,"rt"]+0.1

  #make new xcmsSEt
  test<-new("xcmsSet")
  peaks(test)<-mat
  sampnames(test)<-c("samp1")
  test@filepaths="test"

  # peform CAMERA based annotation
  xa <- CAMERA::xsAnnotate(test)
  xa <- CAMERA::groupFWHM(xa)
  xa <- CAMERA::findIsotopes(xa,intval="into",mzabs=mzTol,filter=F)

  # get spectra an filter values
  spec <- CAMERA::getpspectra(xa,1)[,c("mz","into","isotopes")]

  spec <- spec[-grep("M\\+", spec[,3]),]

  # make Spectrum2 oboject from data
  cleanedSpec <- new("Spectrum2",
                     merged = merged(x),
                     precScanNum = scanIndex(x),
                     precursorMz = precursorMz(x),
                     precursorIntensity = precursorIntensity(x),
                     precursorCharge = precursorCharge(x),
                     mz = spec[,1],
                     intensity = spec[,2],
                     collisionEnergy = collisionEnergy(x),
                     centroided = centroided(x))

  # return isotope cleaned spectrum
  return(cleanedSpec)
}


