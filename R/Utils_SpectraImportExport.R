writeSirius <- function(spectraList, file) {

  ## sanity checks
  if(!class(isoPatternSpectra)[1] == "Spectra") {
    stop("spectraList should be of class Spectra")
  }

  ## get unique cluster ids in the spectraList object
  exportClusters <- unique(spectraList@elementMetadata$cluster_id)

  ## prepare for writing
  con <- file(file, "w")
  on.exit(close(con))

  # custom cat function
  .cat <- function(..., file = con, sep = "", append = TRUE) {
    cat(..., file = file, sep = sep, append = append)
  }

  # iterate over all cluster
  for(exportCluster in exportClusters) {

    print(exportCluster)

    # get all spectra associated with this cluster
    clusterSpectra <- spectraList[which(spectraList@elementMetadata$cluster_id == exportCluster)]

    for(i in 1:length(clusterSpectra)) {

      # write header
      .cat("\nBEGIN IONS\n")
      .cat("PEPMASS=", clusterSpectra[i]@elementMetadata$mz, "\n")

      # dependent on type MSLEVEL=1 or 2
      if(clusterSpectra[i]@elementMetadata$type == "iso") {
        .cat("MSLEVEL=1\n")
      } else if(clusterSpectra[i]@elementMetadata$type == "ms2") {
        .cat("MSLEVEL=2\n")
      }

      .cat("CHARGE=1+\n")
      .cat("NAME=", clusterSpectra[i]@elementMetadata$cluster_id, "\n")
      .cat(paste(mz(clusterSpectra[[i]]), intensity(clusterSpectra[[i]]), collapse = "\n"))
      .cat("\nEND IONS\n")
    }

  }
}


extractMgfSpectrum2Info <- function(mgf, centroided, addFields = NULL) {

  # grep description
  desc.idx <- grep("=", mgf)
  desc <- mgf[desc.idx]
  spec <- mgf[-desc.idx]

  ms <- do.call(rbind, strsplit(spec, "[[:space:]]+"))
  mode(ms) <- "double"

  if(!length(ms)) {
    ms <- matrix(numeric(), ncol = 2L)
  }

  r <- regexpr("=", desc, fixed = TRUE)
  desc <- setNames(substring(desc, r + 1L, nchar(desc)), substring(desc, 1L, r - 1L))
  fdata <- desc

  desc[c("PEPMASSMZ", "PEPMASSINT")] <- strsplit(desc["PEPMASS"], "[[:space:]]+")[[1L]][1:2]

  # select only values of interest and convert to numeric (base fields)
  desc["CHARGE"] <- sub("[+-]", "", desc["CHARGE"])
  voi <- c("RTINSECONDS", "CHARGE", "SCANS", "PEPMASSMZ", "PEPMASSINT")
  desc.base <- setNames(as.numeric(desc[voi]), voi)
  desc.base[is.na(desc.base[voi])] <- 0L
  cat(".")

  if(!is.null(addFields)) {

    desc.add <- setNames(desc[addFields], addFields)
    cat(".")

  } else {
    desc.add <- NA
  }


  # create spectrum
  sp <- Spectrum2_mz_sorted(rt = unname(desc["RTINSECONDS"]),
                            scanIndex = unname(as.integer(desc["SCANS"])),
                            precursorMz = unname(desc["PEPMASSMZ"]),
                            precursorIntensity = unname(desc["PEPMASSINT"]),
                            precursorCharge = unname(as.integer(desc["CHARGE"])),
                            mz = ms[, 1L],
                            intensity = ms[, 2L],
                            fromFile = 1L,
                            centroided = centroided)

  # return values
  return(list(spectrum = sp, fdata = fdata, addData = desc.add))
}
