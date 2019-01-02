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
