# plotting function
#'
#'
#' @import ggplot2
#' @import grid
#' @export
makeMirrorPlot <- function(x, y, align = FALSE, plotIt = FALSE, mzTol = 0.005, treshold = 0.01, title = "Mirrorplot", ...) {

  if(align) {

    #use aligned spectra
    alignedSpectra <- alignSpectra(x, y, mzTol = mzTol)
    commonPeaks <- alignedSpectra[which(alignedSpectra$intensity.top > 0 & alignedSpectra$intensity.bottom > 0),]

    noPeaks_x <- length(mz(x))
    noPeaks_y <- length(mz(y))

    if(nrow(commonPeaks) > 0) {

      p1 <- ggplot() +

        # upper spectrum
        geom_linerange(data = alignedSpectra, aes(x = mz, ymin = 0, ymax = intensity.top), colour = "blue") +
        annotate("text", x = min(alignedSpectra$mz) - 5, y = 120, label = paste0("Precursor m/z: ", round(precursorMz(x), 4), " / Common Peaks: ", nrow(commonPeaks), " / Total Peaks: ", noPeaks_x), hjust = 0) +

        # lower spectrum
        geom_linerange(data = alignedSpectra, aes(x = mz, ymax = 0, ymin = intensity.bottom * -1), colour = "red") +
        annotate("text", x = min(alignedSpectra$mz) - 5, y = -120, label = paste0("Precursor m/z: ", round(precursorMz(y), 4), " / Common Peaks: ", nrow(commonPeaks), " / Total Peaks: ", noPeaks_y), hjust = 0) +

        # common peaks marker
        geom_point(data = commonPeaks, aes(x = mz, y = intensity.top + 5), shape = 25, colour = "black", fill = "blue") +
        geom_point(data = commonPeaks, aes(x = mz, y = intensity.bottom * - 1 -5), shape = 24, colour = "black", fill = "red") +

        # title
        ggtitle(title) +

        # scaling
        scale_x_continuous(limits = c(min(alignedSpectra$mz) - 5, max(alignedSpectra$mz) + 5)) +
        scale_y_continuous(breaks = c(-100, -50, 0, 50, 100)) +

        # axis
        xlab("m/z") + ylab("normalized intensity") +

        # theme
        theme_bw() +
        theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black", linetype = "solid"))

      if(plotIt) {
        plot(p1)
      }


    } else {

      p1 <- ggplot() +

        # upper spectrum
        geom_linerange(data = alignedSpectra, aes(x = mz, ymin = 0, ymax = intensity.top), colour = "blue") +
        annotate("text", x = min(alignedSpectra$mz) - 5, y = 120, label = paste0("Precursor m/z: ", round(precursorMz(x), 4), " / Common Peaks: ", nrow(commonPeaks)), hjust = 0) +

        # lower spectrum
        geom_linerange(data = alignedSpectra, aes(x = mz, ymax = 0, ymin = intensity.bottom * -1), colour = "red") +
        annotate("text", x = min(alignedSpectra$mz) - 5, y = -120, label = paste0("Precursor m/z: ", round(precursorMz(y), 4), " / Common Peaks: ", nrow(commonPeaks)), hjust = 0) +

        # title
        ggtitle(title) +

        # scaling
        scale_x_continuous(limits = c(min(alignedSpectra$mz) - 5, max(alignedSpectra$mz) + 5)) +
        scale_y_continuous(breaks = c(-100, -50, 0, 50, 100)) +

        # axis
        xlab("m/z") + ylab("normalized intensity") +

        # theme
        theme_bw() +
        theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black", linetype = "solid"))

      if(plotIt) {
        plot(p1)
      }

    }
  } else {
    NULL
  }
}

# plotting function
#'
#'
#'
#' @import ggplot2
#' @import grid
#' @export
plotSpectrum <- function(x, highlight = FALSE, highlightMz = NULL, mzTol = 0.005, plotIt = FALSE, ...) {

  # check for class
  if(class(x) == "Spectra") {
    x <- x[[1]]
  }


  if(highlight) {

    #check if m/z values for highlighting are supplied
    if(is.null(highlightMz)) {
      stop("No m/z to highlight supplied!")
    }

    # make spectrum object for highlighting
    highlightSpectrum <- new("Spectrum1",
                             mz = highlightMz,
                             intensity = rep(1, length(highlightMz)))

    #use aligned spectra
    alignedSpectra <- alignSpectra(highlightSpectrum, x, mzTol = mzTol)
    commonPeaks <- alignedSpectra[which(alignedSpectra$intensity.top > 0 & alignedSpectra$intensity.bottom > 0),]

    if(nrow(commonPeaks) > 0) {

      p1 <- ggplot() +

        # upper spectrum
        geom_linerange(data = alignedSpectra, aes(x = mz, ymin = 0, ymax = intensity.bottom), colour = "blue") +

        # common peaks marker
        geom_point(data = commonPeaks, aes(x = mz, y = intensity.bottom + 5), shape = 25, colour = "black", fill = "blue") +

        # title
        ggtitle(paste0("Precursor m/z: ", round(precursorMz(x), 4), " / Common Peaks: ", nrow(commonPeaks))) +

        # scaling
        scale_x_continuous(limits = c(min(alignedSpectra$mz) - 5, max(alignedSpectra$mz) + 5)) +
        scale_y_continuous(breaks = c(0, 50, 100)) +

        # axis
        xlab("m/z") + ylab("normalized intensity") +

        # theme
        theme_bw() +
        theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black", linetype = "solid"))

      if(plotIt) {
        plot(p1)
      }

    } else {

      p1 <- ggplot() +

        # upper spectrum
        geom_linerange(data = alignedSpectra, aes(x = mz, ymin = 0, ymax = intensity.top), colour = "blue") +

        # title
        ggtitle(paste0("Precursor m/z: ", round(precursorMz(x), 4))) +

        # scaling
        scale_x_continuous(limits = c(min(alignedSpectra$mz) - 5, max(alignedSpectra$mz) + 5)) +
        scale_y_continuous(breaks = c(0, 50, 100)) +

        # axis
        xlab("m/z") + ylab("normalized intensity") +

        # theme
        theme_bw() +
        theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black", linetype = "solid"))

      if(plotIt) {
        plot(p1)
      }

    }
  } else {

    alignedSpectra <- data.frame(mz = mz(x),
                                 intensity.top = intensity(x) / max(intensity(x) * 100))

    p1 <- ggplot() +

      # upper spectrum
      geom_linerange(data = alignedSpectra, aes(x = mz, ymin = 0, ymax = intensity.top), colour = "blue") +

      # title
      ggtitle(paste0("Precursor m/z: ", round(precursorMz(x), 4))) +

      # scaling
      scale_x_continuous(limits = c(min(alignedSpectra$mz) - 5, max(alignedSpectra$mz) + 5)) +
      scale_y_continuous(breaks = c(0, 50, 100)) +

      # axis
      xlab("m/z") + ylab("normalized intensity") +

      # theme
      theme_bw() +
      theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black", linetype = "solid"))

    if(plotIt) {
      plot(p1)
    }

  }
}
