# plotting function
#'
#'
#' @export
makeMirrorPlot <- function(x, y, align = FALSE, mzTol = 0.005, treshold = 0.01, title = "Mirrorplot", ...) {

  require(ggplot2)
  require(gridExtra)
  require(grid)

  if(align) {

    #use aligned spectra
    alignedSpectra <- alignSpectra(x, y, mzTol = mzTol, treshold = treshold)
    commonPeaks <- alignedSpectra[which(alignedSpectra$intensity.top > 0 & alignedSpectra$intensity.bottom > 0),]

    #x.precursorMz <- round(precursorMz(x), 4)
    #y.precursorMz <- round(precursorMz(y), 4)

    # df <- data.frame(x = c(min(alignedSpectra$mz), min(alignedSpectra$mz)),
    #                  y = c(110, -100),
    #                  label = c(x.precursorMz, y.precursorMz))
    #
    # print(df)

    if(nrow(commonPeaks) > 0) {

      p1 <- ggplot(alignedSpectra, aes(x = mz, ymin = 0, ymax = intensity.top)) +
        geom_linerange(colour = "blue") +
        geom_linerange(data = alignedSpectra, aes(x = mz, ymax = 0, ymin = intensity.bottom * -1), colour = "red") +
        geom_point(data = commonPeaks, aes(x = mz, y = intensity.top + 5), shape = 25, colour = "black", fill = "blue") +
        geom_point(data = commonPeaks, aes(x = mz, y = intensity.bottom * - 1 -5), shape = 24, colour = "black", fill = "red") +
        ggtitle(title) +
        xlab("m/z") + ylab("normalized intensity") +
        theme_bw()

      p1


    } else {

      p1 <- ggplot(alignedSpectra, aes(x = mz, ymin = 0, ymax = intensity.top)) +
        geom_linerange(colour = "blue") +
        geom_linerange(data = alignedSpectra, aes(x = mz, ymax = 0, ymin = intensity.bottom * -1), colour = "red") +
        ggtitle(title) +
        xlab("m/z") + ylab("normalized intensity") +
        theme_bw()

      p1

    }
  }
}
