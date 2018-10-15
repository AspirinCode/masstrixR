
plot.isoPatternSpectrum <- function(obj, ...) {
  plot.default(obj@mz, obj@intensity, type = "h", ylim = c(0, max(obj@intensity)))
}


