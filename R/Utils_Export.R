writeSirius <- function(isoPattern, ms2spectrum, file) {

  con <- file(file, "w")
  on.exit(close(con))

  # custom cat function
  .cat <- function(..., file = con, sep = "", append = TRUE) {
    cat(..., file = file, sep = sep, append = append)
  }

  # write isopattern (MS1 level)
  .cat("\nBEGIN IONS\n")

  .cat("PEPMASS=", isoPattern@monoMz, "\n")

  .cat("MSLEVEL=1\n",
       "CHARGE=1+\n",
       "NAME=",isoPattern@featureId)

  .cat("\n", paste(mz(isoPattern), intensity(isoPattern), collapse = "\n"))
  .cat("\nEND IONS\n")

  # wirte fragmentation spectrum
  .cat("\nBEGIN IONS\n")

  .cat("PEPMASS=", isoPattern@monoMz, "\n")

  .cat("MSLEVEL=2\n",
       "CHARGE=1+\n",
       "NAME=",isoPattern@featureId)

  .cat("\n", paste(mz(ms2spectrum), intensity(ms2spectrum), collapse = "\n"))
  .cat("\nEND IONS\n")


}
