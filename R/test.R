# load required libraries
library(masstrixR)
library(Rdisop)

# Glucose as example, [M+Na]+
masses <- c(203.052609, 204.056051, 205.057227)
intensities <- c(100, 6.868, 1.434)

# make MS1 spectrum
measured<- new(
  "Spectrum1",
  mz = masses,
  intensity = intensities,
  centroided = TRUE
)

# Decompose all possible formula
molecules <- Rdisop::decomposeIsotopes(masses, intensities, mzabs = 0.005, z = 1, elements = Rdisop::initializeCHNOPSNaK())
cbind(Rdisop::getFormula(molecules), Rdisop::getScore(molecules), Rdisop::getValid(molecules))

formulas <- getFormula(molecules)


df <- data.frame()

for(formula in formulas) {

  print(formula)
  isotopePattern <- predictIsoPattern(formula, charge = 1, plotit = FALSE)
  isotopePattern@intensity <- isotopePattern@intensity / sum(isotopePattern@intensity)

  # compare spectra of enviPat vs Rdisop
  # TODO add comparison here

  #print(isotopePattern)

  # make new molecule from data
  molecule <- list(formula = formula,
                    score = 1,
                    exactmass = max(mz(isotopePattern)),
                    charge = 1,
                    parity = "e",
                    valid = "valid",
                    DBE = 1,
                    isotopes = list(matrix(c(mz(isotopePattern), intensity(isotopePattern)), nrow = 2, byrow = TRUE)))

  p1 <- makeMirrorPlot(measured, isotopePattern, align = TRUE, mzTol = 0.001)
  plot(p1)

  # iso score
  print(isotopeScore(molecule, mz(measured), intensity(measured), z = 1))
  score <- isotopeScore(molecule, mz(measured), intensity(measured), z = 1)

  # make DF for results
  df <- rbind.data.frame(df, cbind.data.frame(formula = as.character(formula),
                                              score = score), stringsAsFactors = FALSE)

}


df$normScore <- df$score / sum(df$score)


