#' This function returns a list with adducts and the required multiplicative and additive part to calculate ion m/z from exact masses. Adduts are taken from http://fiehnlab.ucdavis.edu/staff/kind/metabolomics/ms-adduct-calculator
#'
#' @return Returns a named list with all mulitpliers and additive parts for positive and negative mode
#' @examples
#' adductCalc <- getAdductCalc()
#' @export
getAdductCalc <- function() {

  # create list with multiplicative and additive part plus formula parts to add and substract ("C0" is a fake formula for calculation)
  adductCalc <- list(
<<<<<<< HEAD
    "M+3H" = c(1 / 3, (3 * 1.007276) / 3, "H3", "C0", 3),
    "M+2H+Na" = c(1 / 3, (2 * 1.007276 + 22.989218) / 3, "H2Na", "C0", 3),
    "M+H+2Na" = c(1 / 3, (1.007276 + 2 * 22.989218) / 3, "HNa2", "C0", 3),
    "M+3Na" = c(1 / 3, (3 * 22.989218) / 3, "Na3", "C0", 3),
    "M+2H" = c(1 / 2, (2 * 1.007276) / 2, "H2", "C0", 2),
    "M+H+NH4" = c(1 / 2, (1.007276 + 18.033823) / 2, "NH5", "C0", 2),
    "M+H+K" = c(1 / 2, (1.007276 + 38.963158) / 2, "HK", "C0", 2),
    "M+H+Na" = c(1 / 2, (1.007276 + 22.989218) / 2, "HNa", "C0", 2),
    "M+ACN+2H" = c(1 / 2, (1.007276 + 41.026547) / 2, "C2H5N", "C0", 2),
    "M+2Na" = c(1 / 2, (2 * 22.989218) / 2, "Na2", "C0", 2),
    "M+2ACN+2H" = c(1 / 2, (2 * 41.026547 + 2 * 1.007276) / 2, "C4H8N2", "C0", 2),
    "M+3ACN+2H" = c(1 / 2, (3 * 41.026547 + 2 * 1.007276) / 2, "C6H11N3", "C0", 2),
    "M+H" = c(1, 1.007276, "H", "C0", 1),
    "M+NH4" = c(1, 18.033823, "NH4", "C0", 1),
    "M+Na" = c(1, 22.989218, "Na", "C0", 1),
    "M+CH3OH+H" = c(1, (32.026213 + 1.007276), "CH5O", "C0", 1),
    "M+K" = c(1, 38.963158, "K", "C0", 1),
    "M+ACN+H" = c(1, (41.026547 + 1.007276), "C2H4N", "C0", 1),
    "M+2Na-H" = c(1, (2 * 22.989218 - 1.007276), "Na2", "H1", 1),
    "M+iPrOH+H" = c(1, (60.058064 + 1.007276), "C3H9O", "C0", 1),
    "M+ACN+Na" = c(1, (41.026547 + 22.989218), "C2H3NNa", "C0", 1),
    "M+2K-H" = c(1, (2 * 38.963158 - 1.007276), "K2", "H1", 1),
    "M+DMSO+H" = c(1, (78.013944 + 1.007276), "C2H7OS", "C0", 1),
    "M+2ACN+H" = c(1, (2 * 41.026547 + 1.007276), "C4H7N2", "C0", 1),
    "2M+H" = c(2, 1.007276, "H", "C0", 1),
    "2M+NH4" = c(2, 18.033823, "NH4", "C0", 1),
    "2M+Na" = c(2, 22.989218, "Na", "C0", 1),
    "2M+K" = c(2, 38.963158, "K", "C0", 1),
    "2M+ACN+H" = c(2, (41.026547 + 1.007276), "C2H4N", "C0", 1),
    "2M+ACN+Na" = c(2, (41.026547 + 22.989218), "C2H3NNa", "C0", 1),
    "M" = c(1, 0, "C0", "C0", 0),
    "M-3H" = c(1 / 3, -1.007276, "C0", "H3", -3),
    "M-2H" = c(1 / 2, -1.007276, "C0", "H2", -2),
    "M-H" = c(1, -1.007276, "C0", "H", -1),
    "M+Na-2H" = c(1, (22.989218 - 2 * 1.007276), "Na", "H2", -1),
    "M+Cl" = c(1, 34.969402, "Cl", "C0", -1),
    "M+K-2H" = c(1, (38.963158 - 2 * 1.0072776), "K", "H2", -1),
    "M+FA-H" = c(1, 44.998201, "CHO2", "C0", -1),
    "M+HAc-H" = c(1, 59.013851, "C2H3O2", "C0", -1),
    "M+Br" = c(1, 78.918885, "Br", "C0", -1),
    "M+TFA-H" = c(1, 112.985586, "C2F3O2", "C0", -1),
    "2M-H" = c(2, -1.007276, "C0", "H", -1),
    "2M+FA-H" = c(2, 44.998201, "CHO2", "C0", -1),
    "2M+HAc-H" = c(2, 59.013851, "C2H3O2", "C0", -1),
    "3M-H" = c(3, -1.007276, "C0", "H", -1)
=======
    "[M+3H]3+" = c(1 / 3, (3 * 1.007276) / 3, "H3", "C0", 3),
    "[M+2H+Na]3+" = c(1 / 3, (2 * 1.007276 + 22.989218) / 3, "H2Na", "C0", 3),
    "[M+H+2Na]3+" = c(1 / 3, (1.007276 + 2 * 22.989218) / 3, "HNa2", "C0", 3),
    "[M+3Na]3+" = c(1 / 3, (3 * 22.989218) / 3, "Na3", "C0", 3),
    "[M+2H]2+" = c(1 / 2, (2 * 1.007276) / 2, "H2", "C0", 2),
    "[M+H+NH4]2+" = c(1 / 2, (1.007276 + 18.033823) / 2, "NH5", "C0", 2),
    "[M+H+K]2+" = c(1 / 2, (1.007276 + 38.963158) / 2, "HK", "C0", 2),
    "[M+H+Na]2+" = c(1 / 2, (1.007276 + 22.989218) / 2, "HNa", "C0", 2),
    "[M+ACN+2H]2+" = c(1 / 2, (1.007276 + 41.026547) / 2, "C2H5N", "C0", 2),
    "[M+2Na]2+" = c(1 / 2, (2 * 22.989218) / 2, "Na2", "C0", 2),
    "[M+2ACN+2H]2+" = c(1 / 2, (2 * 41.026547 + 2 * 1.007276) / 2, "C4H8N2", "C0", 2),
    "[M+3ACN+2H]2+" = c(1 / 2, (3 * 41.026547 + 2 * 1.007276) / 2, "C6H11N3", "C0", 2),
    "[M+H]+" = c(1, 1.007276, "H", "C0", 1),
    "[M+NH4]+" = c(1, 18.033823, "NH4", "C0", 1),
    "[M+Na]+" = c(1, 22.989218, "Na", "C0", 1),
    "[M+CH3OH+H]+" = c(1, (32.026213 + 1.007276), "CH5O", "C0", 1),
    "[M+K]+" = c(1, 38.963158, "K", "C0", 1),
    "[M+ACN+H]+" = c(1, (41.026547 + 1.007276), "C2H4N", "C0", 1),
    "[M+2Na-H]+" = c(1, (2 * 22.989218 - 1.007276), "Na2", "H1", 1),
    "[M+iPrOH+H]+" = c(1, (60.058064 + 1.007276), "C3H9O", "C0", 1),
    "[M+ACN+Na]+" = c(1, (41.026547 + 22.989218), "C2H3NNa", "C0", 1),
    "[M+2K-H]+" = c(1, (2 * 38.963158 - 1.007276), "K2", "H1", 1),
    "[M+DMSO+H]+" = c(1, (78.013944 + 1.007276), "C2H7OS", "C0", 1),
    "[M+2ACN+H]+" = c(1, (2 * 41.026547 + 1.007276), "C4H7N2", "C0", 1),
    "[2M+H]+" = c(2, 1.007276, "H", "C0", 1),
    "[2M+NH4]+" = c(2, 18.033823, "NH4", "C0", 1),
    "[2M+Na]+" = c(2, 22.989218, "Na", "C0", 1),
    "[2M+K]+" = c(2, 38.963158, "K", "C0", 1),
    "[2M+ACN+H]+" = c(2, (41.026547 + 1.007276), "C2H4N", "C0", 1),
    "[2M+ACN+Na]+" = c(2, (41.026547 + 22.989218), "C2H3NNa", "C0", 1),
    "[M]" = c(1, 0, "C0", "C0", 0),
    "[M-3H]3-" = c(1 / 3, -1.007276, "C0", "H3", -3),
    "[M-2H]2-" = c(1 / 2, -1.007276, "C0", "H2", -2),
    "[M-H]-" = c(1, -1.007276, "C0", "H", -1),
    "[M+Na-2H]-" = c(1, (22.989218 - 2 * 1.007276), "Na", "H2", -1),
    "[M+Cl]-" = c(1, 34.969402, "Cl", "C0", -1),
    "[M+K-2H]-" = c(1, (38.963158 - 2 * 1.0072776), "K", "H2", -1),
    "[M+FA-H]-" = c(1, 44.998201, "CHO2", "C0", -1),
    "[M+HAc-H]-" = c(1, 59.013851, "C2H3O2", "C0", -1),
    "[M+Br]-" = c(1, 78.918885, "Br", "C0", -1),
    "[M+TFA-H]-" = c(1, 112.985586, "C2F3O2", "C0", -1),
    "[2M-H]-" = c(2, -1.007276, "C0", "H", -1),
    "[2M+FA-H]-" = c(2, 44.998201, "CHO2", "C0", -1),
    "[2M+HAc-H]-" = c(2, 59.013851, "C2H3O2", "C0", -1),
    "[3M-H]-" = c(3, -1.007276, "C0", "H", -1)
>>>>>>> masstrixR_RaMoNA_merge
  )

  # return list
  return(adductCalc)
}

#' This function returns a list with the names of the adducts used to call the fitting formula. Adduts are taken from http://fiehnlab.ucdavis.edu/staff/kind/metabolomics/ms-adduct-calculator
#'
#' @return Returns a list with all adduct names
#' @examples
#' adductCalc <- getAdductCalc()
#' @export
getAdductNames <- function() {
  # return a vector with names of adducts
  return(c(
<<<<<<< HEAD
    "M+3H", "M+2H+Na", "M+H+2Na", "M+3Na", "M+2H", "M+H+NH4", "M+H+K", "M+H+Na", "M+ACN+2H", "M+2Na", "M+2ACN+2H",
    "M+3ACN+2H", "M+H", "M+NH4", "M+Na", "M+CH3OH+H", "M+K", "M+ACN+H", "M+2Na-H", "M+iPrOH+H", "M+ACN+Na", "M+2K-H",
    "M+DMSO+H", "M+2ACN+H", "2M+H", "2M+NH4", "2M+Na", "2M+K", "2M+ACN+H", "2M+ACN+Na",
    "M",
    "M-3H", "M-2H", "M-H", "M+Na-2H", "M+Cl", "M+K-2H", "M+FA-H", "M+HAc-H", "M+Br", "M+TFA-H", "2M-H", "2M+FA-H",
    "2M+HAc-H", "3M-H"
=======
    "[M+3H]3+", "[M+2H+Na]3+", "[M+H+2Na]3+", "[M+3Na]3+", "[M+2H]2+", "[M+H+NH4]2+", "[M+H+K]2+", "[M+H+Na]2+", "[M+ACN+2H]2+", "[M+2Na]2+", "[M+2ACN+2H]2+",
    "[M+3ACN+2H]2+", "[M+H]+", "[M+NH4]+", "[M+Na]+", "[M+CH3OH+H]+", "[M+K]+", "[M+ACN+H]+", "[M+2Na-H]+", "[M+iPrOH+H]+", "[M+ACN+Na]+", "[M+2K-H]+",
    "[M+DMSO+H]+", "[M+2ACN+H]+", "[2M+H]+", "[2M+NH4]+", "[2M+Na]+", "[2M+K]+", "[2M+ACN+H]+", "[2M+ACN+Na]+",
    "[M]",
    "[M-3H]3-", "[M-2H]2-", "[M-H]-", "[M+Na-2H]-", "[M+Cl]-", "[M+K-2H]-", "[M+FA-H]-", "[M+HAc-H]-", "[M+Br]-", "[M+TFA-H]-", "[2M-H]-", "[2M+FA-H]-",
    "[2M+HAc-H]-", "[3M-H]-"
>>>>>>> masstrixR_RaMoNA_merge
  ))
}

#' This function returns a list with typical neutral losses that can be observed as in-source fragments
#' @return Returns a named list with all mulitpliers and additive parts for positive and negative mode
#' @examples
#' neutralLossCalc <- getNeutralLossCalc()
#' @export
getNeutralLossCalc <- function() {

  # create list with multiplicative and additive part
  neutralLossCalc <- list(
    "M-H2O+H" = c(1, -17.003838, "H1", "H2O"),
    "M-H2O-H" = c(1, -19.018390, "C0", "H3O")
  )

  # return list
  return(neutralLossCalc)
}

#' This function returns a list with the names of the adducts used to call the fitting formula. Adduts are taken from http://fiehnlab.ucdavis.edu/staff/kind/metabolomics/ms-adduct-calculator
#'
#' @return Returns a list with all adduct names for the positive ionization mode
#' @examples
#' getAllPosModeAdducts()
#' @export
getAllPosModeAdducts <- function() {
  # return a vector with names of adducts
  return(c(
<<<<<<< HEAD
    "M+3H", "M+2H+Na", "M+H+2Na", "M+3Na", "M+2H", "M+H+NH4", "M+H+K", "M+H+Na", "M+ACN+2H", "M+2Na", "M+2ACN+2H",
    "M+3ACN+2H", "M+H", "M+NH4", "M+Na", "M+CH3OH+H", "M+K", "M+ACN+H", "M+2Na-H", "M+iPrOH+H", "M+ACN+Na", "M+2K-H",
    "M+DMSO+H", "M+2ACN+H", "2M+H", "2M+NH4", "2M+Na", "2M+K", "2M+ACN+H", "2M+ACN+Na"
=======
    "[M+3H]3+", "[M+2H+Na]3+", "[M+H+2Na]3+", "[M+3Na]3+", "[M+2H]2+", "[M+H+NH4]2+", "[M+H+K]2+", "[M+H+Na]2+", "[M+ACN+2H]2+", "[M+2Na]2+", "[M+2ACN+2H]2+",
    "[M+3ACN+2H]2+", "[M+H]+", "[M+NH4]+", "[M+Na]+", "[M+CH3OH+H]+", "[M+K]+", "[M+ACN+H]+", "[M+2Na-H]+", "[M+iPrOH+H]+", "[M+ACN+Na]+", "[M+2K-H]+",
    "[M+DMSO+H]+", "[M+2ACN+H]+", "[2M+H]+", "[2M+NH4]+", "[2M+Na]+", "[2M+K]+", "[2M+ACN+H]+", "[2M+ACN+Na]+"
>>>>>>> masstrixR_RaMoNA_merge
  ))
}
#' This function returns a list with the names of the adducts used to call the fitting formula. Adduts are taken from http://fiehnlab.ucdavis.edu/staff/kind/metabolomics/ms-adduct-calculator
#'
#' @return Returns a list with all adduct names for the negative ionization mode
#' @examples
#' getAllNegModeAdducts
#' @export
getAllNegModeAdducts <- function() {
  # return a vector with names of adducts
  return(c(
<<<<<<< HEAD
    "M-3H", "M-2H", "M-H", "M+Na-2H", "M+Cl", "M+K-2H", "M+FA-H", "M+HAc-H", "M+Br", "M+TFA-H", "2M-H", "2M+FA-H",
    "2M+HAc-H", "3M-H"
=======
    "[M-3H]3-", "[M-2H]2-", "[M-H]-", "[M+Na-2H]-", "[M+Cl]-", "[M+K-2H]-", "[M+FA-H]-", "[M+HAc-H]-", "[M+Br]-", "[M+TFA-H]-", "[2M-H]-", "[2M+FA-H]-",
    "[2M+HAc-H]-", "[3M-H]-"
>>>>>>> masstrixR_RaMoNA_merge
  ))
}
