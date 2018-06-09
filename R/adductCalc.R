#' This function returns a list with adducts and the required multiplicative and additive part to calculate ion m/z from exact masses. Adduts are taken from http://fiehnlab.ucdavis.edu/staff/kind/metabolomics/ms-adduct-calculator
#'
#' @return Returns a named list with all mulitpliers and additive parts for positive and negative mode
#' @examples
#' adductCalc <- getAdductCalc()
getAdductCalc <- function() {

  #create list with multiplicative and additive part
  adductCalc <- list("M+3H" = c(1/3, (3 * 1.007276) / 3),
                     "M+2H+Na" = c(1/3, (2 * 1.007276 + 22.989218)/3),
                     "M+H+2Na" = c(1/3, (1.007276 + 2 * 22.989218)/ 3),
                     "M+3Na" = c(1/3, (3 * 22.989218) / 3),
                     "M+2H" = c(1/2, (2 * 1.007276) / 2),
                     "M+H+NH4" = c(1/2, (1.007276 + 18.033823) / 2),
                     "M+H+K" = c(1/2, (1.007276 + 38.963158) / 2),
                     "M+H+Na" = c(1/2, (1.007276 + 22.989218) / 2),
                     "M+ACN+2H" = c(1/2, (1.007276 + 41.026547) / 2),
                     "M+2Na" = c(1/2, (2 * 22.989218) / 2),
                     "M+2ACN+2H" = c(1/2, (2 * 41.026547 + 2 * 1.007276) / 2),
                     "M+3ACN+2H" = c(1/2, (3 * 41.026547 + 2 * 1.007276) / 2),
                     "M+H" = c(1, 1.007276),
                     "M+NH4" = c(1, 18.033823),
                     "M+Na" = c(1, 22.989218),
                     "M+CH3OH+H" = c(1, (32.026213 + 1.007276)),
                     "M+K" = c(1, 38.963158),
                     "M+ACN+H" = c(1, (41.026547 + 1.007276)),
                     "M+2Na-H" = c(1, (2 * 22.989218 - 1.007276)),
                     "M+iPrOH+H" = c(1, (60.058064 + 1.007276)),
                     "M+ACN+Na" = c(1, (41.026547 + 22.989218)),
                     "M+2K-H" = c(1, (2 * 38.963158 - 1.007276)),
                     "M+DMSO+H" = c(1, (78.013944 + 1.007276)),
                     "M+2ACN+H" = c(1, (2 * 41.026547 + 1.007276)),
                     "2M+H" = c(2, 1.007276),
                     "2M+NH4" = c(2, 18.033823),
                     "2M+Na" = c(2, 22.989218),
                     "2M+K" = c(2, 38.963158),
                     "2M+ACN+H" = c(2, (41.026547 + 1.007276)),
                     "2M+ACN+Na" = c(2, (41.026547 + 22.989218)),
                     "M" = c(1, 0),
                     "M-3H" = c(1/3, -1.007276),
                     "M-2H" = c(1/2, -1.007276),
                     "M-H" = c(1, -1.007276),
                     "M+Na-2H" = c(1, (22.989218 - 2 * 1.007276)),
                     "M+Cl" = c(1, 34.969402),
                     "M+K-2H" = c(1, (38.963158 - 2 * 1.0072776)),
                     "M+FA-H" = c(1, 44.998201),
                     "M+HAc-H" = c(1, 59.013851),
                     "M+Br" = c(1, 78.918885),
                     "M+TFA-H" = c(1, 112.985586),
                     "2M-H" = c(2, -1.007276),
                     "2M+FA-H" = c(2, 44.998201),
                     "2M+HAc-H" = c(2, 59.013851),
                     "3M-H" = c(3, -1.007276))


  #return list
  return(adductCalc)
}

getAdductNames <- function() {
  #return a vector with names of adducts
  return(c("M+3H", "M+2H+Na", "M+H+2Na", "M+3Na", "M+2H","M+H+NH4","M+H+K","M+H+Na","M+ACN+2H", "M+2Na", "M+2ACN+2H",
           "M+3ACN+2H", "M+H", "M+NH4", "M+Na", "M+CH3OH+H", "M+K", "M+ACN+H", "M+2Na-H", "M+iPrOH+H", "M+ACN+Na", "M+2K-H",
           "M+DMSO+H", "M+2ACN+H", "2M+H", "2M+NH4", "2M+Na", "2M+K", "2M+ACN+H", "2M+ACN+Na",
           "M",
           "M-3H", "M-2H", "M-H", "M+Na-2H", "M+Cl", "M+K-2H", "M+FA-H", "M+HAc-H", "M+Br", "M+TFA-H", "2M-H", "2M+FA-H",
           "2M+HAc-H","3M-H"))
}

#' This function returns a list with typical neutral losses that can be observed as in-source fragments
#' @return Returns a named list with all mulitpliers and additive parts for positive and negative mode
#' @examples
#' neutralLossCalc <- getNeutralLossCalc()
getNeutralLossCalc <- function() {

  #create list with multiplicative and additive part
  neutralLossCalc <- list("M-H2O+H" = c(1, -17.003838),
                          "M-H2O-H" = c(1, -19.018390))

  #return list
  return(neutralLossCalc)

}


