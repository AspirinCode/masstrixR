#' This function performs a exact mass search in MoNA using an httr GET request and returns results a list. Data in this list can be access via the functions in parsingFunctions.R.
#'
#' @param exactMass exact mass to search for
#' @param error numeric value for maximum error
#' @param errorType type of error for search, either abs (absolute error in Da) or ppm (ppm error)
#'
#' @return Returns a List of lists with all information available (use parsingFunctions.R for simple data access)
searchByExactMassOffline <- function(exactMass, error = 0.01, errorType = "abs", dbCon) {

  #get absolute variables
  source("R/zzz.R", local = TRUE)

  #check error type and calculate lower and upper search boundaries
  if(errorType == "abs") {
    #get boundaries
    lower <- exactMass - error
    upper <- exactMass + error

  } else if(errorType == "ppm") {
    #get boundaries
    lower <- exactMass - (error / 10e6) * exactMass
    upper <- exactMass + (error / 10e6) * exactMass
  } else {
    stop("Wrong error type, use abs or ppm")
  }

  # make query in DB
  DBI::dbGetQuery(mydb, 'SELECT * FROM massSpec WHERE id = "000001"')
  DBI::dbGetQuery(mydb, 'SELECT * FROM spectra WHERE mz BETWEEN 100 and 200')




  # return a list with AnnotatedSpectrum2 objects
  return(TRUE)

}
