#' This function performs a formula search in MoNA using an httr GET request and returns results a list. Data in this list can be access via the functions in parsingFunctions.R.
#'
#' @param formula Formula to search for
#'
#' @return Returns a List of lists with all information available (use parsingFunctions.R for simple data access)
searchByFormulaOnline <- function(formula) {

  #get absolute variables
  source("R/zzz.R", local = TRUE)

  #GET request
  #perform GET request
  req <- httr::GET(MONASEARCH,
                   httr::add_headers(JSONCONTENT),
                   query = list(query = paste0('compound.metaData=q=\'name=="molecular formula" and value=="', formula, '"\'')))

  #return JSON results parsed into a list
  return(httr::content(req))
}

#' This function performs a exact mass search in MoNA using an httr GET request and returns results a list. Data in this list can be access via the functions in parsingFunctions.R.
#'
#' @param exactMass exact mass to search for
#' @param error numeric value for maximum error
#' @param errorType type of error for search, either abs (absolute error in Da) or ppm (ppm error)
#'
#' @return Returns a List of lists with all information available (use parsingFunctions.R for simple data access)
searchByExactMassOnline <- function(exactMass, error = 0.01, errorType = "abs") {

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

  #GET request
  #return as json
  req <- httr::GET(MONASEARCH,
                   httr::add_headers(JSONCONTENT),
                   query = list(query = paste0('compound.metaData=q=\'name==\"total exact mass\" and value >= ', lower, ' and value <= ', upper, '\'')))


  #return JSON results parsed into a list
  return(httr::content(req))

}

#' This function reconstructs a Ms2 spectrum for all fragments with a correlation value bigger than a set cut-off
#'
#' @param exactMass exact mass to search for
#' @param error numeric value for maximum error
#' @param errorType type of error for search, either abs (absolute error in Da) or ppm (ppm error)
#' @param precursorType Indicates possible adduct ions to search for (currently supported: [M+H]+, [M-H]-)
#'
#' @return Returns a List of lists with all information available (use parsingFunctions.R for simple data access)
#search by precursor and adduct type
searchByPrecursorOnline <- function(precursor, error = 0.01, errorType = "abs", precursorType = "[M+H]+") {

  #get absolute variables
  source("R/zzz.R", local = TRUE)

  #sanity check for precursorType
  if(!(precursorType %in% c("[M+H]+", "[M-H]-"))) {
    stop("Unknown precursor ion type supplied")
  }

  #check error type and calculate lower and upper search boundaries
  if(errorType == "abs") {
    #get boundaries
    lower <- precursor - error
    upper <- precursor + error

  } else if(errorType == "ppm") {
    #get boundaries
    lower <- precursor - (error / 10e6) * precursor
    upper <- precursor + (error / 10e6) * precursor
  } else {
    stop("Wrong error type, use abs or ppm")
  }

  #GET request
  #return as json
  req <- httr::GET(MONASEARCH,
                   httr::add_headers(JSONCONTENT),
                   query = list(query = paste0('metaData=q=\'name==\"precursor m/z\" and value >= ', lower, ' and value <= ', upper, '\' and (metaData=q=\'name==\"precursor type\" and value == \"', precursorType, '\"\')')))

  #return JSON results parsed into a list
  return(httr::content(req))

}
