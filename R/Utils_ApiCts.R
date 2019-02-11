#' This function uses the CTS Webservice from the Fiehnlab to convert IDs
#'
#' @param from format of original identifier (kegg, hmdb, inchikey)
#' @param to format of identifier to convert to
#' @param queryString query that shall be converted
#'
#' @export
getExternalDbIds <- function(from, to, queryString) {

  # sanity checks
  # check from and to
  if(!any(from %in% c("InChIKey", "Chemical Name", "KEGG" ,"Human Metabolome Database", "ChEBI", "BioCyc"))) {
    stop("from not correct, use InChIKey, Chemical Name, KEGG ,Human Metabolome Database, ChEBI or BioCyc")
  }

  if(!any(to %in% c("InChIKey", "Chemical Name", "KEGG" ,"Human Metabolome Database", "ChEBI", "BioCyc"))) {
    stop("to not correct, use InChIKey, Chemical Name, KEGG ,Human Metabolome Database, ChEBI or BioCyc")
  }

  # construct url for GET request
  baseUrl <- "http://cts.fiehnlab.ucdavis.edu/service/convert"
  queryUrl <- paste0(baseUrl, "/", from, "/", to, "/", queryString)

  result <- tryCatch({
    queryResult <- jsonlite::fromJSON(URLencode(queryUrl))

    return(queryResult$results)
  },
  error = function(cond) {
    return("error")
  },
  warning = function(cond) {
    return("error")
  },
  finally = {

  }
  )

  print(result)

  return(result)
}
