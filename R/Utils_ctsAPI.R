#' Function to access the CTS from the FiehnLab
#'
#'
.cts <- function() {

}



getExternalDbIds <- function(inchikey, db) {
  from <- "inchikey"
  to <- db
  queryString <- inchikey

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
