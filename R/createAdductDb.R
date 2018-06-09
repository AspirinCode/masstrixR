#' This function creates an SQLite database based on a given compound list. Minimum input is a metabolite name ($name), an exact mass ($exactMass), a formula ($formula) and an InChIKey ($inchikey)
#'
#' @param compoundList List of compounds that shall be added to DB
#' @param dbName A name for the database file
#' @param adductList Vector with adducts that shall be covered in the DB.
#' @return Returns the filename of the generated database
#' @examples
#' xxx
createDb <- function(compoundList, dbName, adductList) {

  source("R\\adductCalc.R")

  #add some sanity checks here
  # TODO: does the data frame contain the right columns?

  #generate filename
  dbFileName <- paste0(dbFileName, ".sqlite")

  #get adduct calculation list
  adductCalc <- getAdductCalc()

  #create DF for upload
  dbupload <- data.frame()

  #iterate over adducts and create a query list
  for(adduct in adductList) {

    #generate list
    clipboard <- data.frame(metaboliteID = "",
                            adductType = adduct,
                            adductMass = (compoundList$exactMass * adductCalc[[adduct]][1]) + adductCalc[[adduct]][2],
                            neutralMass = compoundList$exactMass,
                            neutralFormula = compoundList$formula,
                            metaboliteName = str_c(compoundList$name, adduct, sep = " "),
                            inchkey = compoundList$inchikey)

    #add to list
    dbupload <- rbind.data.frame(dbupload, clipboard)

  }

  #upload to temp DB
  mydb <- DBI::dbConnect(RSQLite::SQLite(), dbFileName)
  DBI::dbWriteTable(mydb, "adducts", dbupload)

  #disconnect DB
  DBI::dbDisconnect(mydb)

  #return the filename
  return(dbFileName)
}
