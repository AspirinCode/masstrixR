#' This function creates an SQLite database based on a given compound list. Minimum input is a metabolite name ($name), an exact mass ($exactMass), a formula ($formula) and an InChIKey ($inchikey)
#'
#' @param compoundList List of compounds that shall be added to DB
#' @param dbName A name for the database file
#' @param adductList Vector with adducts that shall be covered in the DB.
#' @return Returns the filename of the generated database
#' @examples
#' xxx
createDb <- function(compoundList, dbName, adductList) {

  #add some sanity checks here
  # TODO: does the data frame contain the right columns?

  #generate filename
  dbFileName <- paste0(dbName, ".sqlite")

  #get adduct calculation list
  adductCalc <- getAdductCalc()

  #create DF for upload
  dbupload <- data.frame()

  #iterate over adducts and create a query list
  for(adduct in adductList) {

    #generate list
    clipboard <- data.frame(metaboliteID = compoundList$id,
                            adductType = adduct,
                            adductMass = compoundList$exactmass * as.numeric(adductCalc[[adduct]][1]) + as.numeric(adductCalc[[adduct]][2]),
                            neutralMass = compoundList$exactmass,
                            neutralFormula = compoundList$formula,
                            ionFormula = unlist(lapply(compoundList$formula, calcAdductFormula, adduct = adduct)),
                            metaboliteName = stringr::str_c(compoundList$name, adduct, sep = " "),
                            inchkey = compoundList$inchikey,
                            inchi = compoundList$inchi,
                            smiles = compoundList$smiles)

    #add to list
    dbupload <- rbind.data.frame(dbupload, clipboard)

  }

  #upload to temp DB
  mydb <- DBI::dbConnect(RSQLite::SQLite(), dbFileName)
  DBI::dbWriteTable(mydb, "adducts", dbupload, overwrite = TRUE)

  #disconnect DB
  DBI::dbDisconnect(mydb)

  #return the filename
  return(dbFileName)
}
