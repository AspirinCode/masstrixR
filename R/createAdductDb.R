#' This function creates an SQLite database based on a given compound list. Minimum input is a metabolite name ($name), an exact mass ($exactMass), a formula ($formula) and an InChIKey ($inchikey)
#'
#' @param compoundList List of compounds that shall be added to DB
#' @param dbName A name for the database file
#' @param adductList Vector with adducts that shall be covered in the DB.
#' @return Returns the filename of the generated database
#' @examples
#' xxx
createDb <- function(compoundList, dbName, adductList = NA, rt = FALSE, ccs = FALSE, ccsType = "DT", extId = FALSE) {

  ################################################################################################################################################
  # Sanity Checks on input
  ################################################################################################################################################
  # correct column headers in compound list?
  if(!all(c("id", "smiles", "inchi", "inchikey", "formula", "name", "exactmass") %in% colnames(compoundList))) {
    stop("One or more column header does not match the required headers: id, smiles, inchi, inchikey, formula, name, exactmass.")
  }

  # information on adducts available
  if(!(c("adducts") %in% colnames(compoundList)) && is.na(adductList)) {
    stop("adducts have to be defined!!!")
  }

  # is there an RT column, if RT is true
  if(rt == TRUE & !c("rt") %in% colnames(compoundList)) {
    stop("selected rt option but no rt column defined in compound list")
  }

  # is there an CCS column, if CCS is true
  if(ccs == TRUE & !c("ccs") %in% colnames(compoundList)) {
    stop("selected ccs option but no ccs column defined in compound list")
  }

  # if CCS is true then individual adducts for each metabolite are required
  if(ccs == TRUE & !c("adducts") %in% colnames(compoundList) & all(grepl(",", compoundList$adducts))) {
    stop("ccs option requires individual adduct annotation")
  }

  if(ccs == TRUE & !ccsType %in% c("DT", "TWIMS", "TIMS", "predicted")) {
    stop("Unknown CCS Type")
  }

  # set standard value for fetching
  #fetchIds <- FALSE

  # are external ids part of the compound list
  if(extId == TRUE & !any(c("kegg", "hmdb", "chebi") %in% colnames(compoundList))) {
    message("external ID option was selected, but no KEGG, HMDB or ChEBI ID supplied. IDs will be fetched from CTS.")
    fetchIds <- TRUE
  } else {
    fetchIds <- FALSE
  }

  ################################################################################################################################################
  # get/generate required values
  ################################################################################################################################################
  #generate filename
  dbFileName <- paste0(dbName, ".sqlite")

  # generate data for upload based on input
  dbupload <- .generateSQLiteInput(compoundList, adductList = adductList, rt = rt, ccs = ccs, fetchIds = fetchIds)

  ################################################################################################################################################
  # Generate SQLite DB
  ################################################################################################################################################
  #upload to DB
  mydb <- DBI::dbConnect(RSQLite::SQLite(), dbFileName)
  DBI::dbWriteTable(mydb, "adducts", dbupload, overwrite = TRUE)

  #disconnect DB
  DBI::dbDisconnect(mydb)

  ################################################################################################################################################
  # Return path to SQLite file
  ################################################################################################################################################
  #return the filename
  return(dbFileName)
}






###### testing function for upload
createDbUserDefined <- function(compoundList, dbName, rt = FALSE, ccs = FALSE, ccsType = "DT", extId = FALSE) {

  ################################################################################################################################################
  # Sanity Checks on input
  ################################################################################################################################################
  # correct column headers in compound list?
  if(!all(c("metaboliteID", "adductType", "adductMass", "neutralMass", "neutralFormula", "ionFormula", "metaboliteName", "inchikey", "inchi", "smiles") %in% colnames(compoundList))) {
    stop("One or more column header does not match the required headers: metaboliteID, adductType, adductMass, neutralMass, neutralFormula, ionFormula, metaboliteName, inchikey, inchi, smiles")
  }

  # is there an RT column, if RT is true
  if(rt == TRUE & !c("rt") %in% colnames(compoundList)) {
    stop("selected rt option but no rt column defined in compound list")
  }

  # is there an CCS column, if CCS is true
  if(ccs == TRUE & !c("ccs") %in% colnames(compoundList)) {
    stop("selected ccs option but no ccs column defined in compound list")
  }

  if(ccs == TRUE & !ccsType %in% c("DT", "TWIMS", "TIMS", "predicted")) {
    stop("Unknown CCS Type")
  }

  # check if external id columns exist
  if(extId == TRUE & !all(c("kegg", "hmdb", "chebi") %in% colnames(compoundList))) {
    stop("missing kegg, hmdb and chebi column")
  }


  ################################################################################################################################################
  # add missing columns
  ################################################################################################################################################
  if(rt == FALSE) {
    compoundList$rt <- NA
  }

  if(ccs == FALSE) {
    compoundList$ccs <- NA
  }

  if(extId == FALSE) {
    compoundList$kegg <- NA
    compoundList$hmdb <- NA
    compoundList$chebi <- NA
  }

  ################################################################################################################################################
  # prepare DF for upload
  ################################################################################################################################################
  dbupload <- compoundList

  ################################################################################################################################################
  # Generate SQLite DB
  ################################################################################################################################################
  #upload to DB
  mydb <- DBI::dbConnect(RSQLite::SQLite(), dbFileName)
  DBI::dbWriteTable(mydb, "adducts", dbupload, overwrite = TRUE)

  #disconnect DB
  DBI::dbDisconnect(mydb)

  ################################################################################################################################################
  # Return path to SQLite file
  ################################################################################################################################################
  #return the filename
  return(dbFileName)

}
