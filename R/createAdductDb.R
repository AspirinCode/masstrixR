#' This function creates an SQLite database based on a given compound list. Minimum input is a metabolite name ($name), an exact mass ($exactMass), a formula ($formula) and an InChIKey ($inchikey)
#'
#' @param compoundList List of compounds that shall be added to DB
#' @param dbName A name for the database file
#' @param adductList Vector with adducts that shall be covered in the DB.
#' @return Returns the filename of the generated database
#' @examples
#' xxx
createDb <- function(compoundList, dbName, adductList, rt = TRUE) {

  #add some sanity checks here
  # correct column headers in compound list?
  if(!all(c("id", "smiles", "inchi", "inchikey", "formula", "name", "exactmass") %in% colnames(compoundList))) {
    stop("One or more column header does not match the required headers: id, smiles, inchi, inchikey, formula, name, exactmass.")
  }

  if(rt = TRUE & !c("rt") %in% colnames(compoundList)) {
    stop("selected rt option but no rt column defined in compound list")
  }

  # check if adduct names are correct
  if(!all(adductList %in% getAdductNames())) {
    stop("One or more adduct name does not match the required adduct names.")
  }

  #generate filename
  dbFileName <- paste0(dbName, ".sqlite")

  #get adduct calculation list
  adductCalc <- getAdductCalc()

  #create DF for upload
  dbupload <- data.frame()

  #iterate over adducts and create a query list
  for(adduct in adductList) {

    # add rt if TRUE
    if(rt = TRUE) {
      ritmes <- compoundList$rt
    } else {
      rtimes <- rep(0, nrow(compoundList))
    }

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
                            smiles = compoundList$smiles,
                            rt = rtimes,
                            ccs = NA)

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

#' This function creates an SQLite database based on a given compound list. Minimum input is a metabolite name ($name), an exact mass ($exactMass), a formula ($formula) and an InChIKey ($inchikey)
#'
#' @param compoundList List of compounds that shall be added to DB
#' @param dbName A name for the database file
#' @param adductList Vector with adducts that shall be covered in the DB.
#' @return Returns the filename of the generated database
#' @examples
#' xxx
createDbMiscAdducts <- function(compoundList, dbName, rt = TRUE) {

  # correct column headers in compound list?
  if(!all(c("id", "smiles", "inchi", "inchikey", "formula", "name", "exactmass", "adducts") %in% colnames(compoundList))) {
    stop("One or more column header does not match the required headers: id, smiles, inchi, inchikey, formula, name, exactmass.")
  }

  if(rt = TRUE & !c("rt") %in% colnames(compoundList)) {
    stop("selected rt option but no rt column defined in compound list")
  }

  #generate filename
  dbFileName <- paste0(dbName, ".sqlite")

  #get adduct calculation list
  adductCalc <- getAdductCalc()

  # create DF for upload
  dbupload <- data.frame()

  #iterate over compound list and add to dbupload
  for(i in 1:nrow(compoundList)) {
    if(grepl(",", compoundList$adducts[i])) {
      adductList <- stringr::str_split(compoundList$adducts[i], ",")[[1]]
    } else {
      adductList <- compoundList$adducts[i]
    }

    #iterate over adducts and create a query list
    for(adduct in adductList) {

      # check if adduct names are correct
      if(!all(adductList %in% getAdductNames())) {
        stop("One or more adduct name does not match the required adduct names.")
      }

      #generate list
      clipboard <- data.frame(metaboliteID = compoundList$id[i],
                              adductType = adduct,
                              adductMass = compoundList$exactmass[i] * as.numeric(adductCalc[[adduct]][1]) + as.numeric(adductCalc[[adduct]][2]),
                              neutralMass = compoundList$exactmass[i],
                              neutralFormula = compoundList$formula[i],
                              ionFormula = calcAdductFormula(compoundList$formula[i], adduct),
                              metaboliteName = stringr::str_c(compoundList$name[i], adduct, sep = " "),
                              inchkey = compoundList$inchikey[i],
                              inchi = compoundList$inchi[i],
                              smiles = compoundList$smiles[i],
                              rt = compoundList$rt[i],
                              ccs = NA)

      #add to list
      dbupload <- rbind.data.frame(dbupload, clipboard)

    }
  }

  #upload to temp DB
  mydb <- DBI::dbConnect(RSQLite::SQLite(), dbFileName)
  DBI::dbWriteTable(mydb, "adducts", dbupload, overwrite = TRUE)

  #disconnect DB
  DBI::dbDisconnect(mydb)

  #return the filename
  return(dbFileName)

}

# function to create DB with user defined input
createDbUserList <- function(compoundList, dbName, rt = FALSE, ccs = FALSE) {

  # correct column headers in compound list?
  if(!all(c("metaboliteID", "adductType", "adductMass", "neutralMass", "neutralFormula", "ionFormula", "metaboliteName", "inchikey", "inchi", "smiles") %in% colnames(compoundList))) {
    stop("One or more column header does not match the required headers: id, smiles, inchi, inchikey, formula, name, exactmass.")
  }

  if(rt = TRUE & !c("rt") %in% colnames(compoundList)) {
    stop("rt option selected, but no rt column found")
  }

  if(ccs = TRUE & !c("ccs") %in% colnames(compoundList)) {
    stop("ccs option selected, but no ccs column found")
  }

  # add empty column
  if(rt = FALSE) {
    compoundList$rt <- NA
  }

  if(ccs = FALSE) {
    compoundList$ccs <- NA
  }

  #generate filename
  dbFileName <- paste0(dbName, ".sqlite")

  # list for upload
  dbupload <- compoundList

  #upload to temp DB
  mydb <- DBI::dbConnect(RSQLite::SQLite(), dbFileName)
  DBI::dbWriteTable(mydb, "adducts", dbupload, overwrite = TRUE)

  #disconnect DB
  DBI::dbDisconnect(mydb)

  #return the filename
  return(dbFileName)

}
