#' This function annotates masses in a peak list with possible metabolites from a selected DB. Used adducts are the adducts covered by the respective database
#'
#' @param peakList Data frame containing the peaks of interest, minimum column is called m.z ($m.z)
#' @param dbFileName A .sqlite file containing the DB of interest
#' @param mode Defines if the database shall be used from the disk ("onDisk") or in memory ("inMemory"). The second option enhances performances with very large peaks lists and DBs.
#' @param mzTol Maximum allowed tolerance in Da (ppm will come soon)
#' @return Returns the filename of the generated database
#' @examples
#' xxx
#search in a pre-defined database
mzSearch <-function(peakList, dbFileName, mode = "onDisk", mzTol = 0.005, mzTolType = "abs",
                    rt = FALSE, rtTol = 0.5, rtTolType = "abs",
                    ccs = FALSE, ccsTol = 1, ccsTolType = "rel") {

  #some sanity checks here
  #TODO: add checks on data frame, DB etc...
  if(!all(c("m.z") %in% colnames(peakList)))

  #depended on the chosen mode different connections are required
  if(mode == "onDisk") {

    #connect to DB
    mydb <- DBI::dbConnect(RSQLite::SQLite(), dbFileName)

  } else if(mode == "inMemory") {

    #make connection to DB and copy to memory DB
    tempDb <- DBI::dbConnect(RSQLite::SQLite(), dbFileName)
    mydb <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")

    #copy database
    RSQLite::sqliteCopyDatabase(tempDb, mydb)

    #disconnect not required DB
    DBI::dbDisconnect(tempDb)

  } else {
    stop("Unknown mode")
  }

  #perform MS1 annotation
  ms1annotation <- data.frame()

  for(i in 1:nrow(peakList)) {

    ######################################################
    # m/z boundaries
    ######################################################
    # calc lower and upper mz boundary
    if(mzTolType == "abs") {
      lowerMz <- peakList$m.z[i] - mzTol
      upperMz <- peakList$m.z[i] + mzTol
    } else if(mzTolType == "ppm") {
      #lowerMz <- peakList$m.z[i] - mzTol
      #upperMz <- peakList$m.z[i] + mzTol
    } else {
      stop("Unknown tolType defined!")
    }

    ######################################################
    # RT boundaries
    ######################################################
    # calc lower and upper RT boundary
    # TODO


    ######################################################
    # CCS boundaries
    ######################################################
    # calc lower and upper RT boundary
    # TODO


    ######################################################
    # check database if columns required exist
    ######################################################
    # check DB
    # TODO: write function to check if RT and CCS columns exist.


    #execute query
    resultSet <- DBI::dbSendQuery(mydb, 'SELECT * FROM adducts WHERE adductMass BETWEEN :lowerMz AND :upperMz')
    DBI::dbBind(resultSet, param = list(lowerMz = lowerMz,
                                        upperMz = upperMz))

    #fetch result set into a dataframe and check if an annotation was found
    annotation <- DBI::dbFetch(resultSet)
    if(nrow(annotation) > 0) {
      ms1annotation <- rbind.data.frame(ms1annotation, cbind.data.frame(peakList[i,], annotation))
      #print(cbind.data.frame(peakList[i,], annotation))
    }

    #clear result
    DBI::dbClearResult(resultSet)
  }

  #disconnect DB
  DBI::dbDisconnect(mydb)

  #calculate error
  ms1annotation$absError <- ms1annotation$adductMass - ms1annotation$m.z

  #return results
  return(ms1annotation)

}


#' This function annotates masses in a peak list with possible metabolites from a selected DB. Used adducts are the adducts covered by the respective database
#'
#' @param peakList Data frame containing the peaks of interest, minimum column is called m.z ($m.z)
#' @param dbFileName A .sqlite file containing the DB of interest
#' @param mode Defines if the database shall be used from the disk ("onDisk") or in memory ("inMemory"). The second option enhances performances with very large peaks lists and DBs.
#' @param mzTol Maximum allowed tolerance in Da (ppm will come soon)
#' @return Returns the filename of the generated database
#' @examples
#' xxx
#search in a pre-defined database
mzLookUp <- function(peakList, compoundList, adductList, mzTol = 0.005, tolType = "abs") {

  #add some sanity checks here
  # correct column headers in compound list?
  if(!all(c("id", "smiles", "inchi", "inchikey", "formula", "name", "exactmass") %in% colnames(compoundList))) {
    stop("One or more column header does not match the required headers: id, smiles, inchi, inchikey, formula, name, exactmass.")
  }

  # check if adduct names are correct
  if(!all(adductList %in% getAdductNames())) {
    stop("One or more adduct name does not match the required adduct names.")
  }

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
  mydb <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  DBI::dbWriteTable(mydb, "adducts", dbupload)

  #perform MS1 annotation
  ms1annotation <- data.frame()

  for(i in 1:nrow(peakList)) {
    lower <- peakList$m.z[i] - mzTol
    upper <- peakList$m.z[i] + mzTol

    #execute query
    resultSet <- DBI::dbSendQuery(mydb, 'SELECT * FROM adducts WHERE adductMass BETWEEN :lower AND :upper')
    DBI::dbBind(resultSet, param = list(lower = lower,
                                        upper = upper))

    #fetch result set into a dataframe and check if an annotation was found
    annotation <- DBI::dbFetch(resultSet)
    if(nrow(annotation) > 0) {
      ms1annotation <- rbind.data.frame(ms1annotation, cbind.data.frame(peakList[i,], annotation))
      #print(cbind.data.frame(peakList[i,], annotation))
    }

    #clear result
    DBI::dbClearResult(resultSet)
  }

  #disconnect DB
  DBI::dbDisconnect(mydb)

  #calculate error
  ms1annotation$absError <- ms1annotation$adductMass - ms1annotation$m.z

  #return results
  return(ms1annotation)

}

#' This function annotates masses in a peak list with possible metabolites from a selected DB. Used adducts are the adducts covered by the respective database
#'
#' @param peakList Data frame containing the peaks of interest, minimum column is called m.z ($m.z)
#' @param dbFileName A .sqlite file containing the DB of interest
#' @param mode Defines if the database shall be used from the disk ("onDisk") or in memory ("inMemory"). The second option enhances performances with very large peaks lists and DBs.
#' @param mzTol Maximum allowed tolerance in Da (ppm will come soon)
#' @return Returns the filename of the generated database
#' @examples
#' xxx
#search in a pre-defined database
mzLookUpMiscAdducts <- function(peakList, compoundList, mzTol = 0.005, tolType = "abs") {

  #add some sanity checks here
  # correct column headers in compound list?
  if(!all(c("id", "smiles", "inchi", "inchikey", "formula", "name", "exactmass", "adducts") %in% colnames(compoundList))) {
    stop("One or more column header does not match the required headers: id, smiles, inchi, inchikey, formula, name, exactmass.")
  }

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
                              smiles = compoundList$smiles[i])

      #add to list
      dbupload <- rbind.data.frame(dbupload, clipboard)

    }
  }

  #upload to temp DB
  mydb <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  DBI::dbWriteTable(mydb, "adducts", dbupload)

  #perform MS1 annotation
  ms1annotation <- data.frame()

  for(i in 1:nrow(peakList)) {
    lower <- peakList$m.z[i] - mzTol
    upper <- peakList$m.z[i] + mzTol

    #execute query
    resultSet <- DBI::dbSendQuery(mydb, 'SELECT * FROM adducts WHERE adductMass BETWEEN :lower AND :upper')
    DBI::dbBind(resultSet, param = list(lower = lower,
                                        upper = upper))

    #fetch result set into a dataframe and check if an annotation was found
    annotation <- DBI::dbFetch(resultSet)
    if(nrow(annotation) > 0) {
      ms1annotation <- rbind.data.frame(ms1annotation, cbind.data.frame(peakList[i,], annotation))
      #print(cbind.data.frame(peakList[i,], annotation))
    }

    #clear result
    DBI::dbClearResult(resultSet)
  }

  #disconnect DB
  DBI::dbDisconnect(mydb)

  #calculate error
  ms1annotation$absError <- ms1annotation$adductMass - ms1annotation$m.z

  #return results
  return(ms1annotation)

}
