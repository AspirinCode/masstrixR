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
mzSearch <-function(peakList, dbFileName, mode = "onDisk", mzTol = 0.005, tolType = "abs") {

  #some sanity checks here
  #TODO: add checks on data frame, DB etc...

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
      print(cbind.data.frame(peakList[i,], annotation))
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

  source("R\\adductCalc.R")
  source("R\\formulaUtils.R")

  #add some sanity checks here
  # TODO: does the data frame contain the right columns?

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
                            inchkey = compoundList$inchikey)

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
      print(cbind.data.frame(peakList[i,], annotation))
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
