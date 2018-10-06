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
mzSearch <-function(peakList, dbFileName, mode = "onDisk",
                    mzTol = 0.005, mzTolType = "abs",
                    rt = FALSE, rtTol = 0.5, rtTolType = "abs",
                    ccs = FALSE, ccsTol = 1, ccsTolType = "rel") {

  #some sanity checks here
  #TODO: add checks on data frame, DB etc...
  if(!all(c("m.z") %in% colnames(peakList))) {
    stop("no m/z column defined")
  }

  # check if RT column is supplied
  if(rt == TRUE & !all(c("RT") %in% colnames(peakList))) {
    stop("RT searching is TRUE, but no RT column found in peaklist")
  }

  # check if CCS column is supplied
  if(ccs == TRUE & !all(c("CCS") %in% colnames(peakList))) {
    stop("CCS searching is TRUE, but no CCS column found in peaklist")
  }

  #depended on the chosen mode different connections are required
  if(mode == "onDisk") {

    #connect to DB
    mydb <- DBI::dbConnect(RSQLite::SQLite(), dbFileName)

  } else if(mode == "inMemory") {

    #make connection to DB and copy to memory DB
    tempdb <- DBI::dbConnect(RSQLite::SQLite(), dbFileName)
    mydb <-DBI::dbConnect(RSQLite::SQLite(), ":memory:")

    #copy database
    RSQLite::sqliteCopyDatabase(tempdb, mydb)

    #disconnect not required DB
    DBI::dbDisconnect(tempdb)

  } else {
    stop("Unknown mode")
  }

  ######################################################
  # check database if columns required exist
  ######################################################
  # check DB if correct columns exist
  mydbColumns <- DBI::dbFetch(DBI::dbSendQuery(mydb, "PRAGMA table_info(adducts)"))

  # check for RT
  if(rt == TRUE & !any(grepl("rt", mydbColumns$name))) {
    stop("RT option selected, but selected DB does not contain RT data")
  }

  # check for CCS
  if(ccs == TRUE & !any(grepl("ccs", mydbColumns$name))) {
    stop("CCS option selected, but selected DB does not contain CCS data")
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
      stop("Unknown mzTolType defined!")
    }

    ######################################################
    # RT boundaries
    ######################################################
    # calc lower and upper RT boundary
    if(rt == TRUE & rtTolType == "abs") {
      lowerRt <- peakList$RT[i] - rtTol
      upperRt <- peakList$RT[i] + rtTol
    } else if(rt == TRUE & mzTolType == "rel") {
      lowerRt <- peakList$RT[i] - peakList$RT * rtTolType
      upperRt <- peakList$RT[i] + peakList$RT * rtTolType
    } else if(rt == TRUE) {
      stop("unknown rtTolTyp defined!")
    }

    ######################################################
    # CCS boundaries
    ######################################################
    # calc lower and upper RT boundary
    if(ccs == TRUE & ccsTolType == "abs") {
      lowerCcs <- peakList$CCS[i] - ccsTol
      upperCcs <- peakList$CCS[i] + ccsTol
    } else if(ccs == TRUE & ccsTolType == "rel") {
      lowerCcs <- peakList$CCS[i] - peakList$CCS * ccsTolType
      upperCcs <- peakList$CCS[i] + peakList$CCS * ccsTolType
    } else if(ccs == TRUE) {
      stop("unknown rtTolTyp defined!")
    }

    # generate query
    query <- 'SELECT * FROM adducts WHERE (adductMass BETWEEN :lowerMz AND :upperMz)'
    param <- list(lowerMz = lowerMz,
                  upperMz = upperMz)

    # add part for RT search
    if(rt == TRUE) {
      query <- paste0(query, 'AND (rt BETWEEN :lowerRt AND :upperRt)')
      param <- c(param, list(lowerRt = lowerRt,
                             upperRt = upperRt))
    }

    # add part for CCS search
    if(ccs == TRUE) {
      query <- paste0(query, 'AND (ccs BETWEEN :lowerCcs AND :upperCcs)')
      param <- c(param, list(lowerCcs = lowerCcs,
                             upperCcs = upperCcs))
    }

    #execute query
    resultSet <- DBI::dbSendQuery(mydb, query)
    DBI::dbBind(resultSet, param = param)

    #fetch result set into a dataframe and check if an annotation was found
    annotation <- DBI::dbFetch(resultSet)
    if(nrow(annotation) > 0) {
      ms1annotation <- rbind.data.frame(ms1annotation, cbind.data.frame(peakList[i,], annotation))
    }

    #clear result
    DBI::dbClearResult(resultSet)
  }

  # disconnect DB
  DBI::dbDisconnect(mydb)

  #calculate absolute and relative m/z error
  ms1annotation$mzAbsError <- ms1annotation$m.z - ms1annotation$adductMass
  ms1annotation$mzRelError <- (ms1annotation$m.z - ms1annotation$adductMass) / ms1annotation$adductMass * 1000000

  # calculate RT error
  if(rt == TRUE) {
    ms1annotation$rtAbsError <- ms1annotation$RT - ms1annotation$rt
    ms1annotation$rtRelError <- (ms1annotation$RT - ms1annotation$rt) / ms1annotation$rt
  }

  # calculate RT error
  if(ccs == TRUE) {
    ms1annotation$ccsAbsError <- ms1annotation$CCS - ms1annotation$ccs
    ms1annotation$ccsRelError <- (ms1annotation$CCS - ms1annotation$ccs) / ms1annotation$ccs
  }

  #return results
  return(ms1annotation)

}
