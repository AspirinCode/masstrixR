#' Performing m/z search
#'
#' This function performs annotation of m/z values with putative metabolites using a previously created SQLite DB. Measured m/z values are compared against theoretical values within a defined error range. This error can be either absolute (in Da) or relative (in ppm).
#' If the database contains RT or CCS values results can be a additionally filtered.
#'
#' @param peakList Data frame containing the peaks of interest, minimum column is called m.z ($m.z)
#' @param dbFileName A .sqlite file containing the DB of interest
#' @param mode Defines if the database shall be used from the disk ("onDisk") or in memory ("inMemory"). The second option enhances performances with very large peaks lists and DBs.
#' @param mzTol Maximum allowed tolerance in Da or ppm
#' @param mzTolType Defines the error type for m/z search, "abs" is used for absolute mass error, "ppm" for relative error
#' @param rt Boolean value indicating of RT filtering shall be performed
#' @return Returns the filename of the generated database
#' @examples
#' mzSearch()
#' @export
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
  # rework
  # TODO check for NULL in columns!!!
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
      lowerMz <- peakList$m.z[i] - mzTol / 1e6 * peakList$m.z[i]
      upperMz <- peakList$m.z[i] + mzTol / 1e6 * peakList$m.z[i]
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

#' This function annotates masses in a peak list with possible metabolites from a selected DB. Used adducts are the adducts covered by the respective database
#'
#' @param peakList Data frame containing the peaks of interest, minimum column is called m.z ($m.z)
#' @param dbFileName A .sqlite file containing the DB of interest
#' @param mode Defines if the database shall be used from the disk ("onDisk") or in memory ("inMemory"). The second option enhances performances with very large peaks lists and DBs.
#' @param mzTol Maximum allowed tolerance in Da (ppm will come soon)
#' @return Returns the filename of the generated database
#' @examples
#' mzLookup()
#'
#' @export
mzLookUp <- function(peakList, compoundList, mzTol = 0.005, mzTolType = "abs",
                     rt = FALSE, rtTol = 0.5, rtTolType = "abs",
                     ccs = FALSE, ccsTol = 1, ccsTolType = "rel") {

  ################################################################################################################################################
  # Sanity Checks on input
  ################################################################################################################################################
  # correct column headers in compound list?
  if(!all(c("metaboliteID", "adductType", "adductMass", "neutralMass", "neutralFormula", "ionFormula", "metaboliteName", "inchikey", "inchi", "smiles", "rt", "ccs", "kegg", "hmdb", "chebi") %in% colnames(compoundList))) {
    stop("One or more column header does not match the required headers: id, smiles, inchi, inchikey, formula, name, exactmass.")
  }

  # # information on adducts available
  # if(!(c("adducts") %in% colnames(compoundList)) && is.na(adductList)) {
  #   stop("adducts have to be defined!!!")
  # }

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

  #  # generate data for upload based on input
  dbupload <- compoundList

  #upload to temp DB
  mydb <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  DBI::dbWriteTable(mydb, "adducts", dbupload)

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
