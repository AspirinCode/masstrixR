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
mzLookUp <- function(peakList, compoundList, adductList, mzTol = 0.005, tolType = "abs",
                     rt = FALSE, rtTol = 0.5, rtTolType = "abs",
                     ccs = FALSE, ccsTol = 1, ccsTolType = "rel") {

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

  #  # generate data for upload based on input
  dbupload <- .generateSQLiteInput(compoundList, adductList = adductList, rt = rt, ccs = ccs)

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
