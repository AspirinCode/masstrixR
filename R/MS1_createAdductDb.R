#' Creation of SQLiteDB
#'
#' masstrixR uses SQLite DBs as backend for the annotation process. In order to perform annotation of m/z values with putative metabolites a DB is required.
#' This function creates a SQLite database using a validated compound list containing all required information. Such a list can be created using the \code{prepareCompoundList()} function.
#'
#' @param compoundList List of compounds that shall be added to DB
#' @param dbName A name for the database file
#'
#' @examples
#' xxx
#'
#' @author Michael Witting, \email{michael.witting@@helmholtz-muenchen.de}
#'
#' @seealso \code{\link{validateCompoundList}}
#' @seealso \code{\link{prepareCompoundList}}
#' @export
createDb <- function(compoundList, dbName) {

  ################################################################################################################################################
  # Sanity Checks on input
  ################################################################################################################################################
  # correct column headers in compound list?
  if(!all(c("metaboliteID", "adductType", "adductMass", "neutralMass", "neutralFormula", "ionFormula", "metaboliteName", "inchikey", "inchi", "smiles", "rt", "ccs", "kegg", "hmdb", "chebi") %in% colnames(compoundList))) {
    stop("One or more column header does not match the required headers: id, smiles, inchi, inchikey, formula, name, exactmass.")
  }

  # check correct format of each column
  # TODO

  ################################################################################################################################################
  # get/generate required values
  ################################################################################################################################################
  #generate filename
  dbFileName <- paste0(dbName, ".sqlite")

  # generate data for upload based on input
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

#' Validation of compound list
#'
#' Compound list that are used to generate SQLite DBs require a speficic format prior to generation of the actual DB.
#'
#' This function validates if all required columns are in place and of the right format. Specific columns are added, if required (e.g. KEGG Ids)
#'
#' @param compoundList List of compounds that shall be added to DB
#' @return Returns the filename of the generated database
#' @examples
#' xxx
#'
#' @author Michael Witting, \email{michael.witting@@helmholtz-muenchen.de}
#'
#' @seealso \code{\link{createDb}}
#' @seealso \code{\link{prepareCompoundList}}
#' @export
validateCompoundList <- function(compoundList, rt = FALSE, ccs = FALSE, extId = FALSE) {

  ################################################################################################################################################
  # Sanity Checks on input
  ################################################################################################################################################
  # correct column headers in compound list?
  if(!all(c("metaboliteID", "adductType", "adductMass", "neutralMass", "neutralFormula", "ionFormula", "metaboliteName", "inchikey", "inchi", "smiles") %in% colnames(compoundList))) {
    stop("One or more column header does not match the required headers: metaboliteID, adductType, adductMass, neutralMass, neutralFormula, ionFormula, metaboliteName, inchikey, inchi, smiles")
  }

  # check correct format of each column
  # TODO

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

  validatedCompoundList <- compoundList

  ################################################################################################################################################
  # add missing columns
  ################################################################################################################################################
  if(rt == FALSE) {
    validatedCompoundList$rt <- NA
  }

  if(ccs == FALSE) {
    validatedCompoundList$ccs <- NA
  }

  if(extId == FALSE) {
    validatedCompoundList$kegg <- NA
    validatedCompoundList$hmdb <- NA
    validatedCompoundList$chebi <- NA
  }

  return(validatedCompoundList)

}

#' Preparation of compound list
#'
#' A compound list that can be used with masstrixR can be generated on the fly. Minimum input is a data frame with the following columns: metabolite id ($id), SMILES ($smiles), InChI ($inchi), InChIKey ($inchikey), formula ($formula) metabolite name ($name) and an exact mass ($exactmass).
#' Furthermore, the adducts that shall be covered in the DB have to be defined. This can be either done by using an list of adducts or supplying a adduct definition for each metabolite with an additiona adduct column ($adducts). If it is intended to perform retetion time and collisional cross section matching columns containing this information are required ($rt and $ccs).
#' In case of CCS matching, individual adduct rules are required since each adduct has a different CCS value. Examples for each input are found in the examples in the vignettes.
#'
#' @param compoundList List of compounds that shall be added to DB
#' @param adductList Vector with adducts that shall be covered in the DB.
#' @param rt Boolean value indicating if compound list contains RT data
#' @param ccs Boolean value indicating if compound list contains CCS data
#'
#' @return Returns a data frame suitable for upload to a SQLite DB.
#'
#' @examples
#' xxx
#'
#' @author Michael Witting, \email{michael.witting@@helmholtz-muenchen.de}
#'
#' @seealso \code{\link{validateCompoundList}}
#' @seealso \code{\link{createDb}}
#' @export
prepareCompoundList <- function(compoundList, adductList = NA, rt = FALSE, ccs = FALSE, extId = FALSE, numberOfCores = 1) {

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

  # if(ccs == TRUE & !ccsType %in% c("DT", "TWIMS", "TIMS", "predicted")) {
  #   stop("Unknown CCS Type")
  # }



  ################################################################################################################################################
  # create required variables
  ################################################################################################################################################
  # data frame for final data
  newCompoundList <- data.frame()

  #get adduct calculation list
  adductCalc <- getAdductCalc()

  ###############################################################################################################################################
  # add missing columns to have common input format
  ###############################################################################################################################################
  # add rt column, even if false
  if(rt == FALSE) {
    print("no RT")
    compoundList$rt <- NA
  }

  if(ccs == FALSE) {
    print("no CCS")
    compoundList$ccs <- NA
  }

  ###############################################################################################################################################
  # check if external ids are present or fetch them from CTS
  ###############################################################################################################################################
  # are external ids part of the compound list
  if(extId == TRUE & !all(c("kegg", "hmdb", "chebi") %in% colnames(compoundList))) {
    message("external ID option was selected, but no KEGG, HMDB or ChEBI ID supplied. IDs will be fetched from CTS.")

    ################################################################################################################################################
    # make cluster for CTS
    ################################################################################################################################################
    # Initiate cluster
    cl <- parallel::makeCluster(numberOfCores)

    compoundList$kegg <- NA
    compoundList$hmdb <- NA
    compoundList$chebi <- NA

    # stop parallel cluster
    parallel::stopCluster(cl)

  } else {

    # check the different columns if they exist
    # KEGG
    if(!c("kegg") %in% colnames(compoundList)) {
      compoundList$kegg <- NA
    }

    # HMDB
    if(!c("hmdb") %in% colnames(compoundList)) {
      compoundList$hmdb <- NA
    }

    # ChEBI
    if(!c("chebi") %in% colnames(compoundList)) {
      compoundList$chebi <- NA
    }
  }

  ###############################################################################################################################################
  # check if external ids are present or fetch them from CTS
  ###############################################################################################################################################
  # create adduct
  if(c("adducts") %in% colnames(compoundList)) {

    print("individual adducts")

    #iterate over compound list and add to dbupload
    for(i in 1:nrow(compoundList)) {
      if(grepl(",", compoundList$adducts[i])) {
        adductList <- stringr::str_split(compoundList$adducts[i], ",")[[1]]
      } else {
        adductList <- compoundList$adducts[i]
      }

      #iterate over adducts and create a query list
      for(adduct in adductList) {

        #generate list
        clipboard <- data.frame(metaboliteID = compoundList$id[i],
                                adductType = adduct,
                                adductMass = compoundList$exactmass[i] * as.numeric(adductCalc[[adduct]][1]) + as.numeric(adductCalc[[adduct]][2]),
                                neutralMass = compoundList$exactmass[i],
                                neutralFormula = compoundList$formula[i],
                                ionFormula = calcAdductFormula(compoundList$formula[i], adduct),
                                metaboliteName = stringr::str_c(compoundList$name[i], adduct, sep = " "),
                                inchikey = compoundList$inchikey[i],
                                inchi = compoundList$inchi[i],
                                smiles = compoundList$smiles[i],
                                rt = compoundList$rt[i],
                                ccs = compoundList$ccs[i],
                                kegg = compoundList$kegg[i],
                                hmdb = compoundList$hmdb[i],
                                chebi = compoundList$chebi[i])

        #add to list
        newCompoundList <- rbind.data.frame(newCompoundList, clipboard)

      }
    }

  } else {

    # check if adduct names are correct
    if(!all(adductList %in% getAdductNames())) {
      stop("One or more adduct name does not match the required adduct names.")
    }

    print("commmon adducts")

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
                              inchikey = compoundList$inchikey,
                              inchi = compoundList$inchi,
                              smiles = compoundList$smiles,
                              rt = compoundList$rt,
                              ccs = compoundList$ccs,
                              kegg = compoundList$kegg,
                              hmdb = compoundList$hmdb,
                              chebi = compoundList$chebi)

      #add to list
      newCompoundList <- rbind.data.frame(newCompoundList, clipboard)

    }
  }

  return(newCompoundList)
}


getExternalDbIds <- function(inchikey, db) {

  from <- "inchikey"
  to <- db
  queryString <- inchikey

  #construct url for GET request
  baseUrl <- "http://cts.fiehnlab.ucdavis.edu/service/convert"
  queryUrl <- paste0(baseUrl, "/", from, "/", to, "/", queryString)

  result <- tryCatch(
    {
      queryResult <- jsonlite::fromJSON(URLencode(queryUrl))

      return(queryResult$results)

    },
    error=function(cond) {
      return("error")
    },
    warning=function(cond) {
      return("error")
    },
    finally = {
    }
  )

  print(result)

  return(result)

}

