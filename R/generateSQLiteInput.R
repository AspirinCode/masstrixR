
####
# this function is used to generate the table that is uploaded to the SQLite DB
####
.generateSQLiteInput <- function(compoundList, adductList = NA, rt = FALSE, ccs = FALSE, fetchIds = FALSE) {

  #get adduct calculation list
  adductCalc <- getAdductCalc()

  #create DF for upload
  dbupload <- data.frame()

  print(fetchIds)

  ################################################################################################################################################
  # Create list for upload to SQLite DB
  ################################################################################################################################################
  # check if adducts are predefined or not
  if(c("adducts") %in% colnames(compoundList)) {

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

        if(rt == TRUE) {
          rtime <- compoundList$rt[i]
        } else {
          rtime <- NA
        }

        if(ccs == TRUE) {
          ccsValue <- compoundList$ccs[i]
        } else {
          ccsValue <- NA
        }

        # check which columns with external ids are present
        if(fetchIds) {

          kegg <- .getExternalDbIds(compoundList$inchikey[i], "KEGG")
          hmdb <- .getExternalDbIds(compoundList$inchikey[i], "HMDB")
          chebi <- .getExternalDbIds(compoundList$inchikey[i], "ChEBI")

          # tests
          #print("fetching ids")

        } else {

          # check the different columns if they exist
          # KEGG
          if("kegg" %in% colnames(compoundList)) {
            kegg <- compoundList$kegg[i]
          } else {
            kegg <- NA
          }

          # HMDB
          if("hmdb" %in% colnames(compoundList)) {
            hmdb <- compoundList$hmdb[i]
          } else {
            hmdb <- NA
          }

          # ChEBI
          if("chebi" %in% colnames(compoundList)) {
            chebi <- compoundList$chebi[i]
          } else {
            chebi <- NA
          }

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
                                rt = rtime,
                                ccs = ccsValue,
                                kegg = kegg,
                                hmdb = hmdb,
                                chebi = chebi)

        #add to list
        dbupload <- rbind.data.frame(dbupload, clipboard)

      }
    }

  } else {

    # check if adduct names are correct
    if(!all(adductList %in% getAdductNames())) {
      stop("One or more adduct name does not match the required adduct names.")
    }

    #iterate over adducts and create a query list
    for(adduct in adductList) {

      # add rt if TRUE
      if(rt == TRUE) {
        ritmes <- compoundList$rt
      } else {
        rtimes <- rep(0, nrow(compoundList))
      }

      # check which columns with external ids are present
      if(fetchIds) {

        kegg <- unlist(lapply(compoundList$inchikey, .getExternalDbIds, db = "KEGG"))
        hmdb <- unlist(lapply(compoundList$inchikey, .getExternalDbIds, db = "HMDB"))
        chebi <- unlist(lapply(compoundList$inchikey, .getExternalDbIds, db = "ChEBI"))

      } else {

        # check the different columns if they exist
        # KEGG
        if("kegg" %in% colnames(compoundList)) {
          kegg <- compoundList$kegg
        } else {
          kegg <- NA
        }

        # HMDB
        if("hmdb" %in% colnames(compoundList)) {
          hmdb <- compoundList$hmdb
        } else {
          hmdb <- NA
        }

        # ChEBI
        if("chebi" %in% colnames(compoundList)) {
          chebi <- compoundList$chebi
        } else {
          chebi <- NA
        }

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
                              ccs = NA,
                              kegg = kegg,
                              hmdb = hmdb,
                              chebi = chebi)

      #add to list
      dbupload <- rbind.data.frame(dbupload, clipboard)

    }
  }

  return(dbupload)

}

.getExternalDbIds <- function(inchikey, db) {

  from <- "inchikey"
  to <- db
  queryString <- inchikey

  #construct url for GET request
  baseUrl <- "http://cts.fiehnlab.ucdavis.edu/service/convert"
  queryUrl <- paste0(baseUrl, "/", from, "/", to, "/", queryString)

  result <- tryCatch(
    {
      queryResult <- jsonlite::fromJSON(URLencode(queryUrl))
      print(queryResult$result[1])
      return(queryResult$result[1])
    },
    error=function(cond) {
      return("error")
    },
    warning=function(cond) {
      return("error")
    },
    finally = {
      return("error")
    }
  )

  return(result)

}
