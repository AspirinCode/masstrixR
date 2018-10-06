

####
# this function is used to generate the table that is uploaded to the SQLite DB
####

.generateSQLiteInput <- function(compoundList, adductList = NA, rt = FALSE, ccs = FALSE) {

  #get adduct calculation list
  adductCalc <- getAdductCalc()

  #create DF for upload
  dbupload <- data.frame()

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
                                ccs = ccsValue)

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
  }

  return(dbupload)

}
