
#search in a pre-defined database
mzSearch <-function(peakList, dbFileName, mode = "onDisk", mzTol = 0.005, tolType = "abs") {

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
    annotation <- dbFetch(resultSet)
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


#search in a in memory database created "on-the-fly" from a given compound list
mzLookUp <- function(peakList, compoundList, adductList, mzTol = 0.005, tolType = "abs") {

}
