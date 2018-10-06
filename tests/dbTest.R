#test function
testDB <- function(){
  #read HMDB file and create DB
  compoundList <- data.frame(read.table("DBs\\Txt\\ymdb_20180731.txt", sep = "\t", header = T, stringsAsFactors = F, comment.char = "_"))
  dbPath <- masstrixR::createDb(compoundList, "DBs\\SQLite\\ymdb_20180731_MH_MNa", c("M+H", "M+Na"))

  #connect to DB
  mydb <- DBI::dbConnect(RSQLite::SQLite(), dbPath)

  #execute query
  resultSet <- DBI::dbSendQuery(mydb, 'SELECT * FROM adducts WHERE adductMass BETWEEN 200 AND 250')

  #get length of query
  length <- nrow(DBI::dbFetch(resultSet))
  DBI::dbDisconnect(mydb)

  return(length)
}
