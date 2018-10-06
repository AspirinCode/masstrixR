
################################################################################################################################################
# With adduct definition for all metabolites
################################################################################################################################################

# with adduct definition for all metabolites
compoundList <- data.frame(read.table("DBs\\ExampleTxt\\ymdbExample_20180731.txt", header = T, sep = "\t", comment.char = "_"))

# correct error thrown??
dbFileName <- createDb(compoundList, adductList = c("M+H", "M+Na"),"testDB")

#connect to DB
mydb <- DBI::dbConnect(RSQLite::SQLite(), dbFileName)

# You can fetch all results:
res <- DBI::dbSendQuery(mydb, "SELECT * FROM adducts")
testResult <- DBI::dbFetch(res)

# Clear the result
DBI::dbClearResult(res)

# Disconnect from the database
DBI::dbDisconnect(mydb)

################################################################################################################################################
# With individual adduct definition
################################################################################################################################################
# with individual adduct definition
compoundListAdducts <- data.frame(read.table("DBs\\ExampleTxt\\ymdbExampleAdducts_20180731.txt", header = T, sep = "\t", comment.char = "_"))

# correct error thrown??
dbFileNameAdducts <- createDb(compoundListAdducts, "testDBAdducts")

#connect to DB
mydb <- DBI::dbConnect(RSQLite::SQLite(), dbFileNameAdducts)

# You can fetch all results:
res <- DBI::dbSendQuery(mydb, "SELECT * FROM adducts")
testResultAdducts <- DBI::dbFetch(res)

# Clear the result
DBI::dbClearResult(res)

# Disconnect from the database
DBI::dbDisconnect(mydb)

################################################################################################################################################
# CCS test without individual adduct definition
################################################################################################################################################
# with individual adduct definition
compoundListCcsNoAdduct <- data.frame(read.table("DBs\\ExampleTxt\\ccsExampleNoAdductDef_20180911.txt", header = T, sep = "\t", comment.char = "_"))

# correct error thrown??
#dbFileName <- createDb(compoundListCcsNoAdduct,"testDBCcs", adductList = c("M+H", "M+Na"), ccs = TRUE)

################################################################################################################################################
# CCS test with individual adduct definition
################################################################################################################################################
# with individual adduct definition
compoundListCcs <- data.frame(read.table("DBs\\ExampleTxt\\ccsExample_20180911.txt", header = T, sep = "\t", comment.char = "_"))

# correct error thrown??
dbFileName <- createDb(compoundListCcs,"testDBCcs", ccs = TRUE)

#connect to DB
mydb <- DBI::dbConnect(RSQLite::SQLite(), dbFileName)

# You can fetch all results:
res <- DBI::dbSendQuery(mydb, "SELECT * FROM adducts")
testResultCcs <- DBI::dbFetch(res)

# Clear the result
DBI::dbClearResult(res)

# Disconnect from the database
DBI::dbDisconnect(mydb)
