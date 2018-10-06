library(stringr)
library(DBI)
library(RSQLite)

system.time({

  # where are the MB files
  pathToMBFiles <- "data\\"

  #where should the SQLiteDB be stored?
  pathToSQLite <- ":memory:"

  # get all files
  MBFiles <- list.files(pathToMBFiles, pattern = ".txt$", full.names = TRUE)

  # data frames for different information
  metaData <- data.frame()
  spectra <- data.frame()
  massSpec <- data.frame()

  # for internal id
  i <- 1

  # read each file and add to list
  for(file in MBFiles) {

    # create internal ID
    intID <- str_pad(i, 6, pad = "0")

    # create AnnotatedSpectrum2 object
    mbRecord <- readMassBankFile(file)

    # fill data frame
    metaData <- rbind.data.frame(metaData, cbind(id = intID,
                                                 name = paste0(mbRecord@name, collapse = ";"),
                                                 formula = mbRecord@formula,
                                                 exactMass = mbRecord@exactMass,
                                                 inchi = mbRecord@inchi,
                                                 inchiKey = mbRecord@inchiKey,
                                                 smiles = mbRecord@smiles))

    spectra <- rbind.data.frame(spectra, cbind(id = intID,
                                               mz = mbRecord@mz,
                                               int = mbRecord@intensity))

    massSpec <- rbind.data.frame(massSpec, cbind(id = intID,
                                                 precursorMz = mbRecord@precursorMz,
                                                 precursorIntensity = mbRecord@precursorIntensity,
                                                 msLevel = mbRecord@msLevel,
                                                 peaksCount = mbRecord@peaksCount))


    # for internal id
    i <- i + 1

  }

  # create DB

  #upload to temp DB
  mydb <- DBI::dbConnect(RSQLite::SQLite(), pathToSQLite)
  DBI::dbWriteTable(mydb, "metaData", metaData)
  DBI::dbWriteTable(mydb, "spectra", spectra)
  DBI::dbWriteTable(mydb, "massSpec", massSpec)

})

DBI::dbListTables(mydb)

DBI::dbGetQuery(mydb, 'SELECT * FROM spectra WHERE id = "000001"')
DBI::dbGetQuery(mydb, 'SELECT * FROM metaData WHERE id = "000001"')
DBI::dbGetQuery(mydb, 'SELECT * FROM massSpec WHERE id = "000001"')

system.time({

  ids <- unique(DBI::dbGetQuery(mydb, 'SELECT id FROM spectra WHERE mz BETWEEN 200 and 250'))

  for(id in ids$id) {

    # get spectrum metadata
    metadataRs <- dbSendQuery(mydb, 'SELECT * FROM metadata WHERE "id" = :x')
    dbBind(metadataRs, param = list(x = id))
    metaData <- dbFetch(metadataRs)

    dbClearResult(metadataRs)

    # get spectrum
    spectrumRs <- dbSendQuery(mydb, 'SELECT * FROM spectra WHERE "id" = :x')
    dbBind(spectrumRs, param = list(x = id))
    spectrum <- dbFetch(spectrumRs)

    dbClearResult(spectrumRs)

    # get mass spec details
    massSpecRs <- dbSendQuery(mydb, 'SELECT * FROM massSpec WHERE "id" = :x')
    dbBind(massSpecRs, param = list(x = id))
    massSpec <- dbFetch(massSpecRs)

    dbClearResult(massSpecRs)

    #create new annotated Spectrum2
    resultSpectrum <- new("AnnotatedSpectrum2",
                          merged = 0,
                          precScanNum = as.integer(1),
                          precursorMz = as.numeric(massSpec$precursorMz),
                          precursorIntensity = as.numeric(massSpec$precursorIntensity),
                          precursorCharge = as.integer(1),
                          mz = unlist(as.numeric(spectrum$mz)),
                          intensity = unlist(as.numeric(spectrum$int)),
                          centroided = TRUE,
                          collisionEnergy = 40,
                          name = metaData$name,
                          formula = metaData$formula,
                          exactMass = as.numeric(metaData$exactMass),
                          inchi = metaData$inchi,
                          smiles = metaData$smiles,
                          splash = "")

    print(resultSpectrum)
    plot(resultSpectrum, resultSpectrum)
  }

})

#disconnect DB
DBI::dbDisconnect(mydb)
