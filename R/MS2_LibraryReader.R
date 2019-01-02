#' Function to read a single MassBank record
#' Inspeired from MassBank parser from RMassBank
#' https://github.com/sneumann/RMassBank/blob/master/R/parseMassBank.R
#'
#'
readMassBankFile <- function(pathToMBFile) {

  require(MSnbase)

  #read file
  fileConnection <- file(pathToMBFile)
  record <- readLines(fileConnection)
  close(fileConnection)

  #read the spectrum
  PKStart <- grep('PK$PEAK:',record, fixed = TRUE) + 1
  endslash <- tail(grep('//',record, fixed = TRUE),1)

  if(PKStart < endslash){
    splitted <- strsplit(record[PKStart:(endslash-1)]," ")
    PKPeak <- matrix(nrow = endslash - PKStart, ncol = 3)

    for(k in 1:length(splitted)){
      splitted[[k]] <- splitted[[k]][which(splitted[[k]] != "")]
      PKPeak[k,] <- splitted[[k]]
    }
    PKPeak <- as.data.frame(PKPeak, stringsAsFactors = FALSE)
    PKPeak[] <- lapply(PKPeak, type.convert)
    colnames(PKPeak) <- c("m/z", "int", "rel.int.")
  }


  # get basic information
  id <- substring(grep("ACCESSION:", record, value = TRUE, fixed = TRUE), 12)

  # get information of chemical compound
  chnames <- list()
  chnames <- as.list(substring(grep('CH$NAME:',record, value = TRUE, fixed = TRUE),10))
  formula <- substring(grep("CH$FORMULA:", record, value = TRUE, fixed = TRUE), 13)
  exactMass <- as.numeric(substring(grep("CH$EXACT_MASS:", record, value = TRUE, fixed = TRUE), 16))
  smiles <- substring(grep('CH$SMILES:',record, value = TRUE, fixed = TRUE), 12)
  inchi <- substring(grep('CH$IUPAC:', record, value = TRUE, fixed = TRUE), 11)

  # get instrument information
  instrument <- substring(grep("AC$INSTRUMENT:", record, value = TRUE, fixed = TRUE), 15)
  instrumentType <- substring(grep("AC$INSTRUMENT_TYPE:", record, value = TRUE, fixed = TRUE), 20)

  # get mass spec and ion information
  msType <- substring(grep('AC$MASS_SPECTROMETRY: MS_TYPE', record, value = TRUE, fixed = TRUE), 31)
  ionMode <- substring(grep('AC$MASS_SPECTROMETRY: ION_MODE', record, value = TRUE, fixed = TRUE), 32)
  collisionEnergy <- as.numeric(substring(grep('AC$MASS_SPECTROMETRY: COLLISION_ENERGY', record, value = TRUE, fixed = TRUE),40))

  # get information of precursor
  precursorMz <- as.numeric(substring(grep('MS$FOCUSED_ION: PRECURSOR_M/Z', record, value = TRUE, fixed = TRUE),31))
  adduct <- substring(grep('MS$FOCUSED_ION: PRECURSOR_TYPE', record, value = TRUE, fixed = TRUE), 32)

  # get peak information
  splash <- substring(grep("PK$SPLASH:", record, value = TRUE, fixed = TRUE), 11)
  numPeak <- as.numeric(substring(grep("PK$NUM_PEAK:", record, value = TRUE, fixed = TRUE), 13))


  # make Spectrum2 oboject from data
  ms2spec <- new("Spectrum2",
                 merged = 0,
                 precScanNum = as.integer(1),
                 precursorMz = precursorMz,
                 precursorIntensity = 100,
                 precursorCharge = as.integer(1),
                 mz = unlist(PKPeak["m/z"]),
                 intensity = unlist(PKPeak["int"]),
                 collisionEnergy = collisionEnergy,
                 centroided = TRUE)

  # make new Spectra object
  mbRecord <- MSnbase::Spectra(ms2spec)

  # add annotations
  mcols(mbRecord)$id <- id
  mcols(mbRecord)$name <- paste(unlist(chnames), collapse = ";")
  mcols(mbRecord)$formula <- formula
  mcols(mbRecord)$exactMass <- exactMass
  mcols(mbRecord)$smiles <- smiles
  mcols(mbRecord)$inchi <- inchi
  mcols(mbRecord)$instrument <- instrument
  mcols(mbRecord)$instrumentType <- instrumentType
  mcols(mbRecord)$msType <- msType
  mcols(mbRecord)$ionMode <- ionMode
  mcols(mbRecord)$precursorMz <- precursorMz
  mcols(mbRecord)$precursorType <- adduct
  mcols(mbRecord)$splash <- splash
  mcols(mbRecord)$numPeak <- numPeak

  return(mbRecord)
}

#' Reada a MSP library
#'
# parseMSP_attributes <- function(fileSpectra, progress = FALSE, flexiblePeakList = FALSE, includeIDasRecordSeparator=TRUE, returnEmptySpectra = FALSE){
#   fileLines <- readLines(con = fileSpectra)
#
#   if(!is.na(progress))  if(progress)  incProgress(amount = 0, detail = "MS/MS file: Read file") else print("MS/MS file: Read file")
#
#   #fileLines <- readLines(con = fileSpectra)
#   numberOfFileLines <- length(fileLines)
#
#   ## start with empty lines or not?
#   endOfRecord <- TRUE
#   if(numberOfFileLines > 0)
#     if(nchar(trimws(fileLines[[1]])) > 0)
#       endOfRecord <- FALSE
#
#   ## check for pattern
#   if(!is.na(progress))  if(progress)  incProgress(amount = 0, detail = "MS/MS file: Parse") else print("MS/MS file: Parse")
#   isID      <- grepl(pattern = "^ID:",                                    x = fileLines)
#   isBI      <- grepl(pattern = "^BEGIN IONS$",                            x = fileLines)
#   isName    <- grepl(pattern = "(^Name:)|(^NAME:)",                       x = fileLines)
#   isNAme    <- grepl(pattern = "^NAME=",                                   x = fileLines)
#   isNAME    <- grepl(pattern = "^TITLE=",                                  x = fileLines)
#   #isNumP    <- grepl(pattern = "^Num Peaks:",                                        x = fileLines)
#   #isPeak    <- grepl(pattern = "^\\d+(\\.\\d+)?[ \t]\\d+(\\.\\d+)?$",    x = fileLines)
#   #isPeak    <- grepl(pattern = "^[ \t]*\\d+(\\.\\d+)?[ \t]\\d+(\\.\\d+)?([ \t]+\\d+(\\.\\d+)?[ \t]\\d+(\\.\\d+)?)*[ \t]*$",    x = fileLines)
#
#   if(!includeIDasRecordSeparator) isID <- rep(x = F, times = length(isID))
#
#   numberSmall <- "\\d+(\\.\\d+)?"
#   numberBig   <- "\\d+(\\.\\d+(E\\d+)?)?"
#   mzValueRegEx <- "(\\d+(\\.\\d+)?)"
#   intensityRegex <- "(\\d+((\\.\\d+)?([eE](-)?\\d+)?)?)"
#   annotationRegex <- "\".+\""
#   if(flexiblePeakList){
#     isPeak    <- grepl(pattern = "^[ \t]*\\d+(\\.\\d+([eE](-)?\\d+)?)?([ \t]\\d+(\\.\\d+([eE](-)?\\d+)?)?)*[ \t]*$",    x = fileLines)
#   } else {
#     #isPeak    <- grepl(pattern = "^[ \t]*\\d+(\\.\\d+)?[ \t]\\d+(\\.\\d+([eE](-)?\\d+)?)?([ \t]+((\\d+(\\.\\d+)?[ \t]\\d+(\\.\\d+([eE](-)?\\d+)?)?)|(\".+\")))*[ \t]*$",    x = fileLines)
#     isPeak    <- grepl(pattern = paste("^[ \t]*", mzValueRegEx, "[ \t]", intensityRegex, "([ \t]+((", mzValueRegEx, "[ \t]", intensityRegex, ")|(", annotationRegex, ")))*[ \t]*$", sep = ""),    x = fileLines)
#   }
#   isEmpty    <- nchar(trimws(fileLines)) == 0
#
#   tagVector   <- unlist(lapply(X = str_split(string = fileLines, pattern = "(:)|(=)"), FUN = function(x){x[[1]]}))
#   valueVector <- trimws(substr(x = fileLines, start = nchar(tagVector) + 1 + 1, stop = nchar(fileLines)))
#
#   ## entry line intervals in file
#   entryBorders   <- c(which(isName | isNAME | isNAME | isBI | isID), length(fileLines)+1)
#   entryIntervals <- matrix(data = unlist(lapply(X = seq_len(length(entryBorders) - 1), FUN = function(x){c(entryBorders[[x]], entryBorders[[x+1]] - 1)})), nrow=2)
#
#   ## do it
#   if(!is.na(progress))  if(progress)  incProgress(amount = 0, detail = "MS/MS file: Assemble spectra") else print("MS/MS file: Assemble spectra")
#   suppressWarnings(
#     spectraList <- apply(X = entryIntervals, MARGIN = 2, FUN = function(x){
#       #print(x)
#       fileLines2   <- fileLines  [x[[1]]:x[[2]]]
#       isPeak2         <- isPeak       [x[[1]]:x[[2]]]
#       isEmpty2     <- isEmpty    [x[[1]]:x[[2]]]
#
#       tagVector2   <- tagVector  [x[[1]]:x[[2]]]
#       valueVector2 <- valueVector[x[[1]]:x[[2]]]
#
#       ###################################################################
#       ## built ms set
#       spectrumItem <- list()
#       spectrumItem[tagVector2[!isPeak2 & !isEmpty2]] <- trimws(valueVector2[!isPeak2 & !isEmpty2])
#
#       if(!is.null(spectrumItem$"Num Peaks"))
#         if(spectrumItem$"Num Peaks" == "0" & !returnEmptySpectra)
#           return(NULL)
#
#       #spectrumItem["peaks"] <- paste(fileLines2[isPeak2], collapse = "; ")
#       peakLines <- fileLines2[isPeak2]
#       peakLines <- trimws(gsub(x = peakLines, pattern = "\".*\"", replacement = ""))
#       spectrumItem["peaks"] <- paste(peakLines, collapse = " ")
#       spectrumItem["peaks"] <- trimws(gsub(x = spectrumItem["peaks"], pattern = "  ", replacement = " "))
#
#       ## check peaks
#       tokens <- strsplit(x = spectrumItem[["peaks"]], split = "[ \t]")[[1]]
#       if(length(tokens) == 0){#spectrumItem$"Num Peaks" == "0"){
#         mzs  <- character(0)
#         ints <- character(0)
#       } else {
#         mzs  <- tokens[seq(from=1, to=length(tokens), by = 2)]
#         ints <- tokens[seq(from=2, to=length(tokens), by = 2)]
#       }
#       if(length(mzs) != length(ints)) stop("error in parsing peaks")
#
#       ## handle duplicated tags
#       duplicatedTags    <- unique(names(spectrumItem)[duplicated(names(spectrumItem))])
#       if(length(duplicatedTags) > 0){
#         duplicated <- sapply(X = duplicatedTags, FUN = function(x){
#           unlist(sapply(X = seq_along(spectrumItem), FUN = function(y){ if(names(spectrumItem[y])==x) return(y) }))
#         }, simplify = F)
#
#         indecesToRemove <- vector(mode = "integer", length = 0)
#         for(idx in seq_along(duplicated)){
#           indeces <- duplicated[[idx]]
#           representant <- indeces[[1]]
#           indecesToRemoveHere <- indeces[-1]
#
#           spectrumItem[[representant]] <- paste(spectrumItem[indeces], sep = "; ")
#           indecesToRemove <- c(indecesToRemove, indecesToRemoveHere)
#         }
#         spectrumItem <- spectrumItem[-indecesToRemove]
#       }
#
#       return(spectrumItem)
#     })
#   )## suppressWarnings
#
#   rm(
#     isName,
#     isNAME,
#     isNAme,
#     isBI,
#     isID,
#     isPeak,
#     tagVector,
#     valueVector
#   )
#
#   if(!is.na(progress))  if(progress)  incProgress(amount = 0.1, detail = "MS/MS file: Box") else print("MS/MS file: Box")
#
#   ## remove NULL entries?
#   spectraList[unlist(lapply(X = spectraList, FUN = is.null))] <- NULL
#
#   numberOfSpectra <- length(spectraList)
#
#   ## postprocess
#   if(!is.na(progress))  if(progress)  incProgress(amount = 0.01, detail = paste("MS/MS file postprocessing", sep = "")) else print(paste("MS/MS file postprocessing", sep = ""))
#
#   if(!is.na(progress))  if(progress)  incProgress(amount = 0.01, detail = paste("MS/MS file boxing", sep = "")) else print(paste("MS/MS file boxing", sep = ""))
#   returnObj <- list()
#   returnObj$fileSpectra <- fileSpectra
#   returnObj$spectraList <- spectraList
#   returnObj$numberOfSpectra <- numberOfSpectra
#
#   return(returnObj)
# }
