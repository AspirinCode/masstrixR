#' Function to read a single MassBank record
#' Inspeired from MassBank parser from RMassBank
#' https://github.com/sneumann/RMassBank/blob/master/R/parseMassBank.R
#'
#' @param pathToMBFile File path to a single MassBank Record
#'
#' @examples
#' readMassBankFile()
#'
#' @author Michael Witting, \email{michael.witting@@helmholtz-muenchen.de}
#'
#' @importClassesFrom MSnbase Spectrum2 Spectra
#' @export
readMassBankFile <- function(pathToMBFile) {

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

#' Function to read a single MassBank record
#' Inspeired from MassBank parser from RMassBank
#' https://github.com/sneumann/RMassBank/blob/master/R/parseMassBank.R
#'
#' @param pathToFolder File path to a single MassBank Record
#'
#' @examples
#' readMassBankFolder()
#'
#' @author Michael Witting, \email{michael.witting@@helmholtz-muenchen.de}
#' @importClassesFrom MSnbase Spectra
#' @export
readMassBankFolder <- function(pathToFolder) {

  #read files in folder
  massBankFiles <- list.files(pathToFolder, pattern = ".txt$", full.names = TRUE)

  #make empty spectra
  librarySpectra <- new("Spectra")

  #iterate through files
  for(massBankFile in massBankFiles) {
    spectrum <- readMassBankFile(massBankFile)

    # append to spectra list
    librarySpectra <- append(librarySpectra, spectrum)
  }

  # return list with MassBank Spectra
  return(librarySpectra)

}
