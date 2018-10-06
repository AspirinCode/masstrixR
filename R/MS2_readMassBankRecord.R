#' Function to read a single MassBank record
#' Inspeired from MassBank parser from RMassBank
#' https://github.com/sneumann/RMassBank/blob/master/R/parseMassBank.R
#'
#'
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

  # get chemical names
  chnames <- list()
  chnames <- as.list(substring(grep('CH$NAME:',record, value = TRUE, fixed = TRUE),10))

  # parse adduct type?
  # ToDo

  #create new annotated Spectrum2
  mbRecord <- new("AnnotatedSpectrum2",
                  merged = 0,
                  precScanNum = as.integer(1),
                  precursorMz = as.numeric(substring(grep('MS$FOCUSED_ION: PRECURSOR_M/Z',record, value = TRUE, fixed = TRUE),31)),
                  precursorIntensity = 100,
                  precursorCharge = as.integer(1),
                  mz = unlist(PKPeak["m/z"]),
                  intensity = unlist(PKPeak["int"]),
                  centroided = TRUE,
                  collisionEnergy = as.numeric(substring(grep('AC$MASS_SPECTROMETRY: COLLISION_ENERGY',record, value = TRUE, fixed = TRUE),40)),
                  name = paste(chnames, sep = ";"),
                  formula = substring(grep('CH$FORMULA:',record, value = TRUE, fixed = TRUE),13),
                  exactMass = as.numeric(substring(grep('CH$EXACT_MASS:',record, value = TRUE, fixed = TRUE),16)),
                  inchi = "",
                  smiles = substring(grep('CH$SMILES:',record, value = TRUE, fixed = TRUE),12),
                  splash = "")

  return(mbRecord)
}
