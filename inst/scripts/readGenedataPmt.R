# get file containing peak match table
peakMatchTableFile <-
  system.file("extdata", "MS2/NaAc_Test.pmt", package = "masstrixR")

# read the peak match table
peakMatchTable <- readPmt(peakMatchTable)

# list all mgf files
mgfFiles <- list.files(pattern = ".mgf$")

# get all samples from the PMT
samples <- as.character(unique(peakMatchTable$sampleName))

# remove blanks
samples <- samples[-grep("blank", samples)]

# make new Spectra object
fullMs2SpectraList <- new("Spectra")

# iterate through samples and append MS2 spectra with respective cluster ID
for(sample in samples) {

  print(sample)

  # read mgf file with MS2 spectra and isolate all MS2 spectra a list
  ms2data <- readMgfData(paste0(sample, ".mgf", sep = ""))
  ms2spectra <- spectra(filterMsLevel(ms2data, msLevel = 2))

  # get all the peaks in the respective sample
  filteredPeakMatchTable <-
    peakMatchTable[which(peakMatchTable$sampleName == sample), ]

  # create new Spectra object and add scan index as annotation
  ms2SpectraList <- Spectra(ms2spectra)
  mcols(ms2SpectraList)$scanIndex <-
    unlist(lapply(ms2spectra, scanIndex))

  # add cluster
  mcols(ms2SpectraList)$cluster <-
    as.character(unlist(lapply(mcols(ms2SpectraList)$scanIndex, function(x) {
      return(filteredPeakMatchTable$cluster[which(filteredPeakMatchTable$content == x)])
    })))

  fullMs2SpectraList <- append(fullMs2SpectraList, ms2SpectraList)
}
