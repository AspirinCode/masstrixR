# # new class based on Spectrum2 for storage of library spectra
# setClass("AnnotatedSpectrum2",
#          representation = representation(
#            name = "character",
#            formula = "character",
#            exactMass = "numeric",
#            inchi = "character",
#            inchiKey = "character",
#            smiles = "character",
#            splash = "character"),
#          contains = c("Spectrum2"),
#          prototype = prototype()
# )
#
# #new class for consolidate spectra from samples
# setClass("ConsolidatedSpectrum2",
#          representation = representation(
#            chromPeakId = "character",
#            deltaRt = "numeric",
#            deltaMz = "numeric",
#            featureId = "character"),
#          contains = c("Spectrum2"),
#          prototype = prototype()
# )
#
#
# #new class that stores MS1 and MS2 information in one object, e.g. for SIRIUS export
# setClass("featureSpectra",
#          representation = representation(
#            xxx = "character",
#            yyy = "character",
#            zzz = "character"
#          ),
#          contains = c("Spectrum1", "ConsolidatedSpectrum2"))
