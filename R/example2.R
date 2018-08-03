

exampleGDA <- readGdaFile("example\\NaAcHILICPos_Peaks.gda")

exampleGDA[[3]]

# Cluster 1262
testmz <- c(118.0867417, 119.0896037, 120.0904638)
testInt <- c(412946.66, 25346.621, 1890.015)

measured <- new("Spectrum1",
                mz = testmz,
                intensity = testInt,
                centroided = TRUE)

plot(testmz, testInt, type = "h")


ionFormula <- "C5H12NO2"

theoretical <- masstrixR::generateIsoPattern(ionFormula, "M+H")

compareSpectra(measured, theoretical, fun = "dotproduct")

plot(measured)

