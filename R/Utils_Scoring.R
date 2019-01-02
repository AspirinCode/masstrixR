scoringScheme <- function(annotationResult, scoringType = "masstrixR", scoreTable) {

  # implement here a table with scoring
  # this scoring normalizes all scores from different tools to the same scale

  # MS1 mz score

  # MS1 rt score

  # MS1 isotope score

  # MS1 Sirius

  # MS2 library search score

  # MS2 CSI:FingerID

  # MS2 MetFrag

  # all results shall be normalized on first part of InChIKey


  return(TRUE)
}

plotScores <- function(scoreTable) {



}


formulaConsensusScoring <- function(resultDf) {

  # reformat data frame
  resultDf$formula <- unlist(lapply(resultDf$formula, harmonizeFormula))

  ggplot(resultDf, aes(x = formula, y = tool, fill = score)) + geom_tile() +
    scale_fill_gradient(limits = c(0,100),low = "red", high = "green", space = "Lab", na.value = "grey50", guide = "colourbar", aesthetics = "fill") +
    theme_bw()

}

ms1ConsensusScoring <- function(resultDf) {

}


ms2ConsensusScoring <- function(resultDf) {

}


readScoreTable <- function(pathToFile) {

  # read from file
  scoreTable <- data.frame(read.table(pathToFile, header = FALSE, stringsAsFactors = FALSE))

  # split columns
  scoreTable <- tidyr::separate(scoreTable, "V1", into = c("scoringType", "V1"), sep = "\\.", extra = "merge")
  scoreTable <- tidyr::separate(scoreTable, "V1", into = c("variable", "V1"), sep = "_", extra = "merge")
  scoreTable <- tidyr::separate(scoreTable, "V1", into = c("errorType", "V1"), sep = "_", extra = "merge")
  scoreTable <- tidyr::separate(scoreTable, "V1", into = c("margin", "score"), sep = "=", extra = "merge")

  # convert to correct data type
  scoreTable$score <- as.numeric(scoreTable$score)

  # sanity check on scoreTable
  # TODO


  #return final score table
  return(scoreTable)
}
