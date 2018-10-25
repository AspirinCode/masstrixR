#' This funtion parses a Genedata .gda file and returns three data frames containing the data and row and column annotations
#'
#' @param pathTofile
#' @return Returns a list with the three different dataframes, index 1 is the data, index 2 column annotation, index 3 row annotations
#' @examples
#' xxx
#' @export
readGdaFile <- function(pathToFile) {

  #create data frames for data
  rowAnnoDf <- data.frame(stringsAsFactors = FALSE)
  dataDf <- data.frame()
  clusterNames <- NULL

  # read file line by line
  con <- file(pathToFile, "r")

  while(TRUE) {
    line = readLines(con, n = 1)

    #exit loop
    if(length(line) == 0) {
      break
    }

    #check for different line content
    if(grepl("^Name", line)) {

      #get content splitted
      content <- stringr::str_split(line, "\t")[[1]]

      #isolate sample names and add to ColumnDf
      SampleNames <- as.character(content[2:(length(content) - noOfRowAnno)])
      colAnnoDf <- data.frame(SampleNames = character(length(content) - noOfRowAnno -1))
      colAnnoDf["SampleNames"] <- SampleNames

      #isolate row annotation names
      rowAnnotationNames <- as.character(content[(length(content) - noOfRowAnno + 1):length(content)])

    } else if(grepl("^# Row Annotations:", line)) {

      # split and get number of row annotations
      noOfRowAnno <- as.numeric(stringr::str_extract(stringr::str_split(line, ":")[[1]][2], "\\d+"))
      print(as.numeric(stringr::str_extract(stringr::str_split(line, ":")[[1]][2], "\\d+")))

      #print(line)

    }else if(grepl("\\[[A-Z]\\]", line)) {

      #get content splitted
      content <- stringr::str_split(line, "\t")[[1]]

      #get new annotation
      colAnnotationName <- content[[1]][1]
      colAnnotation <- content[2:(length(SampleNames) + 1)]

      #add to annotation data frame
      colAnnoDf[paste(colAnnotationName)] <- colAnnotation

    } else if(grepl("^(Peak|Cluster)", line)) {

      #get content splitted
      content <- stringr::str_split(line, "\t")[[1]]

      # get cluster names
      clusterNames <- c(clusterNames, content[1])

      #get peak values
      peakValues <- as.numeric(content[2:(length(content) - noOfRowAnno)])
      rowAnnotations <- content[(length(content) - noOfRowAnno + 1):length(content)]

      #add to rowAnnoDf and dataDf
      rowAnnoDf <- rbind.data.frame(rowAnnoDf,
                                    rowAnnotations, stringsAsFactors = FALSE)

      dataDf <- rbind.data.frame(dataDf,
                                 peakValues)
    }
  }

  # close connection to file
  close(con)

  #adjust column headers
  colnames(dataDf) <- SampleNames
  colnames(rowAnnoDf) <- rowAnnotationNames

  row.names(dataDf) <- clusterNames
  row.names(rowAnnoDf) <- clusterNames

  #adjust data type
  for(name in colnames(colAnnoDf)) {
    if(grepl("\\[N|n\\]", name)) {
      print(name)
      colAnnoDf[name] <- as.numeric(unlist(colAnnoDf[name]))
    } else if(grepl("\\[C|c\\]", name))
      colAnnoDf[name] <- as.factor(unlist(colAnnoDf[name]))
  }

  # change the m/z column names
  names(rowAnnoDf)[names(rowAnnoDf) == 'm/z'] <- 'm.z'
  rowAnnoDf$m.z <- as.numeric(rowAnnoDf$m.z)


  return(list(dataDf, colAnnoDf, rowAnnoDf))

}
