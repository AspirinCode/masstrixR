#' This function generates the formula of a pseudomolecular ion from a chemical formula and the type of adduct selected.
#'
#' @return Returns a string with the ion formula
#' @examples
#' ionFormula <- calcAdductFormula("C6H12O6", "M+H")
calcAdductFormula <- function(molecularFormula, adduct) {

  # get adduct calculation list
  adductCalc <- getAdductCalc()

  #print(paste0(molecularFormula, "_", adduct))


  # molecule formula parsed
  if(as.numeric(adductCalc[[adduct]][1]) > 1) {
    molFormulaList <- parseChemFormula(molecularFormula) * as.numeric(adductCalc[[adduct]][1]) # multiplied with factor from adduct calculation
  } else {
    molFormulaList <- parseChemFormula(molecularFormula)
  }


  # get additive and subtractive part for formula
  addFormulaList <- parseChemFormula(adductCalc[[adduct]][3])
  subFormulaList <- parseChemFormula(adductCalc[[adduct]][4])

  # get all involved atoms
  atoms <- unique(c(names(molFormulaList), names(addFormulaList), names(subFormulaList)))

  # create an empty list with all atoms
  adductFormulaList <- sapply(atoms, function(x) NULL)

  for(i in 1:length(adductFormulaList)) {
    #print(names(adductFormulaList[i]))

    atom <- names(adductFormulaList[i])

    mol <- if(!is.na(molFormulaList[atom])) {molFormulaList[atom]} else {0}
    add <- if(!is.na(addFormulaList[atom])) {addFormulaList[atom]} else {0}
    sub <- if(!is.na(subFormulaList[atom])) {subFormulaList[atom]} else {0}

    adductFormulaList[i] <- mol + add - sub
  }

  # generate ion formula
  ionFormula <- generateChemFormula(unlist(adductFormulaList))

  # return generate formula
  return(ionFormula)
}

#' This function checks if one sum formula is contained in another.
#'
#' @return TRUE or FALSE
#' @examples
#' ionFormula <- containsFormula("C6H12O6", "H2O")
containsFormula <- function(targetFormula, queryFormula) {

  # parse both formmula
  targetFormulaList <- parseChemFormula(targetFormula)
  queryFormulaList <- parseChemFormula(queryFormula)

  # get atoms from query formula to check
  atoms <- names(queryFormulaList)

  # return value
  contains <- FALSE

  for(atom in atoms) {

    if(is.na(targetFormulaList[atom])) {
      return(FALSE)
    } else if(targetFormulaList[atom] - queryFormulaList[atom] >= 0) {
      contains <- TRUE
    }
  }

  return(contains)

}

#' This function parses a chemical formula into a named vector
#'
#' @return Named vector with all elements
#' @examples
#' parseChemFormula("C6H12O6")
parseChemFormula <- function(chemFormula) {

  #regex pattern to isolate all elements
  elementPattern <- "([A][cglmrstu]|[B][aehikr]?|[C][adeflmnorsu]?|[D][bsy]|[E][rsu]|[F][elmr]?|[G][ade]|[H][efgos]?|[I][nr]?|[K][r]?|[L][airuv]|[M][cdgnot]|[N][abdehiop]?|[O][gs]?|[P][abdmortu]?|[R][abefghnu]|[S][bcegimnr]?|[T][abcehilms]|[U]|[V]|[W]|[X][e]|[Y][b]?|[Z][nr])([0-9]*)"

  # extract all matching pattern
  regexMatch <- stringr::str_extract_all(chemFormula, elementPattern)

  # get individual elements and their count
  elements <- stringr::str_extract(regexMatch[[1]], "[aA-zZ]+")
  numbers <- as.numeric(stringr::str_extract(regexMatch[[1]], "[0-9]+"))

  # replace NAs with 1 for elements which have a count of one
  numbers[is.na(numbers)] <- 1

  # create named vector for return
  parsedChemFormula <- numbers
  names(parsedChemFormula) <- elements

  return(parsedChemFormula)
}

#' This function generates a chemical formula from a named vector of elemental counts
#'
#' @return Single string with chemical formula
#' @examples
#' parseChemFormula("C6H12O6")
generateChemFormula <- function(parsedChemFormula) {

  # create empty string to append parts of formula
  chemFormulaRecon <- ""

  # first C H N O S P, then elements by alphabetical order
  for(atom in c("C", "H", "N", "O", "S", "P")) {
    if(!is.na(parsedChemFormula[atom])) {
      if(parsedChemFormula[atom] == 1.0) {
        chemFormulaRecon <- paste0(chemFormulaRecon, atom)
      } else {
        chemFormulaRecon <- paste0(chemFormulaRecon, atom, parsedChemFormula[atom])
      }
    }
  }

  # get all remaining elements
  restElements <- names(parsedChemFormula)
  restElements <- restElements[! restElements %in% c("C", "H", "N", "O", "S", "P")]

  # iterate through all remaining elements in alphabetical order
  for(atom in sort(restElements)) {
    if(parsedChemFormula[atom] == 1.0) {
      chemFormulaRecon <- paste0(chemFormulaRecon, atom)
    } else {
      chemFormulaRecon <- paste0(chemFormulaRecon, atom, parsedChemFormula[atom])
    }
  }

  # return formula
  return(chemFormulaRecon)
}


#' This function generates a chemical formula from a named vector of elemental counts
#'
#' @return Single string with chemical formula
#' @examples
#' parseChemFormula("C6H12O6")
harmonizeFormula <- function(chemFormula) {

  # parse and generate
  parsedChemFormula <- parseChemFormula(chemFormula)
  chemFormula <- generateChemFormula(parsedChemFormula)

  return(chemFormula)


}
