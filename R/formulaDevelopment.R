
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



parsedChemFormula <- parseChemFormula("C6H12O6NaCaCl")


generateChemFormula <- function(parsedChemFormula) {

  # create empty string to append parts of formula
  chemFormulaRecon <- ""

  # first C H N O S P, then elements by alphabetical order
  for(atom in c("C", "H", "N", "O", "S", "P")) {
    if(!is.na(parsedChemFormula[atom])) {
      if(parsedChemFormula[atom] == 1) {
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
    if(parsedChemFormula[atom] == 1) {
      chemFormulaRecon <- paste0(chemFormulaRecon, atom)
    } else {
      chemFormulaRecon <- paste0(chemFormulaRecon, atom, parsedChemFormula[atom])
    }
  }

  # return formula
 return(chemFormulaRecon)
}










