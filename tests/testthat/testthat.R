# tests

# http://r-pkgs.had.co.nz/tests.html
# http://r-pkgs.had.co.nz/tests.html
# http://r-pkgs.had.co.nz/description.html


test_that("adduct formula calculation", {
  expect_equal(calcAdductFormula("C6H12O6", "M+Na"), "C6H12O6Na")
})


test_that("isotope mass calculation", {
  expect_equal(calculateIsoLabelMass_Formula("C6H12O6", "C"), 186.0835)
})

test_that("isotope mass calculation", {
  expect_equal(calculateIsoLabelMass_Formula("C6H12O6", "D"), 192.1387)
})

test_that("isotope mass calculation", {
  expect_equal(calculateIsoLabelMass_Formula("C6H12O6", "T"), 204.1621)
})

test_that("isotope mass calculation", {
  expect_equal(calculateIsoLabelMass_Formula("C6H12O6", "C", noOfLabel = 3), 183.0735)
})

test_that("isotope mass calculation", {
  expect_equal(calculateIsoLabelMass_Formula("C6H12O6", "D", noOfLabel = 6), 186.1011)
})

test_that("isotope mass calculation", {
  expect_equal(calculateIsoLabelMass_Formula("C6H12O6", "T", noOfLabel = 2), 184.0798)
})



