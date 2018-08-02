# tests

# http://r-pkgs.had.co.nz/tests.html
# http://r-pkgs.had.co.nz/tests.html
#http://r-pkgs.had.co.nz/description.html


test_that("adduct formula calculation", {
  expect_equal(calcAdductFormula("C6H12O6", "M+Na"), "C6H12O6Na")
})
