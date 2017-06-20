context("rem")

test_that("rem produces a reasonable result", {
  # This is just a very basic regression test
  
  data(hDat)
  grpDat <- split_dat(hDat)
  tm <- 3600
  v <- 1.4
  
  # Test an example output against known outputs and check that the difference
  densities <- lapply(grpDat, rem, tm, v)
  expect_is(densities, "list")
  expect_length(densities, 5)
  
  # "Known" values. Taken from running `lapply(grpDat, rem, tm, v)` and saving
  # the values with `dput`
  known_output <- c(2.24999473853356, 
                    1.07668613384841, 
                    1.86016465692493, 
                    1.45899308473601, 
                    1.52149811808646)
  max_resid <- max(abs(unlist(densities) - known_output))
  
  expect_true(max_resid < 6e-15)
})
