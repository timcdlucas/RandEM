context("rem")

test_that("rem produces a reasonable result", {
  # This is just a very basic regression test
  
  data(hDat)
  hDat$dist <- 1000 * hDat$dist
  grpDat <- split_dat(hDat)
  tm <- 3600
  v <- 1.4
  
  # Test an example output against known outputs and check that the difference
  densities <- lapply(grpDat, rem, tm, v)
  expect_is(densities, "list")
  expect_length(densities, 5)
  
  # "Known" values. Taken from running `lapply(grpDat, rem, tm, v)` and saving
  # the values with `dput`
  known_output <- c(53.99987, 25.84047, 44.64395, 35.01583, 36.51595)
  max_resid <- max(abs(unlist(densities) - known_output))
  
  expect_true(max_resid < 6e-3)
  
})
