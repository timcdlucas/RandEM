context("remBoot")

test_that("remBoot produces a reasonable result", {
  # This is just a very basic regression test
  
  data(hDat)

  # hDat has detection in kilometres, we want metres.
  hDat$dist <- 1000 * hDat$dist
  tm <- 1800
  v <- 0.89
  nboots <- 1000

  output <- remBoot(hDat, tm, v, nboots, error_stat = c("sd"))
  
  expect_is(output, "list")
  expect_length(output, 1)
  expect_named(output, "sd")
  expect_length(output$sd, 5)
  expect_true(all(output$sd < 100))
})
