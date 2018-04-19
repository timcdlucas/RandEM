
context("Test grem functions")


# calcProfileWidth() tests

test_that('Profile width is within 0 and 2*r', {
        alpha <- seq(0, 2*pi, length.out=100)
        theta <- seq(0, 2*pi, length.out=100)
        paras <- expand.grid(alpha, theta)
        paras <- cbind(paras, 1)
        colnames(paras) <- c('theta','alpha','r')
        p <- sapply(1:nrow(paras), function(x) do.call(calcProfileWidth,as.list(paras[x,])))
        expect_that(sum(p>2*pi) == 0, is_true())
        expect_that(sum(p<0) == 0, is_true())
})


test_that('Profile width is 0 when alpha is 0', {
        theta <- seq(0, 2*pi, length.out=100)
        paras <- cbind(theta, 0, 1)
        colnames(paras) <- c('theta','alpha','r')
        p <- sapply(1:nrow(paras), function(x) do.call(calcProfileWidth,as.list(paras[x,])))
        expect_that(all(p==0), is_true())
})


test_that('Profile width is 0 when r is 0', {
        alpha <- seq(0, 2*pi, length.out=100)
        theta <- seq(0, 2*pi, length.out=100)
        paras <- expand.grid(alpha, theta)
        paras <- cbind(paras, 0)
        colnames(paras) <- c('theta','alpha','r')
        p <- sapply(1:nrow(paras), function(x) do.call(calcProfileWidth,as.list(paras[x,])))
        expect_that(all(p==0), is_true())
})


# gremDensity() tests

test_that('Density returns error if p is 0', {
#        expect_that(gremDensity(2, 0.1, 0.1, 0, 2, 2), throws_error())
#        expect_that(gremDensity(2, 0.1, 0.1, 10, 0, 2), throws_error())
#        expect_that(gremDensity(2, 0.1, 0.1, 10, 2, 0), throws_error())
#        expect_that(gremDensity(2, 0, 0.1, 10, 2, 2), throws_error())
})

test_that('count of zero gives error', {
        expect_that(gremDensity(0, 0.1, 0.1, 2, 2, 2), throws_error())
})


# gremAbundance() tests





# Compare rem and gREM

test_that('rem and gREM give consistent results.', {
  
  # Gas model.
  tm <- 50 #hours
  v <- 30 # kilometers per day
  r <- 9 # meters
  count <- 11
  
  # grem requires alpha. Set to 2*pi to make equivalent to REM.
  alpha <- 2*pi
  
  # Try gas model here, other REM below.
  theta <- pi / 3
  
  
  grem_result <- gremDensity(count, alpha, theta, r, v, tm)
    
  
  rem_result <- rem(dat = data.frame(1, count, r, theta), tm, v)
    
  expect_true(abs(grem_result - rem_result) < 1e-6)
  
})
