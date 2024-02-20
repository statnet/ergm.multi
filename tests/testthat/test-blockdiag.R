test_that("Block-diagonal proposal declination if there are off-block-diagonal edges", {
  y <- network.initialize(20, directed=FALSE)
  y %v% "b" <- rep(1:2, each=10)
  y %ergmlhs% "constraints" <- ~ blockdiag("b")
  
  expect_message(simulate(y~edges, coef=0, verbose=TRUE, control=control.simulate.formula(MCMC.prop=~sparse, MCMC.burnin=1)), ".*ergm.multi:MH_blockdiagTNT.*")

  y[1,11] <- 1
  expect_message(simulate(y~edges, coef=0, verbose=TRUE, control=control.simulate.formula(MCMC.prop=~sparse, MCMC.burnin=1)), ".*ergm:MH_TNT.*")
})

