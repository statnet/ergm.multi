nwu <- network.initialize(2, dir=FALSE)
nwu[1,2] <- 1
nwd <- network.initialize(2, dir=TRUE)
nwd[2,1] <- 1

test_that("multilayer heterogeneous directedness summary", {
  local_edition(3)
  lnw <- Layer(nwu, nwd)
  expect_equal(summary(lnw~L(~edges,~`1`)+L(~edges,~`2`)+CMBL), c(2,1,-sum(lchoose(2,as.matrix(nwu)+as.matrix(nwd)))), ignore_attr=TRUE)
})
