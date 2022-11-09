#  File tests/testthat/test-multilayer-MLE.R in package ergm.multi, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2022 Statnet Commons
################################################################################
nw0 <- network.initialize(3, dir=FALSE)
nw1 <- nw0
nw1[1,2] <- 1
nw1[2,3] <- 1

nw2 <- nw0
nw2[1,2] <- 1
nw2[1,3] <- 1

layer_and_MLE <- function(nw1, nw2){
  nd <- network.dyadcount(nw1)
  ne <- summary(nw1&nw2~edges)
  log(3*ne/nd)-log(1-ne/nd)
}

layer_and_Info <- function(nw1, nw2, theta=layer_and_MLE(nw1, nw2)){
  network.dyadcount(nw1)*(3*exp(theta))/(3+exp(theta))^2
}

test_that("layer logic estimation for a single AND layer", {
for(s in 1:1000){
  set.seed(s)
  layer <- ergm(Layer(nw1,nw2)~L(~edges, ~`1`&`2`))
  logic.coef <- layer_and_MLE(nw1,nw2)
  logic.info <- layer_and_Info(nw1,nw2, coef(layer))

  expect_equal(sqrt(c(vcov(layer, sources="model"))), sqrt(1/logic.info), ignore_attr=TRUE, tolerance=0.3)
  expect_lt(abs(layer_and_MLE(nw1,nw2)-coef(layer))/sqrt(vcov(layer, sources="estimation")), 4)
}
})
