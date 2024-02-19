#  File tests/testthat/test-multilayer-heterogeneous-directedness.R in package ergm.multi, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2024 Statnet Commons
################################################################################
nwu <- network.initialize(2, dir=FALSE)
nwu[1,2] <- 1
nwd <- network.initialize(2, dir=TRUE)
nwd[2,1] <- 1

test_that("multilayer heterogeneous directedness summary", {
  lnw <- Layer(nwu, nwd)
  expect_equal(summary(lnw~L(~edges,~`1`)+L(~edges,~`2`)+CMBL), c(2,1,-sum(lchoose(2,as.matrix(nwu)+as.matrix(nwd)))), ignore_attr=TRUE)
})
