#  File tests/testthat/test-multilayer-heterogeneous-directedness.R in package
#  ergm.multi, part of the Statnet suite of packages for network analysis,
#  https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################
test_that("multilayer heterogeneous directedness summary", {
  nwu <- network.initialize(2, directed = FALSE)
  nwu[1,2] <- 1
  nwd <- network.initialize(2, directed = TRUE)
  nwd[2,1] <- 1

  lnw <- Layer(nwu, nwd)
  expect_equal(summary(lnw~L(~edges,~`1`)+L(~edges,~`2`)+CMBL), c(2,1,-sum(lchoose(2,as.matrix(nwu)+as.matrix(nwd)))), ignore_attr=TRUE)
})

test_that("multilayer heterogeneous layers messages", {
  nw1 <- network.initialize(20, dir=FALSE)
  nw12 <- network.initialize(20, dir=FALSE, bipartite=5)
  nw1 %v% "mode" <- rep(1:2,c(5,15))
  nw1 %n% "nattr" <- "abc"

  expect_silent(Layer(nw1, nw12, .active=list(~mode==1, ~TRUE)))

  expect_message(expect_message(
    Layer(nw12, nw1),
    "Vertex attribute\\(s\\) 'mode' are not found in the first layer.*"),
    "Network attribute\\(s\\) 'nattr' are not found in the first layer.*"
  )

  nw12 %v% "mode" <- rep(1:2,c(15,5))
  nw12 %n% "nattr" <- "def"
  expect_message(expect_message(
    Layer(nw1, nw12, .active=list(~mode==1, ~TRUE)),
    "Vertex attribute\\(s\\) 'mode' have values different from those in the first layer.*"),
    "Network attribute\\(s\\) 'nattr' have values different from those in the first layer.*"
  )
})
