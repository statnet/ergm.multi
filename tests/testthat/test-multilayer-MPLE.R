#  File tests/testthat/test-multilayer-MPLE.R in package ergm.multi, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2024 Statnet Commons
################################################################################
logit <- function(p) log(p/(1-p))

nw0 <- network.initialize(10, dir=FALSE)
(nw1 <- simulate(nw0~edges, coef=-.5))
(nw2 <- simulate(nw0~edges, coef=+.5))

test_that("multilayer MPLE for two independent models", {
  (layer <- coef(ergm(Layer(nw1, nw2) ~ L(~edges, ~`1`) + L(~edges, ~`2`))))
  (logic <- c(coef(ergm(nw1~edges)), coef(ergm(nw2~edges))))
  expect_equal(layer, logic, ignore_attr=TRUE, tolerance=1e-7)
})

test_that("multilayer MPLE for two independent models pooled", {
  (layer <- coef(ergm(Layer(nw1, nw2) ~ L(~edges, c(~`1`,~`2`)))))
  (logic <- logit((network.edgecount(nw1)+network.edgecount(nw2))/(network.dyadcount(nw1)+network.dyadcount(nw2))))
  expect_equal(layer, logic, ignore_attr=TRUE)
})
