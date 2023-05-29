#  File tests/testthat/test-multilayer-mutuality.R in package ergm.multi, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2023 Statnet Commons
################################################################################
data(samplk)

test_that("multilayer mutualL term", {
  (layer <- summary(Layer(samplk1, samplk2)~mutualL(Ls=c(~`1`,~`2`))))
  m1 <- as.matrix(samplk1)
  m2 <- as.matrix(samplk2)
  (logic <- (sum(m1*t(m2)+m2*t(m1))/2))
  expect_equal(layer, logic, ignore_attr=TRUE)
})

test_that("multilayer mutualL term with additional layer logic", {
  (layer <- summary(Layer(samplk1, samplk2)~mutualL(Ls=c(~`1`,~`2`&`1`))))
  m1 <- as.matrix(samplk1)
  m2 <- as.matrix(samplk2) * as.matrix(samplk1)
  (logic <- (sum(m1*t(m2)+m2*t(m1))/2))
  expect_equal(layer, logic, ignore_attr=TRUE)
})
