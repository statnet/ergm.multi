#  File tests/testthat/test-multilayer-CMB.R in package ergm.multi, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2024 Statnet Commons
################################################################################
data(samplk)

test_that("3-layer CMBL summary", {
  (layer <- summary(Layer(samplk1, samplk2, samplk3)~CMBL))
  m1 <- as.matrix(samplk1)
  m2 <- as.matrix(samplk2)
  m3 <- as.matrix(samplk3)
  msum <- m1 + m2 + m3
  diag(msum) <- NA
  (logic <- sum(lfactorial(msum) + lfactorial(3-msum) - lfactorial(3), na.rm=TRUE))
  expect_equal(layer, logic, ignore_attr=TRUE)
})
