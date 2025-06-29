#  File tests/testthat/test-multilayer-CMB.R in package ergm.multi, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################
data(samplk)
samplkl <- list(samplk1, samplk2, samplk3)

cmb <- function(l) {
  l <- lapply(l, as.matrix)
  n <- length(l)
  m <- Reduce(`+`, l)
  diag(m) <- NA
  sum(lfactorial(m) + lfactorial(n - m) - lfactorial(n), na.rm = TRUE)
}

test_that("3-layer CMBL summary", {
  layer <- summary(Layer(samplkl) ~ CMBL)
  logic <- cmb(samplkl)
  expect_equal(layer, logic, ignore_attr = TRUE)
})

test_that("6-layer hammingL summary with repeated layers", {
  lsel <- sample.int(3, 6, replace = TRUE)
  layer <- summary(Layer(samplkl[lsel]) ~ CMBL)
  logic <- cmb(samplkl[lsel])
  expect_equal(layer, logic, ignore_attr = TRUE)
})
