#  File tests/testthat/test-multilayer-direction.R in package ergm.multi, part
#  of the Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################
data(samplk)
m1 <- as.matrix(samplk1)
m2 <- as.matrix(samplk2)
m3 <- as.matrix(samplk3)
msum <- m1 + m2 + m3
diag(msum) <- NA

test_that("multilayer direction reversal summary", {
  (layer <- summary(Layer(samplk1, samplk2, samplk3)~
                      L(~edges, ~`1`&t(`1`))+
                      L(~edges, ~t(t(`1`)&`2`))+
                      L(~edges, ~t(`1`&t(`2`)))+
                      L(~edges, ~`1`&t(`2`))+
                      L(~edges, ~t(`1`)&`2`)+
                      L(~edges, ~t(`1`)|`1`)))
  (logic <- c(sum(m1*t(m1)), sum(m1*t(m2)), sum(t(m1)*m2), sum(m1*t(m2)), sum(t(m1)*m2), sum(t(m1)|m1)))
  expect_equal(layer, logic, ignore_attr=TRUE)
})

test_that("multilayer direction reversal with CMBL summary", {
  (layer <- summary(Layer(samplk1, samplk2, samplk3)~CMBL(c(~`1`,~t(`2`),~`3`))))
  msum2r <- m1 + t(m2) + m3
  diag(msum) <- NA
  (logic <- sum(lfactorial(msum2r) + lfactorial(3-msum2r) - lfactorial(3), na.rm=TRUE))
  expect_equal(layer, logic, ignore_attr=TRUE)
})
