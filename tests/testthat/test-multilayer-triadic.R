#  File tests/testthat/test-multilayer-triadic.R in package ergm.multi, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2022 Statnet Commons
################################################################################
data(florentine)

test_that("multilayer triadic summary: dsp", {
  (layer <- summary(Layer(m=flomarriage, b=flobusiness)~ddspL(0:10,Ls.path=c(~b,~b))))
  (logic <- summary(flobusiness~dsp(0:10)))
  expect_equal(layer, logic, ignore_attr=TRUE)
})

test_that("multilayer triadic summary: esp", {
  (layer <- summary(Layer(m=flomarriage, b=flobusiness)~despL(0:10,Ls.path=c(~b,~b),L.base=~b)))
  (logic <- summary(flobusiness~esp(0:10)))
  expect_equal(layer, logic, ignore_attr=TRUE)
})

test_that("multilayer triadic summary: nsp", {
  (layer <- summary(Layer(m=flomarriage, b=flobusiness)~dnspL(0:10,Ls.path=c(~b,~b),L.base=~b)))
  (logic <- summary(flobusiness~nsp(0:10)))
  expect_equal(layer, logic, ignore_attr=TRUE)
})
