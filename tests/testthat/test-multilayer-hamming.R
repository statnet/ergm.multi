#  File tests/testthat/test-multilayer-hamming.R in package
#  ergm.multi, part of the Statnet suite of packages for network
#  analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################
data(samplk)
samplkl <- list(samplk1, samplk2, samplk3)

ham <- function(l) {
  l <- lapply(l, as.matrix)
  sapply(seq_along(l), function(i)
    sapply(seq_len(i - 1L), function(j)
      sum(l[[i]] != l[[j]]))) |>
    unlist() |>
    sum()
}

test_that("3-layer hammingL summary", {
  layer <- summary(Layer(samplkl) ~ hammingL(~.))
  logic <- ham(samplkl)
  expect_equal(layer, logic, ignore_attr = TRUE)
})

test_that("6-layer hammingL summary with repeated layers", {
  lsel <- sample.int(3, 6, replace = TRUE)
  layer <- summary(Layer(samplkl[lsel]) ~ hammingL)
  logic <- ham(samplkl[lsel])
  expect_equal(layer, logic, ignore_attr = TRUE)
})


ent <- function(l) {
  l <- lapply(l, as.matrix)
  sapply(seq_along(l), function(i)
    sapply(seq_len(i - 1L), function(j)
      sum(l[[i]] & l[[j]]))) |>
    unlist() |>
    sum()
}

test_that("3-layer entrainmentL summary", {
  layer <- summary(Layer(samplkl) ~ entrainmentL(~.))
  logic <- ent(samplkl)
  expect_equal(layer, logic, ignore_attr = TRUE)
})

test_that("6-layer entrainmentL summary with repeated layers", {
  lsel <- sample.int(3, 6, replace = TRUE)
  layer <- summary(Layer(samplkl[lsel]) ~ entrainmentL)
  logic <- ent(samplkl[lsel])
  expect_equal(layer, logic, ignore_attr = TRUE)
})
