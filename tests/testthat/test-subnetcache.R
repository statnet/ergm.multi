#  File tests/testthat/test-subnetcache.R in package ergm.multi, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################
data(samplk)

test_that("combine_networks roundtrip", {
  samplkl <- list(samplk1, samplk2, samplk3)
  samplks <- combine_networks(samplkl)
  samplkl2 <- uncombine_network(samplks, populate = TRUE)
  samplkl0 <- uncombine_network(samplks, populate = FALSE)

  expect_equal(lapply(samplkl2, as.edgelist, attrname = "score", output = "tibble"),
               lapply(samplkl, as.edgelist, attrname = "score", output = "tibble"))

  expect_equal(sum(sapply(samplkl0, network.edgecount)), 0)
  expect_equal(sum(sapply(subnetwork_templates(samplks), network.edgecount)), 0)
})
