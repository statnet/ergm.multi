#  File tests/testthat/test-subnetcache.R in package ergm.multi, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2024 Statnet Commons
################################################################################
data(samplk)

test_that("subnet cache is only used if called through subnetwork_templates()", {
  samplks <- combine_networks(list(samplk1, samplk2, samplk3), subnet.cache=TRUE)

  expect_equal(s <- unname(summary(samplks~edges)), sum(sapply(list(samplk1, samplk2, samplk3), network.edgecount)))
  expect_equal(sum(summary(uncombine_network(samplks)~edges)), s)
  expect_equal(sum(sapply(subnetwork_templates(samplks), network.edgecount)), 0)
})
