data(samplk)

test_that("subnet cache is only used if called through .split_constr_network()", {
  samplks <- combine_networks(list(samplk1, samplk2, samplk3), subnet.cache=TRUE)

  expect_equal(s <- unname(summary(samplks~edges)), sum(sapply(list(samplk1, samplk2, samplk3), network.edgecount)))
  expect_equal(sum(summary(uncombine_network(samplks)~edges)), s)
  expect_equal(sum(sapply(.split_constr_network(samplks), network.edgecount)), 0)
})
