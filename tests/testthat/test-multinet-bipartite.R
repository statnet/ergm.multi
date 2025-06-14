#  File tests/testthat/test-multinet-bipartite.R in package ergm.multi, part of
#  the Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################

n <- 10
b1 <- 4
b2 <- n - b1
nw0 <- network.initialize(n, directed = FALSE, bipartite = b1)
nw0 %v% "a" <- runif(n) < 0.5
nwl <- c(unclass(simulate(nw0~edges, coef = -1, nsim = 3)))


logit <- function(p) log(p/(1-p))
times <- 1:3

for(i in times) {
  nwl[[i]] %e% "w" <- floor(runif(network.edgecount(nwl[[i]]), 0, 4))
}

nws <- Networks(nwl)

test_that("Bipartite N() summary with an LM and noncompacted stats", {
  summ.N <- summary(nws~N(~edges+nodematch("a"), ~1+.NetworkID), term.options=list(N.compact_stats=FALSE))
  summ.l <- unlist(lapply(nwl, function(nw) summary(~edges+nodematch("a"), basis = nw)))
  names(summ.l) <- paste0("N#", rep(1:3, each=2), "~", names(summ.l))
  expect_equal(summ.N, summ.l)
})


test_that("Bipartite N() summary with offset and compacted stats", {
  summ.N <- summary(nws~N(~edges, offset=~.NetworkID))
  summ.l <- sapply(nwl, function(nw) summary(~edges, basis = nw))
  summ.l <- c(sum(summ.l), sum(summ.l*1:3))
  expect_equal(summ.N, summ.l, ignore_attr=TRUE)
})

test_that("Bipartite Valued N() summary with an LM and noncompacted stats", {
  summ.N <- summary(nws~N(~sum+nodematch("a", form="sum"), ~1+.NetworkID), term.options=list(N.compact_stats=FALSE), response="w")
  summ.l <- unlist(lapply(nwl, function(nw) summary(~sum+nodematch("a", form="sum"), basis = nw, response="w")))
  expect_equal(summ.N, summ.l, ignore_attr=TRUE)
})


test_that("Bipartite Valued N() summary with offset and compacted stats", {
  summ.N <- summary(nws~N(~sum, offset=~.NetworkID), response="w")
  summ.l <- sapply(nwl, function(nw) summary(~sum, basis = nw, response="w"))
  summ.l <- c(sum(summ.l), sum(summ.l*1:3))
  expect_equal(summ.N, summ.l, ignore_attr=TRUE)
})
