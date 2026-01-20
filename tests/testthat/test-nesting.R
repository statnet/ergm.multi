data(samplk)
data(florentine)

nw <- Networks(Layer(samplk1, samplk2, samplk3),
               Layer(sampdlk1, sampdlk2))

f <- nw ~
  N(~L(~edges, ~`1`) + L(~edges, ~`2`)) +
  N(~L(~edges, ~`3`), subset = 1)

test_that("summary of two multilayer networks with heterogeneous layer counts", {
  expect_equal(unname(summary(f)),
               c(network.edgecount(samplk1) + network.edgecount(sampdlk1),
                 network.edgecount(samplk2) + network.edgecount(sampdlk2),
                 network.edgecount(samplk3)))
})


test_that("simulation of two multilayer networks with heterogeneous layer counts", {
  expect_silent(nw1 <- simulate(f, coef = c(0, 0, 0)))
  expect_error(unLayer(nw1), ".*not a multilayer network at the top level.*")
  unw <- unNetworks(nw1)
  expect_length(unw, 2L)
  expect_identical(unw[[1]]  %ergmlhs% "constraints", statnet.common::base_env(~blockdiag(".LayerID", noncontig = "split")))
 
  uunw1 <- unLayer(unw[[1]])
  uunw2 <- unLayer(unw[[2]])
  expect_length(uunw1, 3L)
  expect_length(uunw2, 2L)

  expect_identical(uunw1[[1]]  %ergmlhs% "constraints", NULL)
})
