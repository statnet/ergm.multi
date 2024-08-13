test_that("Layers could be nested within Networks, but unpacking different combiners in the wrong order is an error", {
  data(samplk)
  Lsamp <- Layer(samplk1, samplk2, samplk3)
  Lsamps <- rep(list(Lsamp), 3)
  for(i in seq_along(Lsamps)) set.network.attribute(Lsamps[[i]], "index", i)
  NLsamps <- Networks(Lsamps)

  expect_equal(summary(NLsamps ~ N(~L(~edges, ~`1`) + L(~edges, ~`2`) + L(~edges, ~`3`), ~0+factor(index))),
               c("N(factor(index)1)~L(1)~edges" = 55, "N(factor(index)2)~L(1)~edges" = 55, "N(factor(index)3)~L(1)~edges" = 55,
                 "N(factor(index)1)~L(2)~edges" = 57, "N(factor(index)2)~L(2)~edges" = 57, "N(factor(index)3)~L(2)~edges" = 57,
                 "N(factor(index)1)~L(3)~edges" = 56, "N(factor(index)2)~L(3)~edges" = 56, "N(factor(index)3)~L(3)~edges" = 56))

  expect_error(summary(NLsamps ~ L(~edges)),
               ".*The LHS was \\(at the top level\\) created by 'Networks\\(\\)' but the term is trying to extract its layers \\(created by 'Layer\\(\\)'\\)\\..*")
})

