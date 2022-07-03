data(samplk)

test_that("multilayer mutualL term", {
  local_edition(3)
  (layer <- summary(Layer(samplk1, samplk2)~mutualL(Ls=c(~`1`,~`2`))))
  m1 <- as.matrix(samplk1)
  m2 <- as.matrix(samplk2)
  (logic <- (sum(m1*t(m2)+m2*t(m1))/2))
  expect_equal(layer, logic, ignore_attr=TRUE)
})

test_that("multilayer mutualL term with additional layer logic", {
  local_edition(3)
  (layer <- summary(Layer(samplk1, samplk2)~mutualL(Ls=c(~`1`,~`2`&`1`))))
  m1 <- as.matrix(samplk1)
  m2 <- as.matrix(samplk2) * as.matrix(samplk1)
  (logic <- (sum(m1*t(m2)+m2*t(m1))/2))
  expect_equal(layer, logic, ignore_attr=TRUE)
})
