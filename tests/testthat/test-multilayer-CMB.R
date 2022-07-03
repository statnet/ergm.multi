data(samplk)

test_that("3-layer CMBL summary", {
  local_edition(3)
  (layer <- summary(Layer(samplk1, samplk2, samplk3)~CMBL))
  m1 <- as.matrix(samplk1)
  m2 <- as.matrix(samplk2)
  m3 <- as.matrix(samplk3)
  msum <- m1 + m2 + m3
  diag(msum) <- NA
  (logic <- sum(lfactorial(msum) + lfactorial(3-msum) - lfactorial(3), na.rm=TRUE))
  expect_equal(layer, logic, ignore_attr=TRUE)
})
