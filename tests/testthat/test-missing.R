library(ergm.multi)

data(samplk)

samplk1[1,] <- NA
samplk2[,2] <- NA

nwl <- list(samplk1, samplk2, samplk3)

truth <- qlogis(sum(sapply(nwl, network.edgecount, na.omit = TRUE)) /
                sum(sapply(nwl, network.dyadcount, na.omit = TRUE)))

test_that("Multilayer network with missing data", {
  expect_equal(coef(ergm(Layer(nwl) ~ edges)), truth, ignore_attr = TRUE)
})

test_that("Network sample with missing data", {
  expect_equal(coef(ergm(Networks(nwl) ~ edges)), truth, ignore_attr = TRUE)
})
