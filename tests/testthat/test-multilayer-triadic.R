data(florentine)

test_that("multilayer triadic summary: dsp", {
  local_edition(3)
  (layer <- summary(Layer(m=flomarriage, b=flobusiness)~ddspL(0:10,Ls.path=c(~b,~b))))
  (logic <- summary(flobusiness~dsp(0:10)))
  expect_equal(layer, logic, ignore_attr=TRUE)
})

test_that("multilayer triadic summary: esp", {
  local_edition(3)
  (layer <- summary(Layer(m=flomarriage, b=flobusiness)~despL(0:10,Ls.path=c(~b,~b),L.base=~b)))
  (logic <- summary(flobusiness~esp(0:10)))
  expect_equal(layer, logic, ignore_attr=TRUE)
})

test_that("multilayer triadic summary: nsp", {
  local_edition(3)
  (layer <- summary(Layer(m=flomarriage, b=flobusiness)~dnspL(0:10,Ls.path=c(~b,~b),L.base=~b)))
  (logic <- summary(flobusiness~nsp(0:10)))
  expect_equal(layer, logic, ignore_attr=TRUE)
})
