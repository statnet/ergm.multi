library(ergm.multi)

data(samplk)

sampdlk1 %n% "time" <- samplk1 %n% "time" <- 1
sampdlk2 %n% "time" <- samplk2 %n% "time" <- 2
sampdlk3 %n% "time" <- samplk3 %n% "time" <- 3

samplk1 %n% "valence" <- samplk2 %n% "valence" <- samplk3 %n% "valence" <- +1
sampdlk1 %n% "valence" <- sampdlk2 %n% "valence" <- sampdlk3 %n% "valence" <- -1

lnw <- Layer(l1 = samplk1, l2 = samplk2, l3 = samplk3,
             d1 = sampdlk1, d2 = sampdlk2, d3 = sampdlk3)


test_that("selecting layers via functions", {
  s <- summary(lnw ~ L(~edges, function(lattr, ...) lattr$time == 2) + L(~edges, c(~l2, ~d2)) +
                 L(~edges, function(lattr, ...) lattr$valence > 0) + L(~edges, 1:3) +
                 L(~edges, function(lattr, ...) rev(lattr$.LayerName)) + L(~edges, ~.))
  expect_equal(s[seq(1, 5, 2)], s[seq(2, 6, 2)], ignore_attr = TRUE)
})
