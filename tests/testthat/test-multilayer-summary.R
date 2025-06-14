#  File tests/testthat/test-multilayer-summary.R in package ergm.multi, part of
#  the Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################
nw0 <- network.initialize(10, dir=FALSE)
(nw1 <- simulate(nw0~edges, coef=-.5))
(nw2 <- simulate(nw0~edges, coef=+.5))

test_that("basic layer logic summary", {
  (layer <- summary(Layer(A=nw1, B=nw2) ~
                      L(~edges, -1.5~.) +
                      L(~edges, ~A) +
                      L(~edges, ~`2`) +
                      L(~edges, 2~A) +
                      L(~edges, c(~`1`,~B)) +
                      L(~edges, c(1.5~`1`,-.5~B)) +
                      L(~edges, c(~`1`,-.5~B)) +
                      L(~density) +
                      L(~meandeg) +
                      L(~edges, ~A&B) +
                      L(~edges, ~`1`||`2`) +
                      L(~edges, ~(!A)&`2`) +
                      L(~edges, ~(`1`!=`2`)) +
                      L(~edges, ~xor(`1`,B))
                    ))
  (logic <- c((summary(nw1~edges)+summary(nw2~edges))*-1.5,
              summary(nw1~edges),
              summary(nw2~edges),
              summary(nw1~edges)*2,
              summary(nw1~edges)+summary(nw2~edges),
              summary(nw1~edges)*1.5+summary(nw2~edges)*-.5,
              summary(nw1~edges)+summary(nw2~edges)*-.5,
              summary(nw1~density)+summary(nw2~density),
              summary(nw1~meandeg)+summary(nw2~meandeg),
              summary((nw1&nw2)~edges),
              summary((nw1|nw2)~edges),
              summary(((!nw1)&nw2)~edges),
              summary(((nw1&!nw2)|(!nw1&nw2))~edges),
              summary(((nw1&!nw2)|(!nw1&nw2))~edges)
              ))

  expect_equal(layer, logic, ignore_attr=TRUE)
})
