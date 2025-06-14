#  File tests/testthat/test-multilayer.R in package ergm.multi, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################
# "Correct" transitivity calculators
dediag <- function(m, x=0) {diag(m) <- x; m}

istar <- function(m1,m2, distinct=TRUE){
  sum(colSums(m1)*colSums(m2) - distinct*colSums(m1*m2))
}

ostar <- function(m1,m2, distinct=TRUE){
  sum(rowSums(m1)*rowSums(m2) - distinct*rowSums(m1*m2))
}

twopath <- function(m1,m2, distinct=TRUE){
  sum(colSums(m1)*rowSums(m2) - distinct*colSums(m1*t(m2)))
}

twostar <- function(m1,m2, distinct=TRUE){
  ostar(m1,m2,distinct)
}

library(purrr)
n <- 5
ctrl <- control.simulate.formula(MCMC.burnin=n^2*2, MCMC.interval=n)

test_that("twostarL statistics for directed networks", {
  nw1 <- nw2 <- network.initialize(n, directed = TRUE)
  lnw <- Layer(nw1, nw2)

  sim <- suppressWarnings(simulate(lnw~
                    twostarL(c(~`1`,~`2`), "out",FALSE)+
                    twostarL(c(~`1`,~`2`), "in",FALSE)+
                    twostarL(c(~`1`,~`2`), "path",FALSE)+
                    twostarL(c(~`1`,~`2`), "out", TRUE)+
                    twostarL(c(~`1`,~`2`), "in", TRUE)+
                    twostarL(c(~`1`,~`2`), "path", TRUE),
                  control=ctrl,
                  nsim=200))

  stats <- sapply(sim,
                  function(nw){
                    n <- network.size(nw)/2
                    m <- as.matrix(nw)
                    m1 <- m[seq_len(n),seq_len(n)]
                    m2 <- m[seq_len(n)+n,seq_len(n)+n]
                    
                    c(
                      ostar(m1,m2,FALSE),
                      istar(m1,m2,FALSE),
                      twopath(m1,m2,FALSE),
                      ostar(m1,m2, TRUE),
                      istar(m1,m2, TRUE),
                      twopath(m1,m2, TRUE)
                    )
                  }) %>% t()

  expect_equal(attr(sim,"stats"), stats, ignore_attr=TRUE)
})

test_that("twostarL statistics for undirected networks", {
  nw1 <- nw2 <- network.initialize(n, directed = FALSE)
  lnw <- Layer(nw1, nw2)

  sim <- suppressWarnings(simulate(lnw~
                    twostarL(c(~`1`,~`2`), distinct = FALSE)+
                    twostarL(c(~`1`,~`2`)),
                  control=ctrl,
                  nsim=200))

  stats <- sapply(sim,
                  function(nw){
                    n <- network.size(nw)/2
                    m <- as.matrix(nw)
                    m1 <- m[seq_len(n),seq_len(n)]
                    m2 <- m[seq_len(n)+n,seq_len(n)+n]

                    c(
                      twostar(m1,m2,FALSE),
                      twostar(m1,m2, TRUE)
                    )
                  }) %>% t()

  expect_equal(attr(sim,"stats"), stats, ignore_attr=TRUE)
})


## ### Heterogeneous directedness.
## nw1 <- nw3 <- network.initialize(n,dir=TRUE)
## nw2 <- network.initialize(n,dir=FALSE)
## lnw <- Layer(nw1,nw2,nw3)

## test_that("Multilayer dgw*sp statistics for heterogeneously directed networks 1", {
##   sim <- suppressWarnings(simulate(lnw~
##                     # desp distinct layers
##                     desp(0:n,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
##                     desp(0:n,type="ITP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
##                     desp(0:n,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
##                     desp(0:n,type="ITP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
##                     desp(0:n,type="OSP",L.base=~`1`,Ls.path=c(~`2`,~`3`))+
##                     desp(0:n,type="ISP",L.base=~`1`,Ls.path=c(~`2`,~`3`))+
##                     # ddsp distinct layers
##                     ddsp(0:n,type="OTP",Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
##                     ddsp(0:n,type="ITP",Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
##                     ddsp(0:n,type="OTP",Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
##                     ddsp(0:n,type="ITP",Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
##                     ddsp(0:n,type="OSP",Ls.path=c(~`2`,~`3`))+
##                     ddsp(0:n,type="ISP",Ls.path=c(~`2`,~`3`))+
##                     # dnsp distinct layers
##                     dnsp(0:n,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
##                     dnsp(0:n,type="ITP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
##                     dnsp(0:n,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
##                     dnsp(0:n,type="ITP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
##                     dnsp(0:n,type="OSP",L.base=~`1`,Ls.path=c(~`2`,~`3`))+
##                     dnsp(0:n,type="ISP",L.base=~`1`,Ls.path=c(~`2`,~`3`))+
##                     # desp base and path distinct
##                     desp(0:n,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
##                     desp(0:n,type="ITP",L.base=~`1`,Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
##                     desp(0:n,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
##                     desp(0:n,type="ITP",L.base=~`1`,Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
##                     desp(0:n,type="OSP",L.base=~`1`,Ls.path=c(~`2`,~`2`))+
##                     desp(0:n,type="ISP",L.base=~`1`,Ls.path=c(~`2`,~`2`))+
##                     # ddsp base and path distinct
##                     ddsp(0:n,type="OTP",Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
##                     ddsp(0:n,type="ITP",Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
##                     ddsp(0:n,type="OTP",Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
##                     ddsp(0:n,type="ITP",Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
##                     ddsp(0:n,type="OSP",Ls.path=c(~`2`,~`2`))+
##                     ddsp(0:n,type="ISP",Ls.path=c(~`2`,~`2`))+
##                     # dnsp base and path distinct
##                     dnsp(0:n,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
##                     dnsp(0:n,type="ITP",L.base=~`1`,Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
##                     dnsp(0:n,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
##                     dnsp(0:n,type="ITP",L.base=~`1`,Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
##                     dnsp(0:n,type="OSP",L.base=~`1`,Ls.path=c(~`2`,~`2`))+
##                     dnsp(0:n,type="ISP",L.base=~`1`,Ls.path=c(~`2`,~`2`))+
##                     # desp base and path same
##                     desp(0:n,type="OTP",L.base=~`2`,Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
##                     desp(0:n,type="ITP",L.base=~`2`,Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
##                     desp(0:n,type="OTP",L.base=~`2`,Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
##                     desp(0:n,type="ITP",L.base=~`2`,Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
##                     desp(0:n,type="OSP",L.base=~`2`,Ls.path=c(~`2`,~`2`))+
##                     desp(0:n,type="ISP",L.base=~`2`,Ls.path=c(~`2`,~`2`))+
##                     # ddsp base and path same
##                     ddsp(0:n,type="OTP",Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
##                     ddsp(0:n,type="ITP",Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
##                     ddsp(0:n,type="OTP",Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
##                     ddsp(0:n,type="ITP",Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
##                     ddsp(0:n,type="OSP",Ls.path=c(~`2`,~`2`))+
##                     ddsp(0:n,type="ISP",Ls.path=c(~`2`,~`2`))+
##                     # dnsp base and path same
##                     dnsp(0:n,type="OTP",L.base=~`2`,Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
##                     dnsp(0:n,type="ITP",L.base=~`2`,Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
##                     dnsp(0:n,type="OTP",L.base=~`2`,Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
##                     dnsp(0:n,type="ITP",L.base=~`2`,Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
##                     dnsp(0:n,type="OSP",L.base=~`2`,Ls.path=c(~`2`,~`2`))+
##                     dnsp(0:n,type="ISP",L.base=~`2`,Ls.path=c(~`2`,~`2`))+
##                     # desp distinct base and one layer
##                     desp(0:n,type="OTP",L.base=~`2`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
##                     desp(0:n,type="ITP",L.base=~`2`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
##                     desp(0:n,type="OTP",L.base=~`2`,Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
##                     desp(0:n,type="ITP",L.base=~`2`,Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
##                     desp(0:n,type="OSP",L.base=~`2`,Ls.path=c(~`2`,~`3`))+
##                     desp(0:n,type="ISP",L.base=~`2`,Ls.path=c(~`2`,~`3`))+
##                     # ddsp distinct base and one layer
##                     ddsp(0:n,type="OTP",Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
##                     ddsp(0:n,type="ITP",Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
##                     ddsp(0:n,type="OTP",Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
##                     ddsp(0:n,type="ITP",Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
##                     ddsp(0:n,type="OSP",Ls.path=c(~`2`,~`3`))+
##                     ddsp(0:n,type="ISP",Ls.path=c(~`2`,~`3`))+
##                     # dnsp distinct base and one layer
##                     dnsp(0:n,type="OTP",L.base=~`2`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
##                     dnsp(0:n,type="ITP",L.base=~`2`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
##                     dnsp(0:n,type="OTP",L.base=~`2`,Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
##                     dnsp(0:n,type="ITP",L.base=~`2`,Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
##                     dnsp(0:n,type="OSP",L.base=~`2`,Ls.path=c(~`2`,~`3`))+
##                     dnsp(0:n,type="ISP",L.base=~`2`,Ls.path=c(~`2`,~`3`))
##                  ,
##                   control=ctrl,
##                   nsim=200))

##   stats <- sapply(sim,
##                   function(nw){
##                     n <- network.size(nw)/3
##                     m <- as.matrix(nw)
##                     m1 <- m[seq_len(n),seq_len(n)]
##                     m2 <- m[seq_len(n)+n,seq_len(n)+n]
##                     m2 <- m2 + t(m2)
##                     m3 <- m[seq_len(n)+n*2,seq_len(n)+n*2]
                    
##                     c(
##                       # desp distinct layers
##                       desp(0:n,OTP,L.base=m1,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
##                       desp(0:n,ITP,L.base=m1,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
##                       desp(0:n,OTP,L.base=m1,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
##                       desp(0:n,ITP,L.base=m1,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
##                       desp(0:n,OSP,L.base=m1,Ls.path1=m2,Ls.path2=m3),
##                       desp(0:n,ISP,L.base=m1,Ls.path1=m2,Ls.path2=m3),
##                       # ddsp distinct layers
##                       ddsp(0:n,OTP,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
##                       ddsp(0:n,ITP,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
##                       ddsp(0:n,OTP,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
##                       ddsp(0:n,ITP,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
##                       ddsp(0:n,OSP,Ls.path1=m2,Ls.path2=m3),
##                       ddsp(0:n,ISP,Ls.path1=m2,Ls.path2=m3),
##                       # dnsp distinct layers
##                       dnsp(0:n,OTP,L.base=m1,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
##                       dnsp(0:n,ITP,L.base=m1,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
##                       dnsp(0:n,OTP,L.base=m1,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
##                       dnsp(0:n,ITP,L.base=m1,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
##                       dnsp(0:n,OSP,L.base=m1,Ls.path1=m2,Ls.path2=m3),
##                       dnsp(0:n,ISP,L.base=m1,Ls.path1=m2,Ls.path2=m3),
##                       # desp base and path distinct
##                       desp(0:n,OTP,L.base=m1,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
##                       desp(0:n,ITP,L.base=m1,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
##                       desp(0:n,OTP,L.base=m1,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
##                       desp(0:n,ITP,L.base=m1,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
##                       desp(0:n,OSP,L.base=m1,Ls.path1=m2,Ls.path2=m2),
##                       desp(0:n,ISP,L.base=m1,Ls.path1=m2,Ls.path2=m2),
##                       # ddsp base and path distinct
##                       ddsp(0:n,OTP,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
##                       ddsp(0:n,ITP,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
##                       ddsp(0:n,OTP,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
##                       ddsp(0:n,ITP,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
##                       ddsp(0:n,OSP,Ls.path1=m2,Ls.path2=m2),
##                       ddsp(0:n,ISP,Ls.path1=m2,Ls.path2=m2),
##                       # dnsp base and path distinct
##                       dnsp(0:n,OTP,L.base=m1,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
##                       dnsp(0:n,ITP,L.base=m1,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
##                       dnsp(0:n,OTP,L.base=m1,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
##                       dnsp(0:n,ITP,L.base=m1,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
##                       dnsp(0:n,OSP,L.base=m1,Ls.path1=m2,Ls.path2=m2),
##                       dnsp(0:n,ISP,L.base=m1,Ls.path1=m2,Ls.path2=m2),
##                       # desp base and path same
##                       desp(0:n,OTP,L.base=m2,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
##                       desp(0:n,ITP,L.base=m2,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
##                       desp(0:n,OTP,L.base=m2,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
##                       desp(0:n,ITP,L.base=m2,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
##                       desp(0:n,OSP,L.base=m2,Ls.path1=m2,Ls.path2=m2),
##                       desp(0:n,ISP,L.base=m2,Ls.path1=m2,Ls.path2=m2),
##                       # ddsp base and path same
##                       ddsp(0:n,OTP,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
##                       ddsp(0:n,ITP,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
##                       ddsp(0:n,OTP,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
##                       ddsp(0:n,ITP,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
##                       ddsp(0:n,OSP,Ls.path1=m2,Ls.path2=m2),
##                       ddsp(0:n,ISP,Ls.path1=m2,Ls.path2=m2),
##                       # dnsp base and path same
##                       dnsp(0:n,OTP,L.base=m2,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
##                       dnsp(0:n,ITP,L.base=m2,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
##                       dnsp(0:n,OTP,L.base=m2,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
##                       dnsp(0:n,ITP,L.base=m2,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
##                       dnsp(0:n,OSP,L.base=m2,Ls.path1=m2,Ls.path2=m2),
##                       dnsp(0:n,ISP,L.base=m2,Ls.path1=m2,Ls.path2=m2),
##                       # desp distinct base and one layer
##                       desp(0:n,OTP,L.base=m2,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
##                       desp(0:n,ITP,L.base=m2,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
##                       desp(0:n,OTP,L.base=m2,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
##                       desp(0:n,ITP,L.base=m2,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
##                       desp(0:n,OSP,L.base=m2,Ls.path1=m2,Ls.path2=m3),
##                       desp(0:n,ISP,L.base=m2,Ls.path1=m2,Ls.path2=m3),
##                       # ddsp distinct base and one layer
##                       ddsp(0:n,OTP,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
##                       ddsp(0:n,ITP,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
##                       ddsp(0:n,OTP,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
##                       ddsp(0:n,ITP,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
##                       ddsp(0:n,OSP,Ls.path1=m2,Ls.path2=m3),
##                       ddsp(0:n,ISP,Ls.path1=m2,Ls.path2=m3),
##                       # dnsp distinct base and one layer
##                       dnsp(0:n,OTP,L.base=m2,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
##                       dnsp(0:n,ITP,L.base=m2,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
##                       dnsp(0:n,OTP,L.base=m2,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
##                       dnsp(0:n,ITP,L.base=m2,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
##                       dnsp(0:n,OSP,L.base=m2,Ls.path1=m2,Ls.path2=m3),
##                       dnsp(0:n,ISP,L.base=m2,Ls.path1=m2,Ls.path2=m3)

##                     )
##                   }) %>% t()

##   expect_equivalent(attr(sim,"stats"), stats)
## })


## ### Undirected.
## nw1 <- nw2 <- nw3 <- network.initialize(n,dir=FALSE)
## lnw <- Layer(nw1,nw2,nw3)

## test_that("Multilayer dgw*sp statistics for undirected networks", {
##   sim <- suppressWarnings(simulate(lnw~
##                     # desp distinct layers
##                     desp(0:n,L.base=~`1`,Ls.path=c(~`2`,~`3`))+
##                     # ddsp distinct layers
##                     ddsp(0:n,Ls.path=c(~`2`,~`3`))+
##                     # dnsp distinct layers
##                     dnsp(0:n,L.base=~`1`,Ls.path=c(~`2`,~`3`))+
##                     # desp base and path distinct
##                     desp(0:n,L.base=~`1`,Ls.path=c(~`2`,~`2`))+
##                     # ddsp base and path distinct
##                     ddsp(0:n,Ls.path=c(~`2`,~`2`))+
##                     # dnsp base and path distinct
##                     dnsp(0:n,L.base=~`1`,Ls.path=c(~`2`,~`2`))+
##                     # desp base and path same
##                     desp(0:n,L.base=~`2`,Ls.path=c(~`2`,~`2`))+
##                     # ddsp base and path same
##                     ddsp(0:n,Ls.path=c(~`2`,~`2`))+
##                     # dnsp base and path same
##                     dnsp(0:n,L.base=~`2`,Ls.path=c(~`2`,~`2`))+
##                     # desp distinct base and one layer
##                     desp(0:n,L.base=~`2`,Ls.path=c(~`2`,~`3`))+
##                     # ddsp distinct base and one layer
##                     ddsp(0:n,Ls.path=c(~`2`,~`3`))+
##                     # dnsp distinct base and one layer
##                     dnsp(0:n,L.base=~`2`,Ls.path=c(~`2`,~`3`))
##                  ,
##                   control=ctrl,
##                   nsim=200))

##   stats <- sapply(sim,
##                   function(nw){
##                     n <- network.size(nw)/3
##                     m <- as.matrix(nw)
##                     m1 <- m[seq_len(n),seq_len(n)]
##                     m2 <- m[seq_len(n)+n,seq_len(n)+n]
##                     m3 <- m[seq_len(n)+n*2,seq_len(n)+n*2]
                    
##                     c(
##                       # esp distinct layers
##                       esp(0:n,L.base=m1,Ls.path1=m2,Ls.path2=m3),
##                       # dsp distinct layers
##                       dsp(0:n,Ls.path1=m2,Ls.path2=m3),
##                       # nsp distinct layers
##                       nsp(0:n,L.base=m1,Ls.path1=m2,Ls.path2=m3),
##                       # esp base and path distinct
##                       esp(0:n,L.base=m1,Ls.path1=m2,Ls.path2=m2),
##                       # dsp base and path distinct
##                       dsp(0:n,Ls.path1=m2,Ls.path2=m2),
##                       # nsp base and path distinct
##                       nsp(0:n,L.base=m1,Ls.path1=m2,Ls.path2=m2),
##                       # esp base and path same
##                       esp(0:n,L.base=m2,Ls.path1=m2,Ls.path2=m2),
##                       # dsp base and path same
##                       dsp(0:n,Ls.path1=m2,Ls.path2=m2),
##                       # nsp base and path same
##                       nsp(0:n,L.base=m2,Ls.path1=m2,Ls.path2=m2),
##                       # esp distinct base and one layer
##                       esp(0:n,L.base=m2,Ls.path1=m2,Ls.path2=m3),
##                       # dsp distinct base and one layer
##                       dsp(0:n,Ls.path1=m2,Ls.path2=m3),
##                       # nsp distinct base and one layer
##                       nsp(0:n,L.base=m2,Ls.path1=m2,Ls.path2=m3)
##                     )
##                   }) %>% t()

##   expect_equivalent(attr(sim,"stats"), stats)
## })

nw1 <- network.initialize(20, dir=FALSE)
nw12 <- network.initialize(20, dir=FALSE, bipartite=5)
nw1 %v% "mode" <- rep(1:2,c(5,15))

test_that("Statistics evaluation for heterogeneously bipartite networks", {
  nw1[1:5,1:5] <- 1
  nw12[1:5,6:20] <- 1

  lnw <- Layer(nw1, nw12, .active=list(~mode==1,~TRUE))

  expect_equal(summary(lnw~edges+L(~edges,~`1`)+L(~edges,~`2`)+L(~edges,~(`1`|`2`))+L(~edges,~(`1`&`2`))),
               c(85,10,75,85,0), ignore_attr=TRUE)
})

test_that("Statistics simulation for heterogeneously bipartite networks", {
  nw1[,] <- 0
  nw12[,] <- 0
  lnw <- Layer(nw1, nw12, .active=list(~mode==1,~TRUE))

  simulate(lnw~edges, coef=1000) -> slnw

  expect_equal(summary(slnw~edges+L(~edges,~`1`)+L(~edges,~`2`)+L(~edges,~(`1`|`2`))+L(~edges,~(`1`&`2`))),
               c(85,10,75,85,0), ignore_attr=TRUE)
})

test_that("L() stops if given a non-multilayer object with a sensible error message.",{
  expect_error(ergm_model(nw1 ~ L(~edges)),
               "In term 'L' in package 'ergm\\.multi': The LHS of the model is not a multilayer 'Layer\\(\\)' construct\\.")
})
