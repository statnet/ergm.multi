library(testthat)
context("test-term-dgwesp-ml.R")

# "Correct" transitivity calculators

dediag <- function(m, x=0) {diag(m) <- x; m}

UTP <- function(m1, m2){
  dediag(m1%*%m2+m2%*%m1-(m1*m2)%*%(m1*m2))
}

OTP <- function(m1, m2, in.order=TRUE){
  if(in.order) dediag(m1%*%m2)
  else dediag(m1%*%m2+m2%*%m1-(m1*m2)%*%(m1*m2))
}

ITP <- function(m1, m2, in.order=TRUE){
  if(in.order) dediag(t(m1%*%m2))
  else dediag(t(m1%*%m2+m2%*%m1-(m1*m2)%*%(m1*m2)))
}

OSP <- function(m1, m2){
  dediag(m1%*%t(m2)+m2%*%t(m1)-(m1*m2)%*%t(m1*m2))
}

ISP <- function(m1, m2){
  dediag(t(m1)%*%m2+t(m2)%*%m1-t(m1*m2)%*%(m1*m2))
}

espL <- function(x, L.base, Ls.path1, Ls.path2=Ls.path1, ...){
  TP <- UTP(Ls.path1, Ls.path2, ...)
  L.base[L.base==0] <- NA # I.e., those with base=0 don't count at all.
  espL <- dediag(L.base*TP, NA)[upper.tri(TP)]
  tabulate(match(espL, x),length(x))
}

dspL <- function(x, Ls.path1, Ls.path2=Ls.path1, ...){
  TP <- UTP(Ls.path1, Ls.path2, ...)
  dspL <- dediag(TP, NA)[upper.tri(TP)]
  tabulate(match(dspL, x),length(x))
}

nspL <- function(x, L.base, Ls.path1, Ls.path2=Ls.path1, ...){
  TP <- UTP(Ls.path1, Ls.path2, ...)
  L.base[L.base==1] <- NA # I.e., those with base=1 don't count at all.
  nspL <- dediag((1-L.base)*TP, NA)[upper.tri(TP)]
  tabulate(match(nspL, x),length(x))
}


despL <- function(x, type, L.base, Ls.path1, Ls.path2=Ls.path1, ...){
  TP <- type(Ls.path1, Ls.path2, ...)
  L.base[L.base==0] <- NA # I.e., those with base=0 don't count at all.
  espL <- dediag(L.base*TP, NA)
  tabulate(match(espL, x),length(x))
}

ddspL <- function(x, type, Ls.path1, Ls.path2=Ls.path1, ...){
  TP <- type(Ls.path1, Ls.path2, ...)
  dspL <- dediag(TP, NA)
  tabulate(match(dspL, x),length(x))
}

dnspL <- function(x, type, L.base, Ls.path1, Ls.path2=Ls.path1, ...){
  TP <- type(Ls.path1, Ls.path2, ...)
  L.base[L.base==1] <- NA # I.e., those with base=1 don't count at all.
  nspL <- dediag((1-L.base)*TP, NA)
  tabulate(match(nspL, x),length(x))
}

GW <- function(decay, n){
  i <- 1:n
  exp(decay)*(1-(1-exp(-decay))^i)
}

dgwespL <- function(decay, n, type, L.base, Ls.path1, Ls.path2=Ls.path1, ...){
  w <- GW(decay,n)
  sp <- despL(1:n, type, L.base, Ls.path1, Ls.path2, ...)
  sum(w*sp)
}

dgwdspL <- function(decay, n, type, Ls.path1, Ls.path2=Ls.path1, ...){
  w <- GW(decay,n)
  sp <- ddspL(1:n, type, Ls.path1, Ls.path2, ...)
  sum(w*sp)
}

dgwnspL <- function(decay, n, type, L.base, Ls.path1, Ls.path2=Ls.path1, ...){
  w <- GW(decay,n)
  sp <- dnspL(1:n, type, L.base, Ls.path1, Ls.path2, ...)
  sum(w*sp)
}

gwespL <- function(decay, n, L.base, Ls.path1, Ls.path2=Ls.path1, ...){
  w <- GW(decay,n)
  sp <- espL(1:n, L.base, Ls.path1, Ls.path2, ...)
  sum(w*sp)
}

gwdspL <- function(decay, n, Ls.path1, Ls.path2=Ls.path1, ...){
  w <- GW(decay,n)
  sp <- dspL(1:n, Ls.path1, Ls.path2, ...)
  sum(w*sp)
}

gwnspL <- function(decay, n, L.base, Ls.path1, Ls.path2=Ls.path1, ...){
  w <- GW(decay,n)
  sp <- nspL(1:n, L.base, Ls.path1, Ls.path2, ...)
  sum(w*sp)
}

library(ergm)
library(purrr)
n <- 5

#### Some code useful for debugging.

## # Construct a transitive triad base-first and then remove the base:
## testseq1 <- list(matrix(c(1,2),ncol=2,byrow=TRUE), # Base
##                  matrix(c(1+n,3+n),ncol=2,byrow=TRUE), # Segment 1
##                  matrix(c(3+2*n,2+2*n),ncol=2,byrow=TRUE), # Segment 2
##                  matrix(c(1,2),ncol=2,byrow=TRUE)) # -Base

## # Construct a transitive triad base-last and then remove the base:
## testseq2 <- list(matrix(c(1+n,3+n),ncol=2,byrow=TRUE), # Segment 1
##                  matrix(c(3+2*n,2+2*n),ncol=2,byrow=TRUE), # Segment 2
##                  matrix(c(1,2),ncol=2,byrow=TRUE), # Base
##                  matrix(c(1,2),ncol=2,byrow=TRUE)) # -Base


## testseq3 <- list(
##   matrix(c(1,2),ncol=2,byrow=TRUE), # Base in L1
##   matrix(c(1+2*n,3+2*n),ncol=2,byrow=TRUE), # Segment 1 in L3
##   matrix(c(3+1*n,2+1*n),ncol=2,byrow=TRUE), # Segment 2 in L2
##   matrix(c(1,2),ncol=2,byrow=TRUE)) # -Base in L1

## # Construct a transitive triad base-last and then remove the base:
## testseq4 <- list(
##   matrix(c(1+2*n,3+2*n),ncol=2,byrow=TRUE), # Segment 1 in L3
##   matrix(c(3+1*n,2+1*n),ncol=2,byrow=TRUE), # Segment 2 in L2
##   matrix(c(1,2),ncol=2,byrow=TRUE), # Base in L1
##   matrix(c(1+2*n,3+2*n),ncol=2,byrow=TRUE) # -Segment 1 in L3
##   )

## ergm.godfather(lnw~edges+
##                  ## despL(1,type="OTP",L.base=~`1`,Ls.path=c(~`3`,~`2`),L.in_order=TRUE)+
##                  ## despL(1,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
##                  despL(1,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=FALSE),
##                changes = testseq4,stats.start=TRUE)

## summary(lnw~despL(1:18,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE))
## ergm.godfather(lnw~edges+despL(1:n,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE),
##                changes = testseq1,stats.start=TRUE)
## ergm.godfather(lnw~edges+despL(1:n,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE),
##                changes = testseq2,stats.start=TRUE)

ctrl <- control.simulate.formula(MCMC.burnin=1, MCMC.interval=1)
decay <- runif(1,0,1)

for(cache.sp in c(FALSE,TRUE)){
  options(ergm.term=list(cache.sp=cache.sp))
  sptxt <- if(cache.sp) "with shared partner caching" else "without shared partner caching"


### Directed.
nw1 <- nw2 <- nw3 <- network.initialize(n,dir=TRUE)
lnw <- Layer(nw1,nw2,nw3)

test_that(paste("Multilayer dgw*sp statistics for homogeneously directed networks",sptxt), {
  sim <- suppressWarnings(simulate(lnw~
                    # despL distinct layers
                    despL(0:n,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    despL(0:n,type="ITP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    despL(0:n,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    despL(0:n,type="ITP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    despL(0:n,type="OSP",L.base=~`1`,Ls.path=c(~`2`,~`3`))+
                    despL(0:n,type="ISP",L.base=~`1`,Ls.path=c(~`2`,~`3`))+
                    # ddspL distinct layers
                    ddspL(0:n,type="OTP",Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    ddspL(0:n,type="ITP",Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    ddspL(0:n,type="OTP",Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    ddspL(0:n,type="ITP",Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    ddspL(0:n,type="OSP",Ls.path=c(~`2`,~`3`))+
                    ddspL(0:n,type="ISP",Ls.path=c(~`2`,~`3`))+
                    # dnspL distinct layers
                    dnspL(0:n,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    dnspL(0:n,type="ITP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    dnspL(0:n,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    dnspL(0:n,type="ITP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    dnspL(0:n,type="OSP",L.base=~`1`,Ls.path=c(~`2`,~`3`))+
                    dnspL(0:n,type="ISP",L.base=~`1`,Ls.path=c(~`2`,~`3`))+
                    # despL base and path distinct
                    despL(0:n,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    despL(0:n,type="ITP",L.base=~`1`,Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    despL(0:n,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    despL(0:n,type="ITP",L.base=~`1`,Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    despL(0:n,type="OSP",L.base=~`1`,Ls.path=c(~`2`,~`2`))+
                    despL(0:n,type="ISP",L.base=~`1`,Ls.path=c(~`2`,~`2`))+
                    # ddspL base and path distinct
                    ddspL(0:n,type="OTP",Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    ddspL(0:n,type="ITP",Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    ddspL(0:n,type="OTP",Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    ddspL(0:n,type="ITP",Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    ddspL(0:n,type="OSP",Ls.path=c(~`2`,~`2`))+
                    ddspL(0:n,type="ISP",Ls.path=c(~`2`,~`2`))+
                    # dnspL base and path distinct
                    dnspL(0:n,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    dnspL(0:n,type="ITP",L.base=~`1`,Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    dnspL(0:n,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    dnspL(0:n,type="ITP",L.base=~`1`,Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    dnspL(0:n,type="OSP",L.base=~`1`,Ls.path=c(~`2`,~`2`))+
                    dnspL(0:n,type="ISP",L.base=~`1`,Ls.path=c(~`2`,~`2`))+
                    # despL base and path same
                    despL(0:n,type="OTP",L.base=~`2`,Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    despL(0:n,type="ITP",L.base=~`2`,Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    despL(0:n,type="OTP",L.base=~`2`,Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    despL(0:n,type="ITP",L.base=~`2`,Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    despL(0:n,type="OSP",L.base=~`2`,Ls.path=c(~`2`,~`2`))+
                    despL(0:n,type="ISP",L.base=~`2`,Ls.path=c(~`2`,~`2`))+
                    # ddspL base and path same
                    ddspL(0:n,type="OTP",Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    ddspL(0:n,type="ITP",Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    ddspL(0:n,type="OTP",Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    ddspL(0:n,type="ITP",Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    ddspL(0:n,type="OSP",Ls.path=c(~`2`,~`2`))+
                    ddspL(0:n,type="ISP",Ls.path=c(~`2`,~`2`))+
                    # dnspL base and path same
                    dnspL(0:n,type="OTP",L.base=~`2`,Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    dnspL(0:n,type="ITP",L.base=~`2`,Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    dnspL(0:n,type="OTP",L.base=~`2`,Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    dnspL(0:n,type="ITP",L.base=~`2`,Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    dnspL(0:n,type="OSP",L.base=~`2`,Ls.path=c(~`2`,~`2`))+
                    dnspL(0:n,type="ISP",L.base=~`2`,Ls.path=c(~`2`,~`2`))+
                    # despL distinct base and one layer
                    despL(0:n,type="OTP",L.base=~`2`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    despL(0:n,type="ITP",L.base=~`2`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    despL(0:n,type="OTP",L.base=~`2`,Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    despL(0:n,type="ITP",L.base=~`2`,Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    despL(0:n,type="OSP",L.base=~`2`,Ls.path=c(~`2`,~`3`))+
                    despL(0:n,type="ISP",L.base=~`2`,Ls.path=c(~`2`,~`3`))+
                    # ddspL distinct base and one layer
                    ddspL(0:n,type="OTP",Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    ddspL(0:n,type="ITP",Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    ddspL(0:n,type="OTP",Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    ddspL(0:n,type="ITP",Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    ddspL(0:n,type="OSP",Ls.path=c(~`2`,~`3`))+
                    ddspL(0:n,type="ISP",Ls.path=c(~`2`,~`3`))+
                    # dnspL distinct base and one layer
                    dnspL(0:n,type="OTP",L.base=~`2`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    dnspL(0:n,type="ITP",L.base=~`2`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    dnspL(0:n,type="OTP",L.base=~`2`,Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    dnspL(0:n,type="ITP",L.base=~`2`,Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    dnspL(0:n,type="OSP",L.base=~`2`,Ls.path=c(~`2`,~`3`))+
                    dnspL(0:n,type="ISP",L.base=~`2`,Ls.path=c(~`2`,~`3`))+
                    # Geometrically weighted
                    # despL distinct layers
                    dgwespL(decay,fixed=TRUE,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    dgwespL(decay,fixed=TRUE,type="ITP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    dgwespL(decay,fixed=TRUE,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    dgwespL(decay,fixed=TRUE,type="ITP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    dgwespL(decay,fixed=TRUE,type="OSP",L.base=~`1`,Ls.path=c(~`2`,~`3`))+
                    dgwespL(decay,fixed=TRUE,type="ISP",L.base=~`1`,Ls.path=c(~`2`,~`3`))+
                    # ddspL distinct layers
                    dgwdspL(decay,fixed=TRUE,type="OTP",Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    dgwdspL(decay,fixed=TRUE,type="ITP",Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    dgwdspL(decay,fixed=TRUE,type="OTP",Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    dgwdspL(decay,fixed=TRUE,type="ITP",Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    dgwdspL(decay,fixed=TRUE,type="OSP",Ls.path=c(~`2`,~`3`))+
                    dgwdspL(decay,fixed=TRUE,type="ISP",Ls.path=c(~`2`,~`3`))+
                    # dnspL distinct layers
                    dgwnspL(decay,fixed=TRUE,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    dgwnspL(decay,fixed=TRUE,type="ITP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    dgwnspL(decay,fixed=TRUE,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    dgwnspL(decay,fixed=TRUE,type="ITP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    dgwnspL(decay,fixed=TRUE,type="OSP",L.base=~`1`,Ls.path=c(~`2`,~`3`))+
                    dgwnspL(decay,fixed=TRUE,type="ISP",L.base=~`1`,Ls.path=c(~`2`,~`3`))## +
                    ## # despL base and path distinct
                    ## dgwespL(decay,fixed=TRUE,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    ## dgwespL(decay,fixed=TRUE,type="ITP",L.base=~`1`,Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    ## dgwespL(decay,fixed=TRUE,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    ## dgwespL(decay,fixed=TRUE,type="ITP",L.base=~`1`,Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    ## dgwespL(decay,fixed=TRUE,type="OSP",L.base=~`1`,Ls.path=c(~`2`,~`2`))+
                    ## dgwespL(decay,fixed=TRUE,type="ISP",L.base=~`1`,Ls.path=c(~`2`,~`2`))+
                    ## # ddspL base and path distinct
                    ## dgwdspL(decay,fixed=TRUE,type="OTP",Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    ## dgwdspL(decay,fixed=TRUE,type="ITP",Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    ## dgwdspL(decay,fixed=TRUE,type="OTP",Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    ## dgwdspL(decay,fixed=TRUE,type="ITP",Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    ## dgwdspL(decay,fixed=TRUE,type="OSP",Ls.path=c(~`2`,~`2`))+
                    ## dgwdspL(decay,fixed=TRUE,type="ISP",Ls.path=c(~`2`,~`2`))+
                    ## # dnspL base and path distinct
                    ## dgwnspL(decay,fixed=TRUE,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    ## dgwnspL(decay,fixed=TRUE,type="ITP",L.base=~`1`,Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    ## dgwnspL(decay,fixed=TRUE,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    ## dgwnspL(decay,fixed=TRUE,type="ITP",L.base=~`1`,Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    ## dgwnspL(decay,fixed=TRUE,type="OSP",L.base=~`1`,Ls.path=c(~`2`,~`2`))+
                    ## dgwnspL(decay,fixed=TRUE,type="ISP",L.base=~`1`,Ls.path=c(~`2`,~`2`))+
                    ## # despL base and path same
                    ## dgwespL(decay,fixed=TRUE,type="OTP",L.base=~`2`,Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    ## dgwespL(decay,fixed=TRUE,type="ITP",L.base=~`2`,Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    ## dgwespL(decay,fixed=TRUE,type="OTP",L.base=~`2`,Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    ## dgwespL(decay,fixed=TRUE,type="ITP",L.base=~`2`,Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    ## dgwespL(decay,fixed=TRUE,type="OSP",L.base=~`2`,Ls.path=c(~`2`,~`2`))+
                    ## dgwespL(decay,fixed=TRUE,type="ISP",L.base=~`2`,Ls.path=c(~`2`,~`2`))+
                    ## # ddspL base and path same
                    ## dgwdspL(decay,fixed=TRUE,type="OTP",Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    ## dgwdspL(decay,fixed=TRUE,type="ITP",Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    ## dgwdspL(decay,fixed=TRUE,type="OTP",Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    ## dgwdspL(decay,fixed=TRUE,type="ITP",Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    ## dgwdspL(decay,fixed=TRUE,type="OSP",Ls.path=c(~`2`,~`2`))+
                    ## dgwdspL(decay,fixed=TRUE,type="ISP",Ls.path=c(~`2`,~`2`))+
                    ## # dnspL base and path same
                    ## dgwnspL(decay,fixed=TRUE,type="OTP",L.base=~`2`,Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    ## dgwnspL(decay,fixed=TRUE,type="ITP",L.base=~`2`,Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    ## dgwnspL(decay,fixed=TRUE,type="OTP",L.base=~`2`,Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    ## dgwnspL(decay,fixed=TRUE,type="ITP",L.base=~`2`,Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    ## dgwnspL(decay,fixed=TRUE,type="OSP",L.base=~`2`,Ls.path=c(~`2`,~`2`))+
                    ## dgwnspL(decay,fixed=TRUE,type="ISP",L.base=~`2`,Ls.path=c(~`2`,~`2`))+
                    ## # despL distinct base and one layer
                    ## dgwespL(decay,fixed=TRUE,type="OTP",L.base=~`2`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    ## dgwespL(decay,fixed=TRUE,type="ITP",L.base=~`2`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    ## dgwespL(decay,fixed=TRUE,type="OTP",L.base=~`2`,Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    ## dgwespL(decay,fixed=TRUE,type="ITP",L.base=~`2`,Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    ## dgwespL(decay,fixed=TRUE,type="OSP",L.base=~`2`,Ls.path=c(~`2`,~`3`))+
                    ## dgwespL(decay,fixed=TRUE,type="ISP",L.base=~`2`,Ls.path=c(~`2`,~`3`))+
                    ## # ddspL distinct base and one layer
                    ## dgwdspL(decay,fixed=TRUE,type="OTP",Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    ## dgwdspL(decay,fixed=TRUE,type="ITP",Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    ## dgwdspL(decay,fixed=TRUE,type="OTP",Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    ## dgwdspL(decay,fixed=TRUE,type="ITP",Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    ## dgwdspL(decay,fixed=TRUE,type="OSP",Ls.path=c(~`2`,~`3`))+
                    ## dgwdspL(decay,fixed=TRUE,type="ISP",Ls.path=c(~`2`,~`3`))+
                    ## # dnspL distinct base and one layer
                    ## dgwnspL(decay,fixed=TRUE,type="OTP",L.base=~`2`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    ## dgwnspL(decay,fixed=TRUE,type="ITP",L.base=~`2`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    ## dgwnspL(decay,fixed=TRUE,type="OTP",L.base=~`2`,Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    ## dgwnspL(decay,fixed=TRUE,type="ITP",L.base=~`2`,Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    ## dgwnspL(decay,fixed=TRUE,type="OSP",L.base=~`2`,Ls.path=c(~`2`,~`3`))+
                    ## dgwnspL(decay,fixed=TRUE,type="ISP",L.base=~`2`,Ls.path=c(~`2`,~`3`))
                 ,
                  control=ctrl,
                  nsim=200))

  stats <- sapply(sim,
                  function(nw){
                    n <- network.size(nw)/3
                    m <- as.matrix(nw)
                    m1 <- m[seq_len(n),seq_len(n)]
                    m2 <- m[seq_len(n)+n,seq_len(n)+n]
                    m3 <- m[seq_len(n)+n*2,seq_len(n)+n*2]
                    
                    c(
                      # despL distinct layers
                      despL(0:n,OTP,L.base=m1,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      despL(0:n,ITP,L.base=m1,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      despL(0:n,OTP,L.base=m1,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      despL(0:n,ITP,L.base=m1,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      despL(0:n,OSP,L.base=m1,Ls.path1=m2,Ls.path2=m3),
                      despL(0:n,ISP,L.base=m1,Ls.path1=m2,Ls.path2=m3),
                      # ddspL distinct layers
                      ddspL(0:n,OTP,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      ddspL(0:n,ITP,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      ddspL(0:n,OTP,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      ddspL(0:n,ITP,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      ddspL(0:n,OSP,Ls.path1=m2,Ls.path2=m3),
                      ddspL(0:n,ISP,Ls.path1=m2,Ls.path2=m3),
                      # dnspL distinct layers
                      dnspL(0:n,OTP,L.base=m1,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      dnspL(0:n,ITP,L.base=m1,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      dnspL(0:n,OTP,L.base=m1,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      dnspL(0:n,ITP,L.base=m1,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      dnspL(0:n,OSP,L.base=m1,Ls.path1=m2,Ls.path2=m3),
                      dnspL(0:n,ISP,L.base=m1,Ls.path1=m2,Ls.path2=m3),
                      # despL base and path distinct
                      despL(0:n,OTP,L.base=m1,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      despL(0:n,ITP,L.base=m1,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      despL(0:n,OTP,L.base=m1,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      despL(0:n,ITP,L.base=m1,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      despL(0:n,OSP,L.base=m1,Ls.path1=m2,Ls.path2=m2),
                      despL(0:n,ISP,L.base=m1,Ls.path1=m2,Ls.path2=m2),
                      # ddspL base and path distinct
                      ddspL(0:n,OTP,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      ddspL(0:n,ITP,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      ddspL(0:n,OTP,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      ddspL(0:n,ITP,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      ddspL(0:n,OSP,Ls.path1=m2,Ls.path2=m2),
                      ddspL(0:n,ISP,Ls.path1=m2,Ls.path2=m2),
                      # dnspL base and path distinct
                      dnspL(0:n,OTP,L.base=m1,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      dnspL(0:n,ITP,L.base=m1,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      dnspL(0:n,OTP,L.base=m1,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      dnspL(0:n,ITP,L.base=m1,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      dnspL(0:n,OSP,L.base=m1,Ls.path1=m2,Ls.path2=m2),
                      dnspL(0:n,ISP,L.base=m1,Ls.path1=m2,Ls.path2=m2),
                      # despL base and path same
                      despL(0:n,OTP,L.base=m2,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      despL(0:n,ITP,L.base=m2,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      despL(0:n,OTP,L.base=m2,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      despL(0:n,ITP,L.base=m2,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      despL(0:n,OSP,L.base=m2,Ls.path1=m2,Ls.path2=m2),
                      despL(0:n,ISP,L.base=m2,Ls.path1=m2,Ls.path2=m2),
                      # ddspL base and path same
                      ddspL(0:n,OTP,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      ddspL(0:n,ITP,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      ddspL(0:n,OTP,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      ddspL(0:n,ITP,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      ddspL(0:n,OSP,Ls.path1=m2,Ls.path2=m2),
                      ddspL(0:n,ISP,Ls.path1=m2,Ls.path2=m2),
                      # dnspL base and path same
                      dnspL(0:n,OTP,L.base=m2,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      dnspL(0:n,ITP,L.base=m2,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      dnspL(0:n,OTP,L.base=m2,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      dnspL(0:n,ITP,L.base=m2,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      dnspL(0:n,OSP,L.base=m2,Ls.path1=m2,Ls.path2=m2),
                      dnspL(0:n,ISP,L.base=m2,Ls.path1=m2,Ls.path2=m2),
                      # despL distinct base and one layer
                      despL(0:n,OTP,L.base=m2,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      despL(0:n,ITP,L.base=m2,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      despL(0:n,OTP,L.base=m2,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      despL(0:n,ITP,L.base=m2,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      despL(0:n,OSP,L.base=m2,Ls.path1=m2,Ls.path2=m3),
                      despL(0:n,ISP,L.base=m2,Ls.path1=m2,Ls.path2=m3),
                      # ddspL distinct base and one layer
                      ddspL(0:n,OTP,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      ddspL(0:n,ITP,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      ddspL(0:n,OTP,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      ddspL(0:n,ITP,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      ddspL(0:n,OSP,Ls.path1=m2,Ls.path2=m3),
                      ddspL(0:n,ISP,Ls.path1=m2,Ls.path2=m3),
                      # dnspL distinct base and one layer
                      dnspL(0:n,OTP,L.base=m2,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      dnspL(0:n,ITP,L.base=m2,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      dnspL(0:n,OTP,L.base=m2,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      dnspL(0:n,ITP,L.base=m2,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      dnspL(0:n,OSP,L.base=m2,Ls.path1=m2,Ls.path2=m3),
                      dnspL(0:n,ISP,L.base=m2,Ls.path1=m2,Ls.path2=m3),
                      # Geometrically weighted
                      # despL distinct layers
                      dgwespL(decay,n,OTP,L.base=m1,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      dgwespL(decay,n,ITP,L.base=m1,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      dgwespL(decay,n,OTP,L.base=m1,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      dgwespL(decay,n,ITP,L.base=m1,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      dgwespL(decay,n,OSP,L.base=m1,Ls.path1=m2,Ls.path2=m3),
                      dgwespL(decay,n,ISP,L.base=m1,Ls.path1=m2,Ls.path2=m3),
                      # ddspL distinct layers
                      dgwdspL(decay,n,OTP,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      dgwdspL(decay,n,ITP,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      dgwdspL(decay,n,OTP,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      dgwdspL(decay,n,ITP,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      dgwdspL(decay,n,OSP,Ls.path1=m2,Ls.path2=m3),
                      dgwdspL(decay,n,ISP,Ls.path1=m2,Ls.path2=m3),
                      # dnspL distinct layers
                      dgwnspL(decay,n,OTP,L.base=m1,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      dgwnspL(decay,n,ITP,L.base=m1,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      dgwnspL(decay,n,OTP,L.base=m1,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      dgwnspL(decay,n,ITP,L.base=m1,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      dgwnspL(decay,n,OSP,L.base=m1,Ls.path1=m2,Ls.path2=m3),
                      dgwnspL(decay,n,ISP,L.base=m1,Ls.path1=m2,Ls.path2=m3)## ,
                      ## # despL base and path distinct
                      ## dgwespL(decay,n,OTP,L.base=m1,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      ## dgwespL(decay,n,ITP,L.base=m1,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      ## dgwespL(decay,n,OTP,L.base=m1,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      ## dgwespL(decay,n,ITP,L.base=m1,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      ## dgwespL(decay,n,OSP,L.base=m1,Ls.path1=m2,Ls.path2=m2),
                      ## dgwespL(decay,n,ISP,L.base=m1,Ls.path1=m2,Ls.path2=m2),
                      ## # ddspL base and path distinct
                      ## dgwdspL(decay,n,OTP,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      ## dgwdspL(decay,n,ITP,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      ## dgwdspL(decay,n,OTP,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      ## dgwdspL(decay,n,ITP,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      ## dgwdspL(decay,n,OSP,Ls.path1=m2,Ls.path2=m2),
                      ## dgwdspL(decay,n,ISP,Ls.path1=m2,Ls.path2=m2),
                      ## # dnspL base and path distinct
                      ## dgwnspL(decay,n,OTP,L.base=m1,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      ## dgwnspL(decay,n,ITP,L.base=m1,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      ## dgwnspL(decay,n,OTP,L.base=m1,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      ## dgwnspL(decay,n,ITP,L.base=m1,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      ## dgwnspL(decay,n,OSP,L.base=m1,Ls.path1=m2,Ls.path2=m2),
                      ## dgwnspL(decay,n,ISP,L.base=m1,Ls.path1=m2,Ls.path2=m2),
                      ## # despL base and path same
                      ## dgwespL(decay,n,OTP,L.base=m2,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      ## dgwespL(decay,n,ITP,L.base=m2,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      ## dgwespL(decay,n,OTP,L.base=m2,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      ## dgwespL(decay,n,ITP,L.base=m2,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      ## dgwespL(decay,n,OSP,L.base=m2,Ls.path1=m2,Ls.path2=m2),
                      ## dgwespL(decay,n,ISP,L.base=m2,Ls.path1=m2,Ls.path2=m2),
                      ## # ddspL base and path same
                      ## dgwdspL(decay,n,OTP,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      ## dgwdspL(decay,n,ITP,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      ## dgwdspL(decay,n,OTP,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      ## dgwdspL(decay,n,ITP,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      ## dgwdspL(decay,n,OSP,Ls.path1=m2,Ls.path2=m2),
                      ## dgwdspL(decay,n,ISP,Ls.path1=m2,Ls.path2=m2),
                      ## # dnspL base and path same
                      ## dgwnspL(decay,n,OTP,L.base=m2,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      ## dgwnspL(decay,n,ITP,L.base=m2,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      ## dgwnspL(decay,n,OTP,L.base=m2,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      ## dgwnspL(decay,n,ITP,L.base=m2,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      ## dgwnspL(decay,n,OSP,L.base=m2,Ls.path1=m2,Ls.path2=m2),
                      ## dgwnspL(decay,n,ISP,L.base=m2,Ls.path1=m2,Ls.path2=m2),
                      ## # despL distinct base and one layer
                      ## dgwespL(decay,n,OTP,L.base=m2,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      ## dgwespL(decay,n,ITP,L.base=m2,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      ## dgwespL(decay,n,OTP,L.base=m2,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      ## dgwespL(decay,n,ITP,L.base=m2,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      ## dgwespL(decay,n,OSP,L.base=m2,Ls.path1=m2,Ls.path2=m3),
                      ## dgwespL(decay,n,ISP,L.base=m2,Ls.path1=m2,Ls.path2=m3),
                      ## # ddspL distinct base and one layer
                      ## dgwdspL(decay,n,OTP,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      ## dgwdspL(decay,n,ITP,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      ## dgwdspL(decay,n,OTP,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      ## dgwdspL(decay,n,ITP,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      ## dgwdspL(decay,n,OSP,Ls.path1=m2,Ls.path2=m3),
                      ## dgwdspL(decay,n,ISP,Ls.path1=m2,Ls.path2=m3),
                      ## # dnspL distinct base and one layer
                      ## dgwnspL(decay,n,OTP,L.base=m2,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      ## dgwnspL(decay,n,ITP,L.base=m2,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      ## dgwnspL(decay,n,OTP,L.base=m2,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      ## dgwnspL(decay,n,ITP,L.base=m2,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      ## dgwnspL(decay,n,OSP,L.base=m2,Ls.path1=m2,Ls.path2=m3),
                      ## dgwnspL(decay,n,ISP,L.base=m2,Ls.path1=m2,Ls.path2=m3)
                    )
                  }) %>% t()

  expect_equivalent(attr(sim,"stats"), stats)
})

### Heterogeneous directedness.
nw1 <- nw3 <- network.initialize(n,dir=TRUE)
nw2 <- network.initialize(n,dir=FALSE)
lnw <- Layer(nw1,nw2,nw3)

test_that(paste("Multilayer dgw*sp statistics for heterogeneously directed networks 1",sptxt), {
  sim <- suppressWarnings(simulate(lnw~
                    # despL distinct layers
                    despL(0:n,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    despL(0:n,type="ITP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    despL(0:n,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    despL(0:n,type="ITP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    despL(0:n,type="OSP",L.base=~`1`,Ls.path=c(~`2`,~`3`))+
                    despL(0:n,type="ISP",L.base=~`1`,Ls.path=c(~`2`,~`3`))+
                    # ddspL distinct layers
                    ddspL(0:n,type="OTP",Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    ddspL(0:n,type="ITP",Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    ddspL(0:n,type="OTP",Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    ddspL(0:n,type="ITP",Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    ddspL(0:n,type="OSP",Ls.path=c(~`2`,~`3`))+
                    ddspL(0:n,type="ISP",Ls.path=c(~`2`,~`3`))+
                    # dnspL distinct layers
                    dnspL(0:n,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    dnspL(0:n,type="ITP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    dnspL(0:n,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    dnspL(0:n,type="ITP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    dnspL(0:n,type="OSP",L.base=~`1`,Ls.path=c(~`2`,~`3`))+
                    dnspL(0:n,type="ISP",L.base=~`1`,Ls.path=c(~`2`,~`3`))+
                    # despL base and path distinct
                    despL(0:n,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    despL(0:n,type="ITP",L.base=~`1`,Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    despL(0:n,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    despL(0:n,type="ITP",L.base=~`1`,Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    despL(0:n,type="OSP",L.base=~`1`,Ls.path=c(~`2`,~`2`))+
                    despL(0:n,type="ISP",L.base=~`1`,Ls.path=c(~`2`,~`2`))+
                    # ddspL base and path distinct
                    ddspL(0:n,type="OTP",Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    ddspL(0:n,type="ITP",Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    ddspL(0:n,type="OTP",Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    ddspL(0:n,type="ITP",Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    ddspL(0:n,type="OSP",Ls.path=c(~`2`,~`2`))+
                    ddspL(0:n,type="ISP",Ls.path=c(~`2`,~`2`))+
                    # dnspL base and path distinct
                    dnspL(0:n,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    dnspL(0:n,type="ITP",L.base=~`1`,Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    dnspL(0:n,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    dnspL(0:n,type="ITP",L.base=~`1`,Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    dnspL(0:n,type="OSP",L.base=~`1`,Ls.path=c(~`2`,~`2`))+
                    dnspL(0:n,type="ISP",L.base=~`1`,Ls.path=c(~`2`,~`2`))+
                    # despL base and path same
                    despL(0:n,type="OTP",L.base=~`2`,Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    despL(0:n,type="ITP",L.base=~`2`,Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    despL(0:n,type="OTP",L.base=~`2`,Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    despL(0:n,type="ITP",L.base=~`2`,Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    despL(0:n,type="OSP",L.base=~`2`,Ls.path=c(~`2`,~`2`))+
                    despL(0:n,type="ISP",L.base=~`2`,Ls.path=c(~`2`,~`2`))+
                    # ddspL base and path same
                    ddspL(0:n,type="OTP",Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    ddspL(0:n,type="ITP",Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    ddspL(0:n,type="OTP",Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    ddspL(0:n,type="ITP",Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    ddspL(0:n,type="OSP",Ls.path=c(~`2`,~`2`))+
                    ddspL(0:n,type="ISP",Ls.path=c(~`2`,~`2`))+
                    # dnspL base and path same
                    dnspL(0:n,type="OTP",L.base=~`2`,Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    dnspL(0:n,type="ITP",L.base=~`2`,Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    dnspL(0:n,type="OTP",L.base=~`2`,Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    dnspL(0:n,type="ITP",L.base=~`2`,Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    dnspL(0:n,type="OSP",L.base=~`2`,Ls.path=c(~`2`,~`2`))+
                    dnspL(0:n,type="ISP",L.base=~`2`,Ls.path=c(~`2`,~`2`))+
                    # despL distinct base and one layer
                    despL(0:n,type="OTP",L.base=~`2`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    despL(0:n,type="ITP",L.base=~`2`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    despL(0:n,type="OTP",L.base=~`2`,Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    despL(0:n,type="ITP",L.base=~`2`,Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    despL(0:n,type="OSP",L.base=~`2`,Ls.path=c(~`2`,~`3`))+
                    despL(0:n,type="ISP",L.base=~`2`,Ls.path=c(~`2`,~`3`))+
                    # ddspL distinct base and one layer
                    ddspL(0:n,type="OTP",Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    ddspL(0:n,type="ITP",Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    ddspL(0:n,type="OTP",Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    ddspL(0:n,type="ITP",Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    ddspL(0:n,type="OSP",Ls.path=c(~`2`,~`3`))+
                    ddspL(0:n,type="ISP",Ls.path=c(~`2`,~`3`))+
                    # dnspL distinct base and one layer
                    dnspL(0:n,type="OTP",L.base=~`2`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    dnspL(0:n,type="ITP",L.base=~`2`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    dnspL(0:n,type="OTP",L.base=~`2`,Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    dnspL(0:n,type="ITP",L.base=~`2`,Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    dnspL(0:n,type="OSP",L.base=~`2`,Ls.path=c(~`2`,~`3`))+
                    dnspL(0:n,type="ISP",L.base=~`2`,Ls.path=c(~`2`,~`3`))+
                    # Geometrically weighted
                    # despL distinct layers
                    dgwespL(decay,fixed=TRUE,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    dgwespL(decay,fixed=TRUE,type="ITP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    dgwespL(decay,fixed=TRUE,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    dgwespL(decay,fixed=TRUE,type="ITP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    dgwespL(decay,fixed=TRUE,type="OSP",L.base=~`1`,Ls.path=c(~`2`,~`3`))+
                    dgwespL(decay,fixed=TRUE,type="ISP",L.base=~`1`,Ls.path=c(~`2`,~`3`))+
                    # ddspL distinct layers
                    dgwdspL(decay,fixed=TRUE,type="OTP",Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    dgwdspL(decay,fixed=TRUE,type="ITP",Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    dgwdspL(decay,fixed=TRUE,type="OTP",Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    dgwdspL(decay,fixed=TRUE,type="ITP",Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    dgwdspL(decay,fixed=TRUE,type="OSP",Ls.path=c(~`2`,~`3`))+
                    dgwdspL(decay,fixed=TRUE,type="ISP",Ls.path=c(~`2`,~`3`))+
                    # dnspL distinct layers
                    dgwnspL(decay,fixed=TRUE,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    dgwnspL(decay,fixed=TRUE,type="ITP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    dgwnspL(decay,fixed=TRUE,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    dgwnspL(decay,fixed=TRUE,type="ITP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    dgwnspL(decay,fixed=TRUE,type="OSP",L.base=~`1`,Ls.path=c(~`2`,~`3`))+
                    dgwnspL(decay,fixed=TRUE,type="ISP",L.base=~`1`,Ls.path=c(~`2`,~`3`))## +
                    ## # despL base and path distinct
                    ## dgwespL(decay,fixed=TRUE,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    ## dgwespL(decay,fixed=TRUE,type="ITP",L.base=~`1`,Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    ## dgwespL(decay,fixed=TRUE,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    ## dgwespL(decay,fixed=TRUE,type="ITP",L.base=~`1`,Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    ## dgwespL(decay,fixed=TRUE,type="OSP",L.base=~`1`,Ls.path=c(~`2`,~`2`))+
                    ## dgwespL(decay,fixed=TRUE,type="ISP",L.base=~`1`,Ls.path=c(~`2`,~`2`))+
                    ## # ddspL base and path distinct
                    ## dgwdspL(decay,fixed=TRUE,type="OTP",Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    ## dgwdspL(decay,fixed=TRUE,type="ITP",Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    ## dgwdspL(decay,fixed=TRUE,type="OTP",Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    ## dgwdspL(decay,fixed=TRUE,type="ITP",Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    ## dgwdspL(decay,fixed=TRUE,type="OSP",Ls.path=c(~`2`,~`2`))+
                    ## dgwdspL(decay,fixed=TRUE,type="ISP",Ls.path=c(~`2`,~`2`))+
                    ## # dnspL base and path distinct
                    ## dgwnspL(decay,fixed=TRUE,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    ## dgwnspL(decay,fixed=TRUE,type="ITP",L.base=~`1`,Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    ## dgwnspL(decay,fixed=TRUE,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    ## dgwnspL(decay,fixed=TRUE,type="ITP",L.base=~`1`,Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    ## dgwnspL(decay,fixed=TRUE,type="OSP",L.base=~`1`,Ls.path=c(~`2`,~`2`))+
                    ## dgwnspL(decay,fixed=TRUE,type="ISP",L.base=~`1`,Ls.path=c(~`2`,~`2`))+
                    ## # despL base and path same
                    ## dgwespL(decay,fixed=TRUE,type="OTP",L.base=~`2`,Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    ## dgwespL(decay,fixed=TRUE,type="ITP",L.base=~`2`,Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    ## dgwespL(decay,fixed=TRUE,type="OTP",L.base=~`2`,Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    ## dgwespL(decay,fixed=TRUE,type="ITP",L.base=~`2`,Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    ## dgwespL(decay,fixed=TRUE,type="OSP",L.base=~`2`,Ls.path=c(~`2`,~`2`))+
                    ## dgwespL(decay,fixed=TRUE,type="ISP",L.base=~`2`,Ls.path=c(~`2`,~`2`))+
                    ## # ddspL base and path same
                    ## dgwdspL(decay,fixed=TRUE,type="OTP",Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    ## dgwdspL(decay,fixed=TRUE,type="ITP",Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    ## dgwdspL(decay,fixed=TRUE,type="OTP",Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    ## dgwdspL(decay,fixed=TRUE,type="ITP",Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    ## dgwdspL(decay,fixed=TRUE,type="OSP",Ls.path=c(~`2`,~`2`))+
                    ## dgwdspL(decay,fixed=TRUE,type="ISP",Ls.path=c(~`2`,~`2`))+
                    ## # dnspL base and path same
                    ## dgwnspL(decay,fixed=TRUE,type="OTP",L.base=~`2`,Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    ## dgwnspL(decay,fixed=TRUE,type="ITP",L.base=~`2`,Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    ## dgwnspL(decay,fixed=TRUE,type="OTP",L.base=~`2`,Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    ## dgwnspL(decay,fixed=TRUE,type="ITP",L.base=~`2`,Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    ## dgwnspL(decay,fixed=TRUE,type="OSP",L.base=~`2`,Ls.path=c(~`2`,~`2`))+
                    ## dgwnspL(decay,fixed=TRUE,type="ISP",L.base=~`2`,Ls.path=c(~`2`,~`2`))+
                    ## # despL distinct base and one layer
                    ## dgwespL(decay,fixed=TRUE,type="OTP",L.base=~`2`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    ## dgwespL(decay,fixed=TRUE,type="ITP",L.base=~`2`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    ## dgwespL(decay,fixed=TRUE,type="OTP",L.base=~`2`,Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    ## dgwespL(decay,fixed=TRUE,type="ITP",L.base=~`2`,Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    ## dgwespL(decay,fixed=TRUE,type="OSP",L.base=~`2`,Ls.path=c(~`2`,~`3`))+
                    ## dgwespL(decay,fixed=TRUE,type="ISP",L.base=~`2`,Ls.path=c(~`2`,~`3`))+
                    ## # ddspL distinct base and one layer
                    ## dgwdspL(decay,fixed=TRUE,type="OTP",Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    ## dgwdspL(decay,fixed=TRUE,type="ITP",Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    ## dgwdspL(decay,fixed=TRUE,type="OTP",Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    ## dgwdspL(decay,fixed=TRUE,type="ITP",Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    ## dgwdspL(decay,fixed=TRUE,type="OSP",Ls.path=c(~`2`,~`3`))+
                    ## dgwdspL(decay,fixed=TRUE,type="ISP",Ls.path=c(~`2`,~`3`))+
                    ## # dnspL distinct base and one layer
                    ## dgwnspL(decay,fixed=TRUE,type="OTP",L.base=~`2`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    ## dgwnspL(decay,fixed=TRUE,type="ITP",L.base=~`2`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    ## dgwnspL(decay,fixed=TRUE,type="OTP",L.base=~`2`,Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    ## dgwnspL(decay,fixed=TRUE,type="ITP",L.base=~`2`,Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    ## dgwnspL(decay,fixed=TRUE,type="OSP",L.base=~`2`,Ls.path=c(~`2`,~`3`))+
                    ## dgwnspL(decay,fixed=TRUE,type="ISP",L.base=~`2`,Ls.path=c(~`2`,~`3`))
                 ,
                  control=ctrl,
                  nsim=200))

  stats <- sapply(sim,
                  function(nw){
                    n <- network.size(nw)/3
                    m <- as.matrix(nw)
                    m1 <- m[seq_len(n),seq_len(n)]
                    m2 <- m[seq_len(n)+n,seq_len(n)+n]
                    m2 <- m2 + t(m2)
                    m3 <- m[seq_len(n)+n*2,seq_len(n)+n*2]
                    
                    c(
                      # despL distinct layers
                      despL(0:n,OTP,L.base=m1,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      despL(0:n,ITP,L.base=m1,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      despL(0:n,OTP,L.base=m1,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      despL(0:n,ITP,L.base=m1,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      despL(0:n,OSP,L.base=m1,Ls.path1=m2,Ls.path2=m3),
                      despL(0:n,ISP,L.base=m1,Ls.path1=m2,Ls.path2=m3),
                      # ddspL distinct layers
                      ddspL(0:n,OTP,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      ddspL(0:n,ITP,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      ddspL(0:n,OTP,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      ddspL(0:n,ITP,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      ddspL(0:n,OSP,Ls.path1=m2,Ls.path2=m3),
                      ddspL(0:n,ISP,Ls.path1=m2,Ls.path2=m3),
                      # dnspL distinct layers
                      dnspL(0:n,OTP,L.base=m1,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      dnspL(0:n,ITP,L.base=m1,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      dnspL(0:n,OTP,L.base=m1,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      dnspL(0:n,ITP,L.base=m1,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      dnspL(0:n,OSP,L.base=m1,Ls.path1=m2,Ls.path2=m3),
                      dnspL(0:n,ISP,L.base=m1,Ls.path1=m2,Ls.path2=m3),
                      # despL base and path distinct
                      despL(0:n,OTP,L.base=m1,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      despL(0:n,ITP,L.base=m1,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      despL(0:n,OTP,L.base=m1,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      despL(0:n,ITP,L.base=m1,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      despL(0:n,OSP,L.base=m1,Ls.path1=m2,Ls.path2=m2),
                      despL(0:n,ISP,L.base=m1,Ls.path1=m2,Ls.path2=m2),
                      # ddspL base and path distinct
                      ddspL(0:n,OTP,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      ddspL(0:n,ITP,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      ddspL(0:n,OTP,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      ddspL(0:n,ITP,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      ddspL(0:n,OSP,Ls.path1=m2,Ls.path2=m2),
                      ddspL(0:n,ISP,Ls.path1=m2,Ls.path2=m2),
                      # dnspL base and path distinct
                      dnspL(0:n,OTP,L.base=m1,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      dnspL(0:n,ITP,L.base=m1,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      dnspL(0:n,OTP,L.base=m1,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      dnspL(0:n,ITP,L.base=m1,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      dnspL(0:n,OSP,L.base=m1,Ls.path1=m2,Ls.path2=m2),
                      dnspL(0:n,ISP,L.base=m1,Ls.path1=m2,Ls.path2=m2),
                      # despL base and path same
                      despL(0:n,OTP,L.base=m2,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      despL(0:n,ITP,L.base=m2,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      despL(0:n,OTP,L.base=m2,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      despL(0:n,ITP,L.base=m2,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      despL(0:n,OSP,L.base=m2,Ls.path1=m2,Ls.path2=m2),
                      despL(0:n,ISP,L.base=m2,Ls.path1=m2,Ls.path2=m2),
                      # ddspL base and path same
                      ddspL(0:n,OTP,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      ddspL(0:n,ITP,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      ddspL(0:n,OTP,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      ddspL(0:n,ITP,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      ddspL(0:n,OSP,Ls.path1=m2,Ls.path2=m2),
                      ddspL(0:n,ISP,Ls.path1=m2,Ls.path2=m2),
                      # dnspL base and path same
                      dnspL(0:n,OTP,L.base=m2,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      dnspL(0:n,ITP,L.base=m2,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      dnspL(0:n,OTP,L.base=m2,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      dnspL(0:n,ITP,L.base=m2,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      dnspL(0:n,OSP,L.base=m2,Ls.path1=m2,Ls.path2=m2),
                      dnspL(0:n,ISP,L.base=m2,Ls.path1=m2,Ls.path2=m2),
                      # despL distinct base and one layer
                      despL(0:n,OTP,L.base=m2,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      despL(0:n,ITP,L.base=m2,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      despL(0:n,OTP,L.base=m2,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      despL(0:n,ITP,L.base=m2,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      despL(0:n,OSP,L.base=m2,Ls.path1=m2,Ls.path2=m3),
                      despL(0:n,ISP,L.base=m2,Ls.path1=m2,Ls.path2=m3),
                      # ddspL distinct base and one layer
                      ddspL(0:n,OTP,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      ddspL(0:n,ITP,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      ddspL(0:n,OTP,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      ddspL(0:n,ITP,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      ddspL(0:n,OSP,Ls.path1=m2,Ls.path2=m3),
                      ddspL(0:n,ISP,Ls.path1=m2,Ls.path2=m3),
                      # dnspL distinct base and one layer
                      dnspL(0:n,OTP,L.base=m2,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      dnspL(0:n,ITP,L.base=m2,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      dnspL(0:n,OTP,L.base=m2,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      dnspL(0:n,ITP,L.base=m2,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      dnspL(0:n,OSP,L.base=m2,Ls.path1=m2,Ls.path2=m3),
                      dnspL(0:n,ISP,L.base=m2,Ls.path1=m2,Ls.path2=m3),
                      # Geometrically weighted
                      # despL distinct layers
                      dgwespL(decay,n,OTP,L.base=m1,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      dgwespL(decay,n,ITP,L.base=m1,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      dgwespL(decay,n,OTP,L.base=m1,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      dgwespL(decay,n,ITP,L.base=m1,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      dgwespL(decay,n,OSP,L.base=m1,Ls.path1=m2,Ls.path2=m3),
                      dgwespL(decay,n,ISP,L.base=m1,Ls.path1=m2,Ls.path2=m3),
                      # ddspL distinct layers
                      dgwdspL(decay,n,OTP,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      dgwdspL(decay,n,ITP,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      dgwdspL(decay,n,OTP,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      dgwdspL(decay,n,ITP,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      dgwdspL(decay,n,OSP,Ls.path1=m2,Ls.path2=m3),
                      dgwdspL(decay,n,ISP,Ls.path1=m2,Ls.path2=m3),
                      # dnspL distinct layers
                      dgwnspL(decay,n,OTP,L.base=m1,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      dgwnspL(decay,n,ITP,L.base=m1,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      dgwnspL(decay,n,OTP,L.base=m1,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      dgwnspL(decay,n,ITP,L.base=m1,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      dgwnspL(decay,n,OSP,L.base=m1,Ls.path1=m2,Ls.path2=m3),
                      dgwnspL(decay,n,ISP,L.base=m1,Ls.path1=m2,Ls.path2=m3)## ,
                      ## # despL base and path distinct
                      ## dgwespL(decay,n,OTP,L.base=m1,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      ## dgwespL(decay,n,ITP,L.base=m1,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      ## dgwespL(decay,n,OTP,L.base=m1,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      ## dgwespL(decay,n,ITP,L.base=m1,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      ## dgwespL(decay,n,OSP,L.base=m1,Ls.path1=m2,Ls.path2=m2),
                      ## dgwespL(decay,n,ISP,L.base=m1,Ls.path1=m2,Ls.path2=m2),
                      ## # ddspL base and path distinct
                      ## dgwdspL(decay,n,OTP,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      ## dgwdspL(decay,n,ITP,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      ## dgwdspL(decay,n,OTP,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      ## dgwdspL(decay,n,ITP,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      ## dgwdspL(decay,n,OSP,Ls.path1=m2,Ls.path2=m2),
                      ## dgwdspL(decay,n,ISP,Ls.path1=m2,Ls.path2=m2),
                      ## # dnspL base and path distinct
                      ## dgwnspL(decay,n,OTP,L.base=m1,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      ## dgwnspL(decay,n,ITP,L.base=m1,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      ## dgwnspL(decay,n,OTP,L.base=m1,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      ## dgwnspL(decay,n,ITP,L.base=m1,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      ## dgwnspL(decay,n,OSP,L.base=m1,Ls.path1=m2,Ls.path2=m2),
                      ## dgwnspL(decay,n,ISP,L.base=m1,Ls.path1=m2,Ls.path2=m2),
                      ## # despL base and path same
                      ## dgwespL(decay,n,OTP,L.base=m2,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      ## dgwespL(decay,n,ITP,L.base=m2,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      ## dgwespL(decay,n,OTP,L.base=m2,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      ## dgwespL(decay,n,ITP,L.base=m2,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      ## dgwespL(decay,n,OSP,L.base=m2,Ls.path1=m2,Ls.path2=m2),
                      ## dgwespL(decay,n,ISP,L.base=m2,Ls.path1=m2,Ls.path2=m2),
                      ## # ddspL base and path same
                      ## dgwdspL(decay,n,OTP,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      ## dgwdspL(decay,n,ITP,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      ## dgwdspL(decay,n,OTP,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      ## dgwdspL(decay,n,ITP,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      ## dgwdspL(decay,n,OSP,Ls.path1=m2,Ls.path2=m2),
                      ## dgwdspL(decay,n,ISP,Ls.path1=m2,Ls.path2=m2),
                      ## # dnspL base and path same
                      ## dgwnspL(decay,n,OTP,L.base=m2,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      ## dgwnspL(decay,n,ITP,L.base=m2,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      ## dgwnspL(decay,n,OTP,L.base=m2,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      ## dgwnspL(decay,n,ITP,L.base=m2,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      ## dgwnspL(decay,n,OSP,L.base=m2,Ls.path1=m2,Ls.path2=m2),
                      ## dgwnspL(decay,n,ISP,L.base=m2,Ls.path1=m2,Ls.path2=m2),
                      ## # despL distinct base and one layer
                      ## dgwespL(decay,n,OTP,L.base=m2,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      ## dgwespL(decay,n,ITP,L.base=m2,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      ## dgwespL(decay,n,OTP,L.base=m2,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      ## dgwespL(decay,n,ITP,L.base=m2,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      ## dgwespL(decay,n,OSP,L.base=m2,Ls.path1=m2,Ls.path2=m3),
                      ## dgwespL(decay,n,ISP,L.base=m2,Ls.path1=m2,Ls.path2=m3),
                      ## # ddspL distinct base and one layer
                      ## dgwdspL(decay,n,OTP,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      ## dgwdspL(decay,n,ITP,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      ## dgwdspL(decay,n,OTP,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      ## dgwdspL(decay,n,ITP,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      ## dgwdspL(decay,n,OSP,Ls.path1=m2,Ls.path2=m3),
                      ## dgwdspL(decay,n,ISP,Ls.path1=m2,Ls.path2=m3),
                      ## # dnspL distinct base and one layer
                      ## dgwnspL(decay,n,OTP,L.base=m2,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      ## dgwnspL(decay,n,ITP,L.base=m2,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      ## dgwnspL(decay,n,OTP,L.base=m2,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      ## dgwnspL(decay,n,ITP,L.base=m2,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      ## dgwnspL(decay,n,OSP,L.base=m2,Ls.path1=m2,Ls.path2=m3),
                      ## dgwnspL(decay,n,ISP,L.base=m2,Ls.path1=m2,Ls.path2=m3)
                    )
                  }) %>% t()

  expect_equivalent(attr(sim,"stats"), stats)
})


### Undirected.
nw1 <- nw2 <- nw3 <- network.initialize(n,dir=FALSE)
lnw <- Layer(nw1,nw2,nw3)

test_that(paste("Multilayer dgw*sp statistics for undirected networks",sptxt), {
  sim <- suppressWarnings(simulate(lnw~
                    # despL distinct layers
                    despL(0:n,L.base=~`1`,Ls.path=c(~`2`,~`3`))+
                    # ddspL distinct layers
                    ddspL(0:n,Ls.path=c(~`2`,~`3`))+
                    # dnspL distinct layers
                    dnspL(0:n,L.base=~`1`,Ls.path=c(~`2`,~`3`))+
                    # despL base and path distinct
                    despL(0:n,L.base=~`1`,Ls.path=c(~`2`,~`2`))+
                    # ddspL base and path distinct
                    ddspL(0:n,Ls.path=c(~`2`,~`2`))+
                    # dnspL base and path distinct
                    dnspL(0:n,L.base=~`1`,Ls.path=c(~`2`,~`2`))+
                    # despL base and path same
                    despL(0:n,L.base=~`2`,Ls.path=c(~`2`,~`2`))+
                    # ddspL base and path same
                    ddspL(0:n,Ls.path=c(~`2`,~`2`))+
                    # dnspL base and path same
                    dnspL(0:n,L.base=~`2`,Ls.path=c(~`2`,~`2`))+
                    # despL distinct base and one layer
                    despL(0:n,L.base=~`2`,Ls.path=c(~`2`,~`3`))+
                    # ddspL distinct base and one layer
                    ddspL(0:n,Ls.path=c(~`2`,~`3`))+
                    # dnspL distinct base and one layer
                    dnspL(0:n,L.base=~`2`,Ls.path=c(~`2`,~`3`))+
                    # Geometrically weighted
                    # despL distinct layers
                    dgwespL(decay,fixed=TRUE,L.base=~`1`,Ls.path=c(~`2`,~`3`))+
                    # ddspL distinct layers
                    dgwdspL(decay,fixed=TRUE,Ls.path=c(~`2`,~`3`))+
                    # dnspL distinct layers
                    dgwnspL(decay,fixed=TRUE,L.base=~`1`,Ls.path=c(~`2`,~`3`))## +
                    ## # despL base and path distinct
                    ## dgwespL(decay,fixed=TRUE,L.base=~`1`,Ls.path=c(~`2`,~`2`))+
                    ## # ddspL base and path distinct
                    ## dgwdspL(decay,fixed=TRUE,Ls.path=c(~`2`,~`2`))+
                    ## # dnspL base and path distinct
                    ## dgwnspL(decay,fixed=TRUE,L.base=~`1`,Ls.path=c(~`2`,~`2`))+
                    ## # despL base and path same
                    ## dgwespL(decay,fixed=TRUE,L.base=~`2`,Ls.path=c(~`2`,~`2`))+
                    ## # ddspL base and path same
                    ## dgwdspL(decay,fixed=TRUE,Ls.path=c(~`2`,~`2`))+
                    ## # dnspL base and path same
                    ## dgwnspL(decay,fixed=TRUE,L.base=~`2`,Ls.path=c(~`2`,~`2`))+
                    ## # despL distinct base and one layer
                    ## dgwespL(decay,fixed=TRUE,L.base=~`2`,Ls.path=c(~`2`,~`3`))+
                    ## # ddspL distinct base and one layer
                    ## dgwdspL(decay,fixed=TRUE,Ls.path=c(~`2`,~`3`))+
                    ## # dnspL distinct base and one layer
                    ## dgwnspL(decay,fixed=TRUE,L.base=~`2`,Ls.path=c(~`2`,~`3`))
                 ,
                  control=ctrl,
                  nsim=200))

  stats <- sapply(sim,
                  function(nw){
                    n <- network.size(nw)/3
                    m <- as.matrix(nw)
                    m1 <- m[seq_len(n),seq_len(n)]
                    m2 <- m[seq_len(n)+n,seq_len(n)+n]
                    m3 <- m[seq_len(n)+n*2,seq_len(n)+n*2]
                    
                    c(
                      # espL distinct layers
                      espL(0:n,L.base=m1,Ls.path1=m2,Ls.path2=m3),
                      # dspL distinct layers
                      dspL(0:n,Ls.path1=m2,Ls.path2=m3),
                      # nspL distinct layers
                      nspL(0:n,L.base=m1,Ls.path1=m2,Ls.path2=m3),
                      # espL base and path distinct
                      espL(0:n,L.base=m1,Ls.path1=m2,Ls.path2=m2),
                      # dspL base and path distinct
                      dspL(0:n,Ls.path1=m2,Ls.path2=m2),
                      # nspL base and path distinct
                      nspL(0:n,L.base=m1,Ls.path1=m2,Ls.path2=m2),
                      # espL base and path same
                      espL(0:n,L.base=m2,Ls.path1=m2,Ls.path2=m2),
                      # dspL base and path same
                      dspL(0:n,Ls.path1=m2,Ls.path2=m2),
                      # nspL base and path same
                      nspL(0:n,L.base=m2,Ls.path1=m2,Ls.path2=m2),
                      # espL distinct base and one layer
                      espL(0:n,L.base=m2,Ls.path1=m2,Ls.path2=m3),
                      # dspL distinct base and one layer
                      dspL(0:n,Ls.path1=m2,Ls.path2=m3),
                      # nspL distinct base and one layer
                      nspL(0:n,L.base=m2,Ls.path1=m2,Ls.path2=m3),
                      # Geometrically weighted
                      # espL distinct layers
                      gwespL(decay,n,L.base=m1,Ls.path1=m2,Ls.path2=m3),
                      # dspL distinct layers
                      gwdspL(decay,n,Ls.path1=m2,Ls.path2=m3),
                      # nspL distinct layers
                      gwnspL(decay,n,L.base=m1,Ls.path1=m2,Ls.path2=m3)## ,
                      ## # espL base and path distinct
                      ## gwespL(decay,n,L.base=m1,Ls.path1=m2,Ls.path2=m2),
                      ## # dspL base and path distinct
                      ## gwdspL(decay,n,Ls.path1=m2,Ls.path2=m2),
                      ## # nspL base and path distinct
                      ## gwnspL(decay,n,L.base=m1,Ls.path1=m2,Ls.path2=m2),
                      ## # espL base and path same
                      ## gwespL(decay,n,L.base=m2,Ls.path1=m2,Ls.path2=m2),
                      ## # dspL base and path same
                      ## gwdspL(decay,n,Ls.path1=m2,Ls.path2=m2),
                      ## # nspL base and path same
                      ## gwnspL(decay,n,L.base=m2,Ls.path1=m2,Ls.path2=m2),
                      ## # espL distinct base and one layer
                      ## gwespL(decay,n,L.base=m2,Ls.path1=m2,Ls.path2=m3),
                      ## # dspL distinct base and one layer
                      ## gwdspL(decay,n,Ls.path1=m2,Ls.path2=m3),
                      ## # nspL distinct base and one layer
                      ## gwnspL(decay,n,L.base=m2,Ls.path1=m2,Ls.path2=m3)
                    )
                  }) %>% t()

  expect_equivalent(attr(sim,"stats"), stats)
})
}
