#  File tests/testthat/test-multinet.R in package ergm.multi, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################
data(samplk)

logit <- function(p) log(p/(1-p))
times <- 1:3

samplk1%e%"w" <- floor(runif(network.edgecount(samplk1), 0, 4))
samplk2%e%"w" <- floor(runif(network.edgecount(samplk2), 0, 4))
samplk3%e%"w" <- floor(runif(network.edgecount(samplk3), 0, 4))
samplkl <- list(samplk1, samplk2, samplk3)
samplks <- Networks(samplk1, samplk2, samplk3)

test_that("N() summary with an LM and noncompacted stats", {
  summ.N <- summary(samplks~N(~edges+nodematch("cloisterville"), ~1+.NetworkID), term.options=list(N.compact_stats=FALSE))
  summ.l <- unlist(lapply(samplkl, function(nw) summary(~edges+nodematch("cloisterville"), basis = nw)))
  names(summ.l) <- paste0("N#", rep(1:3, each=2), "~", names(summ.l))
  expect_equal(summ.N, summ.l)
})


test_that("N() summary with offset and compacted stats", {
  summ.N <- summary(samplks~N(~edges, offset=~.NetworkID))
  summ.l <- sapply(samplkl, function(nw) summary(~edges, basis = nw))
  summ.l <- c(sum(summ.l), sum(summ.l*1:3))
  expect_equal(summ.N, summ.l, ignore_attr=TRUE)
})

test_that("Valued N() summary with an LM and noncompacted stats", {
  summ.N <- summary(samplks~N(~sum+nodematch("cloisterville", form="sum"), ~1+.NetworkID), term.options=list(N.compact_stats=FALSE), response="w")
  summ.l <- unlist(lapply(samplkl, function(nw) summary(~sum+nodematch("cloisterville", form="sum"), basis = nw, response="w")))
  expect_equal(summ.N, summ.l, ignore_attr=TRUE)
})


test_that("Valued N() summary with offset and compacted stats", {
  summ.N <- summary(samplks~N(~sum, offset=~.NetworkID), response="w")
  summ.l <- sapply(samplkl, function(nw) summary(~sum, basis = nw, response="w"))
  summ.l <- c(sum(summ.l), sum(summ.l*1:3))
  expect_equal(summ.N, summ.l, ignore_attr=TRUE)
})


pl <- lapply(samplkl, function(nw) ergmMPLE(~edges+nodematch("cloisterville"), basis = nw))

for(N.compact_stats in c(FALSE,TRUE)){
  testlab <- if(N.compact_stats) "compacted stats" else "noncompacted stats"
  
  test_that(paste("N() estimation with an LM and", testlab),{
    # Three networks, jointly.
    ergmfit <- ergm(samplks~N(~edges+nodematch("cloisterville"), ~1+.NetworkID), control=control.ergm(term.options=list(N.compact_stats=N.compact_stats)))
  
    nr <- sapply(lapply(pl, `[[`, "response"),length)
  
    y <- unlist(lapply(pl, `[[`, "response"))
    x <- do.call(rbind,lapply(pl, `[[`, "predictor"))
    x <- as.data.frame(cbind(x,.NetworkID=rep(times, nr)))
    w <- unlist(lapply(pl, `[[`, "weights"))
    glmfit <- glm(y~.NetworkID*nodematch.cloisterville,data=x,weights=w,family="binomial")
    expect_equal(coef(glmfit),coef(ergmfit), ignore_attr=TRUE, tolerance=1e-7)
  })
  
  test_that(paste("N() estimation with an LM, subsetting, and", testlab),{
    # Ignore second (test subset).
    ergmfit <- ergm(samplks~N(~edges+nodematch("cloisterville"), ~1+.NetworkID, subset=quote(.NetworkID!=2)), control=control.ergm(term.options=list(N.compact_stats=N.compact_stats)))
    
    pl2 <- pl[-2]
    times2 <- times[-2]
    
    nr <- sapply(lapply(pl2, `[[`, "response"),length)
    
    y <- unlist(lapply(pl2, `[[`, "response"))
    x <- do.call(rbind,lapply(pl2, `[[`, "predictor"))
    x <- as.data.frame(cbind(x,.NetworkID=rep(times2, nr)))
    w <- unlist(lapply(pl2, `[[`, "weights"))
    glmfit <- glm(y~.NetworkID*nodematch.cloisterville,data=x,weights=w,family="binomial")
    expect_equal(coef(glmfit),coef(ergmfit), ignore_attr=TRUE, tolerance=1e-7)
  })
    
  test_that(paste("N() estimation with an LM, subsetting down to only 1 network, and", testlab),{
    # Ignore first and third (test subset).
    expect_warning(ergmfit <- ergm(samplks~N(~edges+nodematch("cloisterville"), ~1+.NetworkID, subset=quote(.NetworkID==2)), control=control.ergm(term.options=list(N.compact_stats=N.compact_stats))), ".*model is nonidentifiable.*")
    
    pl2 <- pl[2]
    times2 <- times[2]
    
    nr <- sapply(lapply(pl2, `[[`, "response"),length)
    
    y <- unlist(lapply(pl2, `[[`, "response"))
    x <- do.call(rbind,lapply(pl2, `[[`, "predictor"))
    x <- as.data.frame(cbind(x,.NetworkID=rep(times2, nr)))
    w <- unlist(lapply(pl2, `[[`, "weights"))
    glmfit <- glm(y~.NetworkID*nodematch.cloisterville,data=x,weights=w,family="binomial")
    # Note that unlike glmift, ergm cannot detect nonidentifiability at
    # this time.
    if(N.compact_stats){
      expect_equal(na.omit(coef(glmfit)),na.omit(coef(ergmfit)), ignore_attr=TRUE)
    }else{
      expect_equal(na.omit(coef(glmfit)),
                        c(matrix(c(1,2,0,0,
                                   0,0,1,2),2,4,byrow=TRUE)%*%coef(ergmfit)),
                   ignore_attr=TRUE, tolerance=1e-7)
    }
  })
}


test_that("N() warns on parameter name mismatch",{
  samplk1%v%"x" <- 1:2
  samplk2%v%"x" <- 3:4
  expect_warning(summary(Networks(samplk1,samplk2)~N(~nodemix("x"))),
                 "different parameter names")
})


test_that("N() stops if given a non-multi-net object with a sensible error message.",{
  expect_error(ergm_model(samplk1 ~ N(~edges)),
               "In term 'N' in package 'ergm\\.multi': The LHS of the model is not a multi-network 'Networks\\(\\)' construct\\.")
})

test_that("gofN() stops if given a non-multi-net fit with a sensible error message.",{
  expect_error(gofN(ergm(samplk1 ~ edges)),
               "The LHS of the model is not a multi-network 'Networks\\(\\)' construct\\.")
})
