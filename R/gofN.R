.mean_var <- function(x, ng){
  .Call("mean_var_wrapper", x, length(x)/ng, PACKAGE="ergm.multi")
}

.col_var <- function(x){
  .Call("vars_wrapper", x, nrow(x), PACKAGE="ergm.multi")
}

.update.list <- function(l, v){
  l[names(v)]<-v
  l
}

# l must be a list with three elements (in this specific order):
# * running sample size
# * running mean
# * running sum squared deviations
Welford_update <- function(l, x){
  if(is.numeric(x)){ # Either a vector or a matrix with statistics in rows.
    xm <- rbind(x)
    for(r in seq_len(nrow(xm))){
      x <- xm[r,]
      n.prev <- l[[1]]
      l[[1]] <- n.new <- n.prev + 1
      m.prev <- l[[2]]
      l[[2]] <- m.new <- m.prev + (x-m.prev)/n.new
      l[[3]] <- l[[3]] + (x-m.prev)*(x-m.new)
    }
  }else{ # Multielement update: both l and x are lists
    l.n <- l[[1]]; x.n <- x[[1]]; l[[1]] <- n.new <- l.n + x.n
    l.m <- l[[2]]; x.m <- x[[2]]
    d <- x.m - l.m
    # In our application, n for x and n for l are going to be similar, so we use weighted average.
    l[[2]] <- (l.n*l.m + x.n*x.m)/n.new
    l[[3]] <- l[[3]] + x[[3]] + d*d*l.n*x.n/n.new
  }
  l
}

#' Linear model diagnostics for multinetwork linear models
#'
#' @param object an [`ergm`] object.
#' @param GOF a one-sided [`ergm`] formula specifying network
#'   statistics whose goodness of fit to test, or [`NULL`]; if `NULL`,
#'   uses the original model.
#' @param subset argument for the [`N`][ergm-terms] term.
#' @param \dots additional arguments to functions ([simulate.ergm()]
#'   and [summary.ergm_model()] for the constructor, [plot()],
#'   [qqnorm()], and [qqline()] for the plotting method).
#' @param control See [control.gofN.ergm()].
#' @param save_stats If `TRUE`, save the simulated network statistics;
#'   defaults to `FALSE` to save memory and disk space.
#'
#' @return An object of class `gofN`: a named list containing a list
#'   for every statistic in the specified `GOF` formula with the
#'   following elements vectors of length equal to the number of
#'   subnetworks:
#'
#' \item{observed}{For completely observed networks, their value of
#' the statistic. For partially observed networks, the expected value
#' of their imputations under the model.}
#'
#' \item{fitted}{Expected value of the statistic under the model.}
#'
#' \item{var}{Variance of the statistic under the model.}
#'
#' \item{var.obs}{Conditional variance under imputation statistic.}
#' 
#' \item{pearson}{The Pearson residual computed from the above.}
#'
#' \item{stats,stats.obs}{If `save_stats` control parameter is `TRUE`, the simulated statistics.}
#' 
#' In addition, the following [`attr`]-style attributes are included:
#'
#' \item{nw}{The observed multinetwork object.}
#' 
#' \item{subset}{A logical vector giving the subset of networks that were used.}
#' 
#' \item{control}{Control parameters passed.}
#'
#' @examples
#' data(samplk)
#' monks <- Networks(samplk1, samplk2, samplk3,samplk1, samplk2, samplk3,samplk1, samplk2, samplk3)
#' fit <- ergm(monks~N(~edges+nodematch("group")))
#' fit.gof <- gofN(fit) # GOF = original model
#' summary(fit.gof)
#' plot(fit.gof)
#' fit.gof <- gofN(fit, GOF=~triangles)
#' summary(fit.gof)
#' plot(fit.gof)
#'
#' samplk1[1,]<-NA
#' samplk2[,2]<-NA
#' monks <- Networks(samplk1, samplk2, samplk3,samplk1, samplk2, samplk3,samplk1, samplk2, samplk3)
#' fit <- ergm(monks~N(~edges+nodematch("group")))
#' fit.gof <- gofN(fit) # GOF = original model
#' summary(fit.gof)
#' plot(fit.gof)
#' fit.gof <- gofN(fit, GOF=~triangles)
#' summary(fit.gof)
#' plot(fit.gof)
#' 
#' # Default is good enough in this case, but sometimes, we might want to set it higher. E.g.,
#' \dontrun{
#' fit.gof <- gofN(fit, GOF=~edges, control=control.gofN.ergm(nsim=400))
#' }
#' 
#' @exportc
gofN <- function(object, GOF=NULL, subset=TRUE, control=control.gofN.ergm(), save_stats=FALSE, ...){
  check.control.class(c("gofN.ergm"), "gofN")
  if(control$obs.twostage && control$nsim %% control$obs.twostage !=0) stop("Number of imputation networks specified by obs.twostage control parameter must divide the nsim control parameter evenly.")
  max_elts <- if(save_stats) Inf else control$array.max*1024^2/8 # A numeric is 8 bytes per element.

  nw <- object$network
  nnets <- length(unique(.peek_vattrv(nw, ".NetworkID")))

  if(is.numeric(subset)) subset <- unwhich(subset, nnets)
  subset <- rep(subset, length.out=nnets)

  stats <- stats.obs <- NULL

  cl <- ergm.getCluster(control)
  nthreads <- nthreads(control) # Fix this, so as not to become confused inside a clusterCall().

  message("Constructing simulation model(s).")
  if(!is.null(object$constrained.obs)){

    # First, make sure that the sample size is a multiple of the number of threads.
    if(control$obs.twostage && control$obs.twostage%%nthreads!=0){
      obs.twostage.new <- ceiling(control$obs.twostage/nthreads)*nthreads
      message(sQuote("control$obs.twostage"), " must be a multiple of number of parallel threads; increasing it to ", obs.twostage.new, ".")
      control$nsim <- round(obs.twostage.new/control$obs.twostage*control$nsim)
      control$obs.twostage <- obs.twostage.new
    }

    sim.m.obs_settings <- simulate(object, monitor=NULL, observational=TRUE, nsim=control$nsim, control=set.control.class("control.simulate.ergm",control), basis=nw, output="stats", response = object$response, ..., do.sim=FALSE)
  }else control$obs.twostage <- FALSE # Ignore two-stage setting if no observational process.

  sim.m_settings <- simulate(object, monitor=NULL, nsim=control$nsim, control=set.control.class("control.simulate.ergm",control), basis=nw, output="stats", response = object$response, ..., do.sim=FALSE)

  message("Constructing GOF model.")
  NVL(GOF) <- if(length(object$formula)==3) object$formula[-2] else object$formula
  pernet.m <- ergm_model(~ByNetDStats(GOF), nw=nw, response = object$response, ...)
  nmonitored <- nparam(pernet.m, canonical=TRUE)
  nstats <- nmonitored/nnets
  cn <- param_names(pernet.m)[seq_len(nstats)] %>% sub(".*?:","", .)

  # Indices of monitored elements.
  monitored <- nparam(sim.m_settings$object, canonical=TRUE) + seq_len(nmonitored)

  # The two-stage sample, taken marginally, *is* an unconstrained
  # sample.
  if(!control$obs.twostage){
    message("Simulating unconstrained sample.")
    args <- .update.list(sim.m_settings, list(monitor=pernet.m, do.sim=FALSE))
    sim.s_settings <- do.call(simulate, args)
    sim <- sim_stats_piecemeal(sim.s_settings, monitored,  max_elts, save_stats=save_stats)
    rm(sim.s_settings)
    SST <- attr(sim, "SST")

    if(!save_stats) suppressWarnings(rm(sim))
  }

  # TODO: Make this adaptive: start with a small simulation,
  # increase on fail; or perhaps use a pilot sample.
  if(!is.null(object$constrained.obs)){

    # Construct a simulate.ergm_state() call list for constrained simulation.
    args <- .update.list(sim.m.obs_settings,
                         list(monitor=pernet.m, nsim=control$nsim/control$obs.twostage,
                              do.sim=FALSE))
    sim.s.obs_settings <- do.call(simulate, args)
    suppressWarnings(rm(sim.m.obs_settings, pernet.m))

    if(control$obs.twostage){
      message("Simulating imputed networks.", appendLF=FALSE)
      # Construct a simulate.ergm_state() call list for unconstrained simulation.
      args <- .update.list(sim.m_settings, list(do.sim=FALSE))
      sim.s_settings <- do.call(simulate, args)
      suppressWarnings(rm(sim.m_settings, args))

      #' @importFrom parallel clusterCall
      if(!is.null(cl)){
        sim <- clusterCall(cl, gen_obs_imputation_series, sim.s_settings, sim.s.obs_settings, control, nthreads, monitored, save_stats)
        SST <- lapply(sim, attr, "SST")
        MV <- lapply(sim, attr, "MV")
        sim <- unlist(sim, recursive=FALSE)
      }else{
        sim <- gen_obs_imputation_series(sim.s_settings, sim.s.obs_settings, control, nthreads, monitored, save_stats)
        SST <- list(attr(sim, "SST"))
        MV <- list(attr(sim, "MV"))
      }
      message("")
      if(save_stats) sim <- do.call(rbind, sim) else rm(sim)
      
      SST <- Reduce(Welford_update, SST)
      MV <- colMeans(do.call(rbind, MV))
      message("")
    }
    suppressWarnings(rm(sim.s_settings))
    message("Simulating constrained sample.")
    sim.obs <- sim_stats_piecemeal(.update.list(sim.s.obs_settings,
                                                list(nsim=control$nsim)), monitored, max_elts, save_stats)
    rm(sim.s.obs_settings)
    SST.obs <- attr(sim.obs, "SST")
    if(!save_stats) rm(sim.obs)
  }else{
    SST.obs <- list(0, summary(pernet.m, object$network, response = object$response, ...))
    SST.obs[[3]] <- numeric(length(SST.obs[[2]]))
    if(save_stats) sim.obs <- matrix(SST.obs[[2]], control$nsim, nparam(pernet.m, canonical=TRUE), byrow=TRUE)
    suppressWarnings(rm(pernet.m))
  }
  message("Collating the simulations.")

  if(save_stats){
    stats <- array(c(sim), c(control$nsim, nstats, nnets))
    dimnames(stats) <- list(Iterations=NULL, Statistic=cn, Network=NULL)
    stats.obs <- array(c(sim.obs), c(control$nsim, nstats, nnets))
    dimnames(stats.obs) <- list(Iterations=NULL, Statistic=cn, Network=NULL)
  }

  # Calculate variances for each network and statistic.
  v <- SST[[3]]/(SST[[1]]-1)
  vo <- if(control$obs.twostage) MV else SST.obs[[3]]/(SST.obs[[1]]-1)
  # If any statistic for the network has negative variance estimate, stop with an error.
  remain <- any(v>0 & v-vo<=0)
  if(any(remain))
    stop(sum(remain), " network statistics have bad simulations after permitted number of retries. Rerun with higher nsim= control parameter.")

  m <- SST[[2]]
  mo <- SST.obs[[2]]

  suppressWarnings(rm(sim, sim.obs))

  # Reshape into matrices:
  v <- matrix(v, nstats, nnets, dimnames=list(Statistic=cn,Network=NULL))
  vo <- matrix(vo, nstats, nnets, dimnames=list(Statistic=cn,Network=NULL))
  m <- matrix(m, nstats, nnets, dimnames=list(Statistic=cn,Network=NULL))
  mo <- matrix(mo, nstats, nnets, dimnames=list(Statistic=cn,Network=NULL))

  message("Summarizing.")
  o <- setNames(lapply(seq_len(nstats), function(i){
    #' @importFrom tibble lst
    l <- list()
    l$var <- ifelse(v[i,]>0, v[,i], NA)
    l$var.obs <- ifelse(v[i,]>0, vo[i,], NA)
    l$observed <- ifelse(v[i,]>0, mo[i,], NA)
    l$fitted <- ifelse(v[i,]>0, m[i,], NA)
    l$pearson <- ifelse(v[i,]>0, (mo[i,]-m[i,])/sqrt(v[i,]-vo[i,]), NA)
    if(save_stats){
      s <- stats[,i,]
      so <- stats.obs[,i,]
      l$stats <- s
      l$stats.obs <- so
    }
    l
  }), cn)

  structure(o, nw=nw, subset=subset, control=control, class="gofN")
}

# Helper functions for gofN.
sim_stats_piecemeal <- function(sim.s_settings, monitored, max_elts, save_stats=FALSE){
  nthreads <- nthreads(sim.s_settings$control)
  # This needs to be a multiple of nthreads:
  chunk_nsim <- if(is.finite(max_elts)) ((max_elts %/% nparam(sim.s_settings$object, canonical=TRUE)) %/% nthreads)*nthreads else sim.s_settings$nsim
  nreruns <- ceiling(sim.s_settings$nsim / chunk_nsim) - 1
  first_nsim <- sim.s_settings$nsim - nreruns*chunk_nsim

  sim <- if(save_stats) vector("list", nreruns+1) else list()

  nstats <- switch(mode(monitored), # Vector of the correct length.
                    numeric = length(monitored),
                    logical = sum(monitored))

  # Welford's algorithm setup
  SST <- list(0L, numeric(nstats), numeric(nstats))

  sim.s_settings$control$MCMC.samplesize <- first_nsim
  o <- with(sim.s_settings, ergm_MCMC_sample(object, control, coef))
  state <- sim.s_settings$object <- o$networks
  sim1 <- as.matrix(o$stats)[,monitored,drop=FALSE]
  if(save_stats) sim[[1L]] <- sim1
  
  SST <- Welford_update(SST, sim1)

  for(rerun in seq_len(nreruns)){
    sim.s_settings$control$MCMC.samplesize <- chunk_nsim
    sim.s_settings$control$MCMC.burnin <- sim.s_settings$control$MCMC.interval
    o <- with(sim.s_settings, ergm_MCMC_sample(object, control, coef))
    state <- sim.s_settings$object <- o$networks # Make sure ext.state and nw0 are reconciled.
    sim1 <- as.matrix(o$stats)[,monitored,drop=FALSE]
    if(save_stats) sim[[rerun+1L]] <- sim1

    SST <- Welford_update(SST, sim1)
  }

  if(save_stats) sim <- do.call(rbind, sim)

  structure(sim, SST=SST, state=state)
}
gen_obs_imputation_series <- function(sim.s_settings, sim.s.obs_settings, control, nthreads, monitored, save_stats=FALSE){
  n_cond <- control$obs.twostage/nthreads
  sim <- if(save_stats) vector("list", n_cond) else list()

  sim.s_settings$control$parallel <- 0
  sim.s_settings$control$MCMC.samplesize <- 1

  sim.s.obs_settings$control$parallel <- 0
  sim.s.obs_settings$control$MCMC.samplesize <- control$nsim/control$obs.twostage

  nstats <- switch(mode(monitored), # Vector of the correct length.
                    numeric = length(monitored),
                    logical = sum(monitored))
  MV <- numeric(nstats)
  # Welford's algorithm setup
  SST <- list(0L, numeric(nstats), numeric(nstats))

  for(i in seq_len(n_cond)){
    # First, simulate a realisation of the unconstrained network.
    o <- with(sim.s_settings, ergm_MCMC_sample(object, control, coef))
    state <- sim.s_settings$object <- update(o$networks[[1]]) # Make sure ext.state and nw0 are reconciled.
    if(i==1) sim.s_settings$control$MCMC.burnin <- sim.s_settings$control$MCMC.interval # After the first iteration, shorter intervals.

    # Next, simulate a realisation of the constrained network conditional on the unconstrained.
    sim.s.obs_settings$object <- update(sim.s.obs_settings$object, el = state$el, nw0 = state$nw0) # Replace network state with that of unconstrained sample.

    sim1 <- with(sim.s.obs_settings, as.matrix(ergm_MCMC_sample(object, control, coef)$stats)[,monitored,drop=FALSE])
    if(save_stats) sim[[i]] <- sim1

    SST <- Welford_update(SST, sim1)
    MV <- MV + (.col_var(sim1) - MV)/i

    message(".", appendLF=FALSE)
  }
  structure(sim, SST=SST, MV=MV)
}

#' @describeIn gofN Extract a subset of statistics for which goodness-of-fit had been computed.
#' @param drop whether the indexing operator should drop attributes and return simply a list.
#' @export
`[.gofN` <- function(x, i, ..., drop = FALSE){
  y <- NextMethod("[")

  if(drop) y
  else{
    mostattributes(y) <- attributes(x)
    # This is necessary because i might be a character index.
    names(y) <- setNames(nm=names(x))[i]
    y
  }
}

#' @describeIn gofN A plotting method, making residual and scale-location plots.
#'
#' @param x a [`gofN`] object.
#' @param against vector of values, network attribute, or a formula whose RHS gives an expression in terms of network attributes to plot against; if `NULL` (default), plots against fitted values.
#' @param col,pch,cex vector of values (wrapped in [I()]), network attribute, or a formula whose RHS gives an expression in terms of network attributes to plot against.
#' @param which which to plot (`1` for residuals plot, `2` for \eqn{\sqrt{|R_i|}}{sqrt(|R_i|)} scale plot, and `3` for normal quantile-quantile plot).
#' @param ask whether the user should be prompted between the plots.
#' @param id.n Number of extreme points to label explicitly.
#' @param main A template for the plots' titles; these use [glue()]'s templating, with `{type}` replaced with the type of plot and `{name}` replaced with the statistic.
#' @param xlab Horizontal axis label; defaults to a character representation of `against`.
#' 
#' @importFrom grDevices dev.interactive devAskNewPage
#' @importFrom graphics abline panel.smooth plot text
#' @importFrom methods is
#' @importFrom glue glue
#' @export
plot.gofN <- function(x, against=NULL, which=1:2, col=1, pch=1, cex=1, ..., ask = length(which)>1 && dev.interactive(TRUE), id.n=3, main="{type} for {sQuote(name)}", xlab=NULL){
  if(ask){
    prev.ask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(prev.ask))
  }

  if(any(sapply(list(against, col, pch, cex),
                function(x) is.character(x) || is(x,"formula"))))
    nattrs <- get_multinet_nattr_tibble(attr(x,"nw"))[attr(x,"subset"),]
  
  xlab <- NVL(xlab,
              switch(class(against),
                     character = against,
                     formula = despace(deparse(against[[length(against)]])),
                     `NULL` = "Predicted value",
                     despace(deparse(substitute(against)))))
  againstval <- switch(class(against),
                       character = nattrs[[against]],
                       formula = eval(against[[length(against)]], envir = nattrs, enclos = environment(against)),
                       against)

  for(gpar in c("col", "pch", "cex")){
    a <- get(gpar)
    a <- switch(class(a),
                AsIs = a,
                character = nattrs[[a]],
                formula = eval(a[[length(a)]], envir = nattrs, enclos = environment(a)),
                a)
    assign(gpar, a)
  }

  for(name in names(x)){
    summ <- x[[name]]
    
    nn <- sum(!is.na(summ$pearson))
    ez <- qnorm((nn+.5)/(nn+1)) # Extreme standard normal quantile appropriate to the sample size.
    ei <- !is.na(summ$pearson) & rank(-abs(summ$pearson), ties.method="min")<=id.n & abs(summ$pearson)>ez

    if(1L %in% which){
      plot(NVL(againstval,summ$fitted), summ$pearson, col=col, pch=pch, cex=cex,..., main = glue(main, type="Residuals vs. Fitted"), xlab=xlab, ylab="Pearson residual",type="n")
      panel.smooth(NVL(againstval,summ$fitted), summ$pearson, col=col, pch=ifelse(ei, NA, pch), cex=cex, ...)
      if(any(ei)) text(NVL(againstval,summ$fitted)[ei], summ$pearson[ei], col=col[ei], label=seq_along(summ$pearson)[ei], cex=cex[ei], ...)
      abline(h=0, lty=3, col="gray")
    }
    
    if(2L %in% which){
      plot(NVL(againstval,summ$fitted), sqrt(abs(summ$pearson)), col=col, pch=pch, cex=cex,..., main = glue(main, type="Scale-location"), xlab=xlab, ylab=expression(sqrt(abs("Pearson residual"))), type="n")
      panel.smooth(NVL(againstval,summ$fitted), sqrt(abs(summ$pearson)), col=col, pch=ifelse(ei, NA, pch), cex=cex, ...)
      if(any(ei)) text(NVL(againstval,summ$fitted)[ei], sqrt(abs(summ$pearson))[ei], col=col[ei], label=seq_along(summ$pearson)[ei], cex=cex[ei], ...)
      abline(h=0, lty=3, col="gray")
    }

    if(3L %in% which){
      qqnorm(summ$pearson, col=col, pch=pch, cex=cex,..., main = glue(main, type="Normal Q-Q"))
      qqline(summ$pearson, ...)
    }
    
  }
}

#' @describeIn gofN A simple summary function.
#' @param by a numeric or character vector, or a formula whose RHS gives an expression in terms of network attributes, used as a grouping variable for summarizing the values.
#' @export
summary.gofN <- function(object, by=NULL, ...){
  cns <- names(object)
  if(is.null(by)){
    list(`Observed/Imputed values` = object %>% map("observed") %>% as_tibble %>% summary,
         `Fitted values` = object %>% map("fitted") %>% as_tibble %>% summary,
         `Pearson residuals`  = object %>% map("pearson") %>% as_tibble %>% summary,
         `Variance of Pearson residuals` = object %>% map("pearson") %>% map(var,na.rm=TRUE),
         `Std. dev. of Pearson residuals` = object %>% map("pearson") %>% map(sd,na.rm=TRUE))
  }else{
    if(is(by,"formula"))
      nattrs <- get_multinet_nattr_tibble(attr(object,"nw"))[attr(object,"subset"),]
    
    byname <- switch(class(by),
                     formula = despace(deparse(by[[length(by)]])),
                     despace(deparse(substitute(by))))
    byval <- switch(class(by),
                    formula = eval(by[[length(by)]], envir = nattrs, enclos = environment(by)),
                    by)

    list(`Observed/Imputed values` = object %>% map("observed") %>% as_tibble %>% split(byval) %>% map(summary),
         `Fitted values` = object %>% map("fitted") %>% as_tibble %>% split(byval) %>% map(summary),
         `Pearson residuals`  = object %>% map("pearson") %>% as_tibble %>% split(byval) %>% map(summary),
         `Variance of Pearson residuals` = object %>% map("pearson") %>% as_tibble %>% split(byval) %>% map(~apply(.,2,var,na.rm=TRUE)),
         `Std. dev. of Pearson residuals` = object %>% map("pearson") %>% as_tibble %>% split(byval) %>% map(~apply(.,2,sd,na.rm=TRUE)))
  }
}

#' Auxiliary for Controlling Multinetwork ERGM Linear Goodness-of-Fit Evaluation
#'
#' Auxiliary function as user interface for fine-tuning ERGM Goodness-of-Fit
#' Evaluation.
#' 
#' @param obs.twostage Either `FALSE` or an integer. This parameter
#'   only has an effect if the network has missing data or
#'   observational process. For such networks, evaluating the Pearson
#'   residual requires simulating the expected value of the
#'   conditional variance under the observation process. If `FALSE`,
#'   the simulation is performed conditional on the observed
#'   network. However, a more accurate estimate can be obtained via a
#'   two-stage process: \enumerate{
#' 
#' \item Sample networks from the model without the observational
#' constraint.
#'
#' \item Conditional on each of those networks, sample with the
#' observational constraint, estimating the variance within each
#' sample and then averaging over the first-stage sample.
#'
#' }
#'
#' Then, `obs.twostage` specifies the number of unconstrained networks
#' to simulate from, which should divide the [control.gofN.ergm()]'s
#' `nsim` argument evenly.
#'
#' @param nsim Number of networks to be randomly drawn using Markov chain Monte
#' Carlo.  This sample of networks provides the basis for comparing the model
#' to the observed network.
#'
#' @param array.max Try to avoid creating arrays larger in size (in
#'   megabytes) than this. Is ignored if `save_stats` is passed.
#'
#' @param MCMC.burnin Number of proposals before any MCMC sampling is done. It
#' typically is set to a fairly large number.
#' @param MCMC.interval Number of proposals between sampled statistics.
#' @param MCMC.prop.weights Specifies the proposal distribution used in the
#' MCMC Metropolis-Hastings algorithm.  Possible choices are \code{"TNT"} or
#' \code{"random"}; the \code{"default"} is one of these two, depending on the
#' constraints in place (as defined by the \code{constraints} argument of the
#' \code{\link{ergm}} function), though not all weights may be used with all
#' constraints.  The \code{TNT} (tie / no tie) option puts roughly equal weight
#' on selecting a dyad with or without a tie as a candidate for toggling,
#' whereas the \code{random} option puts equal weight on all possible dyads,
#' though the interpretation of \code{random} may change according to the
#' constraints in place.  When no constraints are in place, the default is TNT,
#' which appears to improve Markov chain mixing particularly for networks with
#' a low edge density, as is typical of many realistic social networks.
#' @param MCMC.prop.args An alternative, direct way of specifying additional
#' arguments to proposal.
#' @param MCMC.init.maxedges Maximum number of edges expected in network.
#' @param MCMC.runtime.traceplot Logical: If TRUE, plot traceplots of the MCMC
#' sample after every MCMC MLE iteration.
#' @param network.output R class with which to output networks. The options are
#' "network" (default) and "edgelist.compressed" (which saves space but only
#' supports networks without vertex attributes)
#' @template control_MCMC_parallel
#' @template seed
#' @template control_MCMC_packagenames
#' 
#' @description `control.gofN.ergm` (or its alias, `control.gofN`) is
#'   intended to be used with [gofN()] specifically and will "inherit"
#'   as many control parameters from [`ergm`] fit as possible().
#'  
#' @export control.gofN.ergm
control.gofN.ergm<-function(nsim=100,
                            obs.twostage=nsim/2,
                            array.max=128,

                       MCMC.burnin=NULL,
                       MCMC.interval=NULL,
                       MCMC.prop.weights=NULL,
                       MCMC.prop.args=NULL,
                       
                       MCMC.init.maxedges=NULL,
                       MCMC.packagenames=NULL,
                       
                       MCMC.runtime.traceplot=FALSE,
                       network.output="network",
                       
                       seed=NULL,
                       parallel=0,
                       parallel.type=NULL,
                       parallel.version.check=TRUE,
                       parallel.inherit.MT=FALSE){
  control<-list()
  for(arg in names(formals(sys.function())))
    control[arg]<-list(get(arg))

  set.control.class("control.gofN.ergm")
}

#' @rdname control.gofN.ergm
#' @export control.gofN
control.gofN <- control.gofN.ergm

#' Fit a linear model to the residuals in a gofN object.
#'
#' This non-method runs a properly weighted linear model on the raw
#' residuals of a [`gofN`] simulation for a multi-network ERGM fit.
#'
#' @param formula an [`lm`]-style formula. See Details for
#'   interpretation.
#'
#' @param data a [`gofN`] object.
#'
#' @param ... additional arguments to [lm()], excluding `weights`.
#'
#' @details The `formula`'s RHS is evaluated in an environment
#'   comprising the network statistics used in the [gofN()] call
#'   (which refer to the raw residuals for the corresponding
#'   statistic) and the network attributes.
#'
#'   The LHS is handled in a nonstandard manner, designed to make it
#'   easier to reference the usually lengthy network statistics:
#'   first, it is evaluated in the formula's environment. If the
#'   evaluation is successful and the result is numeric, these numbers
#'   are used as indices of the statistics in the [`gofN`] object to
#'   use on the RHS. If it is a character vector, it is treated as
#'   names of these statistics.
#'
#' @return A list of [`lm`] objects, one for each element of the
#'   vector on the LHS.
#'
#' @seealso [gofN()] and related methods.
#'
#' @examples
#' data(samplk)
#' # Add time indices:
#' samplk1 %n% "t" <- 1
#' samplk2 %n% "t" <- 2
#' samplk3 %n% "t" <- 3
#'
#' monks <- Networks(samplk1, samplk2, samplk3)
#'
#' fit <- ergm(monks~N(~edges+nodematch("group")))
#' fit.gof <- gofN(fit) # GOF = original model
#'
#' # Is there a time effect we should incorporate?
#' fit.gof.lm <- lm.gofN((1:2)~t, data=fit.gof)
#'
#' lapply(fit.gof.lm, summary)
#'
#' @importFrom statnet.common ERRVL eval_lhs.formula
#' @export lm.gofN
lm.gofN <- function(formula, data, ...){
  if(!is(data, "gofN")) stop("lm.gofN() requires a gofN object as the data= argument.")
  if("weights" %in% names(list(...))) stop("lm.gofN() does not accept a weights= argument at this time.")

  nattrs <- get_multinet_nattr_tibble(attr(data,"nw"))[attr(data,"subset"),]
  residuals <- data %>% map(~.$observed-.$fitted) %>% do.call(cbind, .)

  data.all <- cbind(residuals, nattrs)

  ws <- data %>% map(~.$var-.$var.obs) %>% map(`^`, -1)

  lhsl <- ERRVL(try({
    lhsl <- eval_lhs.formula(formula)
    if(is.numeric(lhsl)) lhsl <- names(data)[lhsl]
    as.character(lhsl)
  },silent=TRUE), as.character(formula[[2]]))

  lapply(lhsl, function(lhs, ...){
    formula[[2]] <- as.name(lhs)
    data.all$.weights <- ws[[lhs]]
    lm(formula, data=data.all, weights=.weights, ...)
  }, ...) %>% set_names(lhsl)
}
