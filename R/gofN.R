#  File R/gofN.R in package ergm.multi, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2022 Statnet Commons
################################################################################
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

#' Linear model diagnostics for multinetwork linear models
#'
#' [gofN()] performs a simulation to obtain Pearson residuals for the
#' multivariate linear model for ERGM parameters, which can then be
#' used for a variety of diagnostics and diagnostic plots developed by
#' \insertCite{KrCo22t;textual}{ergm.multi}.
#'
#' @param object an [`ergm`] object.
#' @param x a `gofN` object.
#' @param GOF a one-sided [`ergm`] formula specifying network
#'   statistics whose goodness of fit to test, or [`NULL`]; if `NULL`,
#'   uses the original model.
#' @param subset argument for the [`N`][N-ergmTerm] term.
#' @param \dots additional arguments to functions ([simulate.ergm()]
#'   and [summary.ergm_model()]) for the constructor.
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
#' @seealso [plot.gofN()] and [autoplot.gofN()] for plotting `gofN`
#'   objects to make residual plots; [ergm::gof()] for single-network
#'   goodness-of-fit simulations in \CRANpkg{ergm}
#'
#' @references \insertAllCited{}
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
#' \donttest{
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
#' plot(fit.gof, against=~log(.fitted)) # Plot against transformed fitted values.
#' }
#'
#' ### If 'ggplot2' and 'ggrepel' are installed, illustrate the autoplot() method.
#' if(require("ggplot2") && requireNamespace("ggrepel")){
#'   autoplot(fit.gof)
#' }
#'
#' # Default is good enough in this case, but sometimes, we might want to set it higher. E.g.,
#' \dontrun{
#' fit.gof <- gofN(fit, GOF=~edges, control=control.gofN.ergm(nsim=400))
#' }
#' 
#' @export
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
  if(is.na(object)){

    # First, make sure that the sample size is a multiple of the number of threads.
    if(control$obs.twostage && control$obs.twostage%%nthreads!=0){
      obs.twostage.new <- ceiling(control$obs.twostage/nthreads)*nthreads
      message(sQuote("control$obs.twostage"), " must be a multiple of number of parallel threads; increasing it to ", obs.twostage.new, ".")
      control$nsim <- round(obs.twostage.new/control$obs.twostage*control$nsim)
      control$obs.twostage <- obs.twostage.new
    }

    sim.m.obs_settings <- simulate(object, monitor=NULL, observational=TRUE, nsim=control$nsim, control=control$obs.simulate, basis=nw, output="stats", ..., return.args="ergm_model")
  }else control$obs.twostage <- FALSE # Ignore two-stage setting if no observational process.

  sim.m_settings <- simulate(object, monitor=NULL, nsim=control$nsim, control=control$simulate, basis=nw, output="stats", ..., return.args="ergm_model")

  message("Constructing GOF model.")
  NVL(GOF) <- if(length(object$formula)==3) object$formula[-2] else object$formula
  pernet.m <- ergm_model(~ByNetDStats(GOF), nw=nw, ...)
  nmonitored <- nparam(pernet.m, canonical=TRUE)
  nstats <- nmonitored/nnets
  cn <- param_names(pernet.m)[seq_len(nstats)] %>% sub(".*?~","", .)

  # Indices of monitored elements.
  monitored <- nparam(sim.m_settings$object, canonical=TRUE) + seq_len(nmonitored)

  # The two-stage sample, taken marginally, *is* an unconstrained
  # sample.
  if(!control$obs.twostage){
    message("Simulating unconstrained sample.")
    args <- .update.list(sim.m_settings, list(monitor=pernet.m, return.args="ergm_state"))
    sim.s_settings <- do.call(simulate, args)
    sim <- sim_stats_piecemeal(sim.s_settings, monitored,  max_elts, save_stats=save_stats)
    rm(sim.s_settings)
    SST <- attr(sim, "SST")

    if(!save_stats) suppressWarnings(rm(sim))
  }

  # TODO: Make this adaptive: start with a small simulation,
  # increase on fail; or perhaps use a pilot sample.
  if(is.na(object)){

    # Construct a simulate.ergm_state() call list for constrained simulation.
    args <- .update.list(sim.m.obs_settings,
                         list(monitor=pernet.m, nsim=control$nsim/control$obs.twostage,
                              return.args="ergm_state"))
    sim.s.obs_settings <- do.call(simulate, args)
    suppressWarnings(rm(sim.m.obs_settings, pernet.m))

    if(control$obs.twostage){
      message("Simulating imputed networks.", appendLF=FALSE)
      # Construct a simulate.ergm_state() call list for unconstrained simulation.
      args <- .update.list(sim.m_settings, list(return.args="ergm_state"))
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
      
      SST <- Reduce(update, SST)
      MV <- colMeans(do.call(rbind, MV))
      message("")
    }
    suppressWarnings(rm(sim.s_settings))
    message("Simulating constrained sample.")
    sim.obs <- sim_stats_piecemeal(.update.list(sim.s.obs_settings,
                                                list(nsim=control$nsim, return.args=NULL)), monitored, max_elts, save_stats)
    rm(sim.s.obs_settings)
    SST.obs <- attr(sim.obs, "SST")
    if(!save_stats) rm(sim.obs)
  }else{
    SST.obs <- Welford(2, summary(pernet.m, object$network, ...), numeric(nparam(pernet.m, canonical=TRUE)))
    if(save_stats) sim.obs <- matrix(SST.obs$means, control$nsim, nparam(pernet.m, canonical=TRUE), byrow=TRUE)
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
  v <- SST$vars
  vo <- if(control$obs.twostage) MV else SST.obs$vars
  # If any statistic for the network has negative variance estimate, stop with an error.
  remain <- v>0 & v-vo<=0
  if(any(remain))
    stop(sum(remain), " network statistics have bad simulations after permitted number of retries. Rerun with higher nsim= control parameter.")

  m <- SST$means
  mo <- SST.obs$means

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
    l$var <- ifelse(v[i,]>0, v[i,], NA)
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

  SST <- Welford(nstats)

  sim.s_settings$control$MCMC.samplesize <- first_nsim
  o <- with(sim.s_settings, ergm_MCMC_sample(object, control, coef))
  state <- sim.s_settings$object <- o$networks
  sim1 <- as.matrix(o$stats)[,monitored,drop=FALSE]

  if(save_stats) sim[[1L]] <- sim1
  SST <- update(SST, sim1)
  rm(sim1)

  for(rerun in seq_len(nreruns)){
    sim.s_settings$control$MCMC.samplesize <- chunk_nsim
    sim.s_settings$control$MCMC.burnin <- sim.s_settings$control$MCMC.interval
    o <- with(sim.s_settings, ergm_MCMC_sample(object, control, coef))
    state <- sim.s_settings$object <- o$networks # Make sure ext.state and nw0 are reconciled.
    sim1 <- as.matrix(o$stats)[,monitored,drop=FALSE]

    if(save_stats) sim[[rerun+1L]] <- sim1
    SST <- update(SST, sim1)
    rm(sim1)
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
  SST <- Welford(nstats)

  for(i in seq_len(n_cond)){
    # First, simulate a realisation of the unconstrained network.
    o <- with(sim.s_settings, ergm_MCMC_sample(object, control, coef))
    state <- sim.s_settings$object <- update(o$networks[[1]]) # Make sure ext.state and nw0 are reconciled.
    if(i==1) sim.s_settings$control$MCMC.burnin <- sim.s_settings$control$MCMC.interval # After the first iteration, shorter intervals.

    # Next, simulate a realisation of the constrained network conditional on the unconstrained.
    sim.s.obs_settings$object <- update(sim.s.obs_settings$object, el = state$el, nw0 = state$nw0) # Replace network state with that of unconstrained sample.
    sim.s.obs_settings$object <- update(sim.s.obs_settings$object, stats = summary(sim.s.obs_settings$object)) # Find the current network statistic.

    sim1 <- with(sim.s.obs_settings, as.matrix(ergm_MCMC_sample(object, control, coef)$stats)[,monitored,drop=FALSE])
    if(save_stats) sim[[i]] <- sim1

    SST <- update(SST, sim1)
    MV <- MV + (.col_var(sim1) - MV)/i

    message(".", appendLF=FALSE)
  }
  structure(sim, SST=SST, MV=MV)
}

#' @describeIn gofN Extract a subset of statistics for which goodness-of-fit had been computed.
#' @param i for the indexing operator, index of statistics to be kept in the subset.
#' @param j for the indexing operator, index of networks to be kept in the subset.
#' @param drop whether the indexing operator should drop attributes and return simply a list.
#' @export
`[.gofN` <- function(x, i, j, ..., drop = FALSE){
  y <- if(!missing(i)) unclass(x)[i] else x

  if(!missing(j)) y <- lapply(y, function(.y) lapply(.y, `[`, j))

  if(drop) y
  else{
    mostattributes(y) <- attributes(x)
    if(!missing(j)) attr(y, "subset") <- replace(logical(length(attr(x, "subset"))),
                                                 which(attr(x, "subset"))[j], TRUE)
    # This is necessary because i might be a character index.
    names(y) <- setNames(nm=names(x))[i]
    y
  }
}

#' Plotting methods for [`gofN`], making residual and scale-location plots.
#'
#' The [plot()] method uses \R graphics.
#'
#' @param x a [`gofN`] object.
#' @param against what the residuals should be plotted against. Note that different methods use different formats: see Details. Categorical ([`factor`] and [`ordered`]) values are visualised using boxplots, with [`ordered`] values also adding a smoothing line like the quantitative. Defaults to the fitted values.
#' @param col,pch,cex,bg vector of values (wrapped in [I()]), network attribute, or a formula whose RHS gives an expression in terms of network attributes to plot against.
#' @param which which to plot (`1` for residuals plot, `2` for \eqn{\sqrt{|R_i|}}{sqrt(|R_i|)} scale plot, and `3` for normal quantile-quantile plot).
#' @param ... additional arguments to [plot()], [qqnorm()], and [qqline()], and others.
#' @param ask whether the user should be prompted between the plots.
#' @param id.n maximum number of extreme points to label explicitly.
#' @param id.label specification for how extreme points are to be labeled, defaulting to network's index in the combined network.
#' @param main a template for the plots' titles; these use [glue()]'s templating, with `{type}` replaced with the type of plot and `{name}` replaced with the statistic.
#' @param xlab horizontal axis label; defaults to a character representation of `against`.
#' @param ylim vertical range for the plots, interpreted as in [graphics::plot()]; can be specified as a list with 3 elements, giving the range for the corresponding plot according to the plot numbers for the `which=` argument, and can be used to ensure that, e.g., diagnostic plots for different models are on the same scale.
#' @param cex.id scaling factor for characters used to label extreme points; see [plot.lm()].
#'
#' @details For the `plot()` method, `against` and `id.label` can be vectors of values (enclosed in [I()] to be used as is), a character string identifying a network attribute, or a formula whose RHS gives an expression in terms of network attributes to plot against. The `against` formula may also contain a `.fitted` variable which will be substituted with the fitted values.
#' 
#' @importFrom grDevices dev.interactive devAskNewPage adjustcolor
#' @importFrom graphics abline panel.smooth plot text axis boxplot lines points
#' @importFrom methods is
#' @importFrom glue glue
#'
#' @seealso [gofN()] for examples, [plot.lm()], [graphics::plot()] for regression diagnostic plots and their parameters.
#' @export
plot.gofN <- function(x, against=NULL, which=1:2, col=1, pch=1, cex=1, bg=0, ..., ask = length(which)>1 && dev.interactive(TRUE), id.n=3, id.label=NULL, main="{type} for {sQuote(name)}", xlab=NULL, ylim=NULL, cex.id=0.75){
  if(ask){
    prev.ask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(prev.ask))
  }

  if(!is.list(ylim)) ylim <- list(ylim)
  ylim <- rep_len(ylim, 3)

  if(any(sapply(list(against, col, pch, cex, id.label),
                function(x) is.character(x) || is(x,"formula"))))
    nattrs <- as_tibble(attr(x,"nw"), unit="networks")[attr(x,"subset"),]
  
  xlab <- NVL(xlab,
              switch(class(against),
                     character = against,
                     formula = deparse(do.call(substitute, list(ult(against), list(.fitted=as.name("Fitted values")))),width.cutoff=500L),
                     `NULL` = "Fitted values",
                     despace(deparse(substitute(against),width.cutoff=500L))))

  np <- sum(attr(x,"subset"))
  for(gpar in c("col", "bg", "pch", "cex", "id.label")){
    a <- get(gpar)
    a <- switch(class(a),
                AsIs = a,
                character = nattrs[[a]],
                formula = eval(a[[length(a)]], envir = nattrs, enclos = environment(a)),
                a)
    a <- rep(a, length.out=np)
    assign(gpar, a)
  }

  add.multismooth <- function(x, y, w, col, n){
    for(c in unique(col)){
      csel <- col == c
      lofit <- suppressWarnings(loess(y~x, weights=w, subset=csel, ...))
      xnew <- seq(from=min(x[csel],na.rm=TRUE), to=max(x[csel],na.rm=TRUE), length.out=ceiling(n*diff(range(x[csel],na.rm=TRUE))/diff(range(x,na.rm=TRUE))))
      lines(xnew, predict(lofit, newdata=data.frame(x=xnew)), col=c)
    }
  }

  points.smooth <- function(x, y, w, col, bg, pch, cex, cex.id, main, xlab, ylab, ylim, col.smooth=col, id=NA, n=101, ...){
    id <- rep_len(id, length(x))
    toID <- !is.na(id)
    col <- rep_len(col, length(x))
    plot(x, y, pch=pch, col=ifelse(toID, NA, col), bg=bg, cex=cex, main=main, xlab=xlab, ylab=ylab, ylim=ylim, type="p", ...)
    if(any(toID)) text(x[toID], y[toID], col=col[toID], label=id[toID], cex=cex[toID]*cex.id, ...)

    add.multismooth(x, y, w, col.smooth, n)

    abline(h=0, lty=3, col="gray")
  }

  boxplot.smooth <- function(x, y, w, col, bg, pch, cex, cex.id, main, xlab, ylab, ylim, col.bp=unique(col), id=NA, n=101, ...){
    id <- rep_len(id, length(x))
    toID <- !is.na(id)
    col <- rep_len(col, length(x))
    colID <- match(col, col.bp)
    nc <- length(col.bp)
    nx <- length(unique(x))
    ys <- split(y, list(colID,x))
    ymap <- as.numeric(interaction(list(colID,x)))
    xs <- rep(seq_len(nx)-1, each=nc) * (nc*1.5) + rep(seq_len(nc), nx)
    boxplot(ys, varwidth=TRUE, border=rep_len(col.bp, length(ys)), lty=1, col=rep_len(adjustcolor(col.bp,alpha.f=0.5), length(ys)), pch=pch, at=xs, main=main, xlab=xlab, ylab=ylab, ylim=ylim, xaxt="n")
    axis(1, at=tapply(xs, rep(seq_len(nx), each=nc), mean), labels=levels(x))

    if(any(toID)){
      # Ugly hack: overplot the points with white before adding text.
      points(xs[ymap[toID]], y[toID], col="white", pch=19, cex=1.1)
      text(xs[ymap[toID]], y[toID], col=col[toID], label=id[toID], cex=cex[toID]*cex.id)
    }

    xscl <- tapply(xs, rep(seq_len(nx), each=nc), mean)[match(x, sort(unique(x)))]
    add.multismooth(xscl, y, w, col, n)

    abline(h=0, lty=3, col="gray")
  }

  for(name in names(x)){
    summ <- x[[name]]
    againstval <- switch(class(against),
                         AsIs = {nattrs$.fitted <- summ$fitted; against},
                         character = {nattrs$.fitted <- summ$fitted; nattrs[[against]]},
                         formula = {nattrs$.fitted <- summ$fitted; eval(ult(against), envir = nattrs, enclos = environment(against))},
                         against)

    nn <- sum(!is.na(summ$pearson))
    ez <- qnorm((nn+.5)/(nn+1)) # Extreme standard normal quantile appropriate to the sample size.
    ei <- !is.na(summ$pearson) & rank(-abs(summ$pearson), ties.method="min")<=id.n & abs(summ$pearson)>ez
    NVL(id.label) <- seq_along(summ$pearson)

    NVL(againstval) <- summ$fitted
    resid.plotter <- if(is.factor(againstval)) boxplot.smooth else points.smooth

    w <- 1/(summ$var - summ$var.obs)

    if(1L %in% which){
      resid.plotter(againstval, summ$pearson, w=w, col=col, bg=bg, pch=pch, cex=cex, cex.id=cex.id, id=ifelse(ei, id.label, NA), main = glue(main, type="Residuals vs. Fitted"), xlab=xlab, ylab="Std. Pearson resid.", ylim=ylim[[1L]], ...)
    }
    
    if(2L %in% which){
      resid.plotter(againstval, sqrt(abs(summ$pearson)), w=w, col=col, bg=bg, pch=pch, cex=cex, cex.id=cex.id, id=ifelse(ei, id.label, NA), main = glue(main, type="Scale-location"), xlab=xlab, ylab=expression(sqrt(abs("Std. Pearson resid."))), ylim=ylim[[2L]], ...)
    }

    if(3L %in% which){
      qqnorm(summ$pearson, col=col, pch=pch, cex=cex,..., main = glue(main, type="Normal Q-Q"), ylim=ylim[[3L]])
      qqline(summ$pearson, ...)
    }
    
  }
}


#' @rdname plot.gofN
#'
#' @description The [ggplot2::autoplot()] method uses \CRANpkg{ggplot2} and \CRANpkg{ggrepel}.
#'
#' @param mappings a named list of lists of mappings constructed by [ggplot2::aes()] overriding the defaults. See Details below.
#' @param geom_args a named list of lists of arguments overriding the defaults for the individual geoms. See Details below.
#'
#' @details For `autoplot.gofN()`, `against` and `id.label` are interpreted as
#'   expressions in terms of network attributes and values generated by
#'   [augment.gofN()], included `.fitted` for the fitted values.
#'
#' @section Customising `autoplot.gofN()`:
#'
#' `autoplot.gofN()` constructs the plots out of [ggplot2::ggplot()],
#' [ggplot2::geom_point()] (for numeric `against`), [ggplot2::geom_boxplot()] for
#' categorical or ordinal `against`), and [ggplot2::geom_smooth()] (for numeric
#' or ordinal `against`), and [ggrepel::geom_text_repel()]. Mappings and
#' arguments passed through `mappings` and `geom_args` override the
#' respective defaults. They may have elements `default` (for
#' `ggplot()`), `point` (for `geom_point()` and `geom_boxplot()`),
#' `smooth` (for `geom_smooth`), and `text` (for `geom_text_repel()`).
#'
#' @return `autoplot.gofN()` returns a list of `ggplot` objects that
#'   if printed render to diagnostic plots. If there is only one, the
#'   object itself is returned.
#'
#' @method autoplot gofN
#' @rawNamespace S3method(ggplot2::autoplot, gofN)
autoplot.gofN <- function(x, against=.fitted, which=1:2,
                          mappings=list(),
                          geom_args=list(),
                          id.n=3, id.label=NULL){
  if(!requireNamespace("ggplot2") || !requireNamespace("ggrepel")) stop(sQuote("autoplot()"), " method for ", sQuote("gofN"), " objects requires packages ", paste.and(sQuote(c("ggplot2", "ggrepel"))), ".")

  against <- substitute(against)
  id.label <- substitute(id.label)

  nattrs <- generics::augment(x)
  xlab <- deparse1(do.call(substitute, list(against, list(.fitted=as.name("Fitted values")))))

  plots <- list()

  for(name in names(x)){
    a <- nattrs[nattrs$.stat_name==name,]

    nn <- sum(!is.na(a$.pearson))
    ez <- qnorm((nn+.5)/(nn+1)) # Extreme standard normal quantile appropriate to the sample size.
    ei <- !is.na(a$.pearson) & rank(-abs(a$.pearson), ties.method="min")<=id.n & abs(a$.pearson)>ez

    a$.against <- eval(against, envir = a, enclos = parent.frame())
    a$.rownames <- NVL(eval(id.label, envir = a, enclos = parent.frame()), seq_along(ei))

    mode_call <- function(f, custom_map, custom_args, def_map=ggplot2::aes(), def_args=list())
      do.call(f, c(list(mapping=modifyList(def_map, as.list(custom_map))), modifyList(def_args, as.list(custom_args))))

    resid_plot <- function(mapping, ylab){
      o <- mode_call(ggplot2::ggplot, mappings$default, geom_args$default, mapping, list(data=a))
      if(is.factor(a$.against)){
        o <- o + mode_call(ggplot2::geom_boxplot, mappings$point, geom_args$point)
        if(is.ordered(a$.against)) o <- o +  mode_call(ggplot2::geom_smooth, mappings$smooth, geom_args$smooth, ggplot2::aes(x=as.numeric(.against)), list(se=FALSE))
      }
      else o <- o + mode_call(ggplot2::geom_point, mappings$point, geom_args$point) +
             mode_call(ggplot2::geom_smooth, mappings$smooth, geom_args$smooth, def_args=list(se=FALSE))
      o <- o + ggplot2::xlab(xlab) + ggplot2::ylab(ylab) +
        mode_call(ggrepel::geom_text_repel, mappings$text, geom_args$text, ggplot2::aes(label=ifelse(ei, .rownames, "")))
      o
    }

    if(1L %in% which){
      o <- resid_plot(ggplot2::aes(x=.against, y=.pearson, weight=.weight),
                      "Std. Pearson resid.")
      plots <- c(plots, list(o))
    }

    if(2L %in% which){
      o <- resid_plot(ggplot2::aes(x=.against, y=sqrt(abs(.pearson)), weight=.weight),
                      expression(sqrt(abs("Std. Pearson resid."))))
      plots <- c(plots, list(o))
    }
  }

  if(length(plots)==1) plots[[1]] else plots
}


#' @describeIn gofN a method for constructing a [`tibble`] of network attributes augmented with goodness of fit information. Columns include:\describe{
#' \item{network attributes}{the attributes of each of the networks}
#'
#' \item{`.stat_name`}{name of the simulated statistic}
#'
#' \item{`.stat_id`}{index of the simulated statistic in the `gofN` object}
#'
#' \item{`.network_id`}{index of the network in the networks for which `gofN` was run (excluding those not in the subset)}
#'
#' \item{`.fitted`}{predicted value for the statistic}
#'
#' \item{`.observed`}{either the observed (for completely observed networks) or the predicted conditional on observed (for partially observed networks) value of the statistic}
#'
#' \item{`.pearson`}{the standardised Pearson residual}
#'
#' \item{`.var`, `.var.obs`}{estimated unconditional and average conditional variance of the statistic}
#'
#' \item{`.weight`}{inverse of the variance of the residual}
#' }
#'
#' @examples
#'
#' ### If 'generics' is installed, illustrate the augment() method.
#' if(require("generics")){
#'   augment(fit.gof)
#' }
#'
#' @method augment gofN
#' @rawNamespace S3method(generics::augment, gofN)
augment.gofN <- function(x, ...){
  if(!requireNamespace("generics")) stop(sQuote("augment()"), " method for ", sQuote("gofN"), " objects requires package ", paste.and(sQuote(c("generics"))), ".")

  nattrs <- as_tibble(attr(x,"nw"), unit="networks")[attr(x,"subset"),]

  seq_along(x) %>%
    map(~dplyr::bind_cols(.stat_name = names(x)[.],
                          .stat_id = .,
                          .network_id = seq_len(nrow(nattrs)),
                          nattrs,
                          .fitted = x[[.]]$fitted,
                          .pearson = x[[.]]$pearson,
                          .observed = x[[.]]$observed,
                          .var = x[[.]]$var,
                          .var.obs = x[[.]]$var.obs,
                          .weight = 1/(x[[.]]$var - x[[.]]$var.obs))) %>%
    dplyr::bind_rows()
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
      nattrs <- as_tibble(attr(object,"nw"), unit="networks")[attr(object,"subset"),]
    
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
#' @param simulate,obs.simulate Control lists produced by
#'   [control.simulate.ergm()] or equivalent for unconstrained and
#'   constrained simulation, respectively. Parameters are inherited
#'   from the model fit and can be overridden here.
#'
#' @template control_MCMC_parallel
#' 
#' @description `control.gofN.ergm` (or its alias, `control.gofN`) is
#'   intended to be used with [gofN()] specifically and will "inherit"
#'   as many control parameters from [`ergm`] fit as possible().
#'  
#' @export control.gofN.ergm
control.gofN.ergm<-function(nsim=100,
                            obs.twostage=nsim/2,
                            array.max=128,

                            simulate = control.simulate.ergm(),
                            obs.simulate = control.simulate.ergm(),

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

  nattrs <- as_tibble(attr(data,"nw"), unit="networks")[attr(data,"subset"),]
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
    .weights <- data.all$.weights <- ws[[lhs]]
    lm(formula, data=data.all, weights=.weights, ...)
  }, ...) %>% set_names(lhsl)
}
