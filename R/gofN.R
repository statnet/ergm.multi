.mean_var <- function(x, ng){
  .Call("mean_var_wrapper", x, length(x)/ng, PACKAGE="ergm.multi")
}

.update.list <- function(l, v){
  l[names(v)]<-v
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
#' @export
gofN <- function(object, GOF=NULL, subset=TRUE, control=control.gofN.ergm(), ...){
  check.control.class(c("gofN.ergm"), "gofN")
  if(control$obs.twostage && control$nsim %% control$obs.twostage !=0) stop("Number of imputation networks specified by obs.twostage control parameter must divide the nsim control parameter evenly.")
  
  nw <- object$network
  nnets <- length(unique(.peek_vattrv(nw, ".NetworkID")))

  if(is.numeric(subset)) subset <- unwhich(subset, nnets)
  subset <- rep(subset, length.out=nnets)
  remain <- subset

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

    sim.obs_settings <- simulate(object, monitor=NULL, observational=TRUE, nsim=control$nsim, control=set.control.class("control.simulate.ergm",control), basis=nw, output="stats", response = object$response, ..., do.sim=FALSE)
  }else control$obs.twostage <- FALSE # Ignore two-stage setting if no observational process.

  sim_settings <- simulate(object, monitor=NULL, nsim=control$nsim, control=set.control.class("control.simulate.ergm",control), basis=nw, output="stats", response = object$response, ..., do.sim=FALSE)

  prev.remain <- NULL

  for(attempt in seq_len(control$retry_bad_nets + 1)){
    if(attempt!=1) message(sum(remain), " networks (", paste(which(remain),collapse=", "), ") have bad simulations; rerunning.")
    message("Constructing GOF model.")
    NVL(GOF) <- if(length(object$formula)==3) object$formula[-2] else object$formula
    if(any(NVL(prev.remain,FALSE)!=remain))
      pernet.m <- ergm_model(~ByNetDStats(GOF, subset=remain), nw=nw, response = object$response, ...)
    prev.remain <- remain
    # Subset of GOF statistics (as opposed to those used to drive the simulation):
    gof.stats <- c(rep(FALSE, nparam(sim_settings$object, canonical=TRUE)), # sim_settings$object is the ergm_model.
                   rep(TRUE, nparam(pernet.m, canonical=TRUE)))

    nstats <- nparam(pernet.m, canonical=TRUE)/sum(remain)

    # TODO: Simulations can be rerun only on the networks in the subset.
    
    # The two-stage sample, taken marginally, *is* an unconstrained
    # sample.
    if(!control$obs.twostage){
      message("Simulating unconstrained sample.")
      sim <- do.call(simulate, .update.list(sim_settings, list(monitor=pernet.m)))
      sim <- sim[,gof.stats,drop=FALSE]
    }

    # TODO: Make this adaptive: start with a small simulation,
    # increase on fail; or perhaps use a pilot sample.
    if(!is.null(object$constrained.obs)){
      sim <-
        if(control$obs.twostage){
          message("Simulating imputed networks.", appendLF=FALSE)
          sim.net <- sim_settings$basis
          genseries <- function(){
            sim <- list()
            for(i in seq_len(control$obs.twostage/nthreads)){
              args <- .update.list(sim_settings,
                                 list(
                                   basis=sim.net,
                                   control=.update.list(sim_settings$control,
                                                      list(parallel=0,MCMC.burnin=if(i==1)sim_settings$control$MCMC.burnin else sim_settings$control$MCMC.interval)),
                                   output="pending_update_network", nsim=1))
              sim.net <- do.call(simulate, args)
              args <- .update.list(sim.obs_settings,
                                 list(basis=sim.net, monitor=pernet.m, nsim=control$nsim/control$obs.twostage,
                                      control=.update.list(sim.obs_settings$control,
                                                         list(parallel=0))))
              sim[[i]] <- do.call(simulate, args)
              sim[[i]] <- sim[[i]][,gof.stats,drop=FALSE]
              message(".", appendLF=FALSE)
            }
            sim
          }
          #' @importFrom parallel clusterCall
          sim <- if(!is.null(cl)) unlist(clusterCall(cl, genseries),recursive=FALSE) else genseries()
          message("")
          do.call(rbind, sim)
        }
      message("Simulating constrained sample.")
      sim.obs <- do.call(simulate, .update.list(sim.obs_settings,
                                                list(nsim=control$nsim,
                                                     monitor=pernet.m)))
      sim.obs <- sim.obs[,gof.stats,drop=FALSE]
    }else{
      sim.obs <- matrix(summary(pernet.m, object$network, response = object$response, ...), nrow(sim), ncol(sim), byrow=TRUE)
    }

    message("Collating the simulations.")
    cn <- colnames(sim)[seq_len(nstats)] %>% sub(".*?:","", .)
    
    statarray <- array(c(sim), c(control$nsim, nstats, sum(remain)))
    dimnames(statarray) <- list(Iterations=NULL, Statistic=cn, Network=NULL)
    statarray.obs <- array(c(sim.obs), c(control$nsim, nstats, sum(remain)))
    dimnames(statarray.obs) <- list(Iterations=NULL, Statistic=cn, Network=NULL)

    if(is.null(stats)){
      stats <- statarray
      stats.obs <- statarray.obs
    }else{
      stats[,,remain[subset]] <- statarray
      stats.obs[,,remain[subset]] <- statarray.obs
    }

    # Calculate variances for each network and statistic.
    v <- apply(statarray, 2:3, var)
    vo <- if(control$obs.twostage) apply(statarray, 2:3, .mean_var, control$obs.twostage) else apply(statarray.obs, 2:3, var)
    # If any statistic for the network has negative variance estimate, rerun it.
    remain[remain] <- apply(v>0 & v-vo<=0, 2, any)
    if(!any(remain)) break;
  }
  if(any(remain))
    stop(sum(remain), " networks (", paste(which(remain),collapse=", "), ") have bad simulations after permitted number of retries. Rerun with higher nsim= control parameter.")

  message("Summarizing.")
  o <- setNames(lapply(seq_along(cn), function(i){
    s <- stats[,i,]
    so <- stats.obs[,i,]
    #' @importFrom tibble lst
    v <- apply(s,2,var)
    l <- list()
    l$var <- ifelse(v>0, v, NA)
    l$var.obs <- ifelse(v>0, if(control$obs.twostage) apply(s, 2, .mean_var, control$obs.twostage) else apply(so, 2, var), NA)
    l$observed <- ifelse(v>0, colMeans(so), NA)
    l$fitted <- ifelse(v>0, colMeans(s), NA)
    l$pearson <- ifelse(v>0, (l$observed-l$fitted)/sqrt(l$var-l$var.obs), NA)
    if(control$save_stats){
      l$stats <- s
      l$stats.obs <- so
    }
    l
  }), cn)

  structure(o, nw=nw, subset=subset, control=control, class="gofN")
}

#' @describeIn gofN Extract a subset of statistics for which goodness-of-fit had been computed.
#' @param drop whether the indexing operator should drop attributes and return simply a list.
#' @export
`[.gofN` <- function(x, i, ..., drop = FALSE){
  y <- NextMethod("[")

  if(drop) y
  else{
    mostattributes(y) <- attributes(x)
    names(y) <- names(x)[i]
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
#' 
#' @importFrom grDevices dev.interactive devAskNewPage
#' @importFrom graphics abline panel.smooth plot
#' @importFrom methods is
#' @importFrom glue glue
#' @export
plot.gofN <- function(x, against=NULL, which=1:2, col=1, pch=1, cex=1, ..., ask = length(which)>1 && dev.interactive(TRUE), id.n=3, main="{type} for {sQuote(name)}"){
  if(ask){
    prev.ask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(prev.ask))
  }

  if(any(sapply(list(against, col, pch, cex),
                function(x) is.character(x) || is(x,"formula"))))
    nattrs <- get_multinet_nattr_tibble(attr(x,"nw"))[attr(x,"subset"),]
  
  againstname <- switch(class(against),
                        character = against,
                        formula = despace(deparse(against[[length(against)]])),
                        `NULL` = "Predicted value",
                        despace(deparse(substitute(against))))
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
      plot(NVL(againstval,summ$fitted), summ$pearson, col=col, pch=pch, cex=cex,..., main = glue(main, type="Residuals vs. Fitted"), xlab=againstname, ylab="Pearson residual",type="n")
      panel.smooth(NVL(againstval,summ$fitted), summ$pearson, col=col, pch=ifelse(ei, NA, pch), cex=cex, ...)
      if(any(ei)) text(NVL(againstval,summ$fitted)[ei], summ$pearson[ei], col=col[ei], label=seq_along(summ$pearson)[ei], cex=cex[ei], ...)
      abline(h=0, lty=3, col="gray")
    }
    
    if(2L %in% which){
      plot(NVL(againstval,summ$fitted), sqrt(abs(summ$pearson)), col=col, pch=pch, cex=cex,..., main = glue(main, type="Scale-location"), xlab=againstname, ylab=expression(sqrt(abs("Pearson residual"))), type="n")
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
#' @param retry_bad_nets This parameter only has an effect if the
#'   network has missing data or observational process. It gives the
#'   number of times to retry simulating networks for which the
#'   estimated constrained variance is higher than unconstrained. Note
#'   that setting it `>0` is likely to bias estimates: the simulation
#'   should instead be rerun with a larger `nsim`.
#'
#' @param save_stats If `TRUE`, save the simulated network statistics;
#'   defaults to `FALSE` to save memory and disk space.
#'
#' @param nsim Number of networks to be randomly drawn using Markov chain Monte
#' Carlo.  This sample of networks provides the basis for comparing the model
#' to the observed network.
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
                            retry_bad_nets=0,
                            save_stats=FALSE,
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
