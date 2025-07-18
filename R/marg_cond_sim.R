#  File R/marg_cond_sim.R in package ergm.multi, part of the Statnet suite of
#  packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################
#' Calculate gofN()-style Pearson residuals for arbitrary statistics
#'
#' This function is to be considered experimental. Do NOT rely on
#' it. It may, eventually, be moved to `ergm`, perhaps integrated into
#' the `simulate` methods.
#'
#' @param object an [`ergm`] object.
#' @param nsim number of realizations.
#' @param obs.twostage,GOF,save_stats see [gofN()].
#' @param control a control list returned by [control.gofN.ergm()]; note that `nsim` and `obs.twostage` parameters in the control list are ignored in favor of those passed to the function directly.
#' @param negative_info how to handle the situation in which the constrained variance exceeds the unconstrained: the corresponding action will be taken.
#' @param ... additional arguments to [ergm_model()], [simulate.ergm()], and [summary.ergm_model()].
#'
#' @return an object of similar structure as that returned by [gofN()].
#' @export
marg_cond_sim <- function(object, nsim=1, obs.twostage=nsim/2, GOF=NULL, control=control.gofN.ergm(), save_stats=FALSE, negative_info=c("error","warning","message","ignore"), ...){
  check.control.class(c("gofN.ergm"), "marg_cond_sim")
  negative_info <- match.arg(negative_info)
  control$obs.twostage <- obs.twostage
  control$nsim <- nsim
  if(control$obs.twostage && nsim %% control$obs.twostage !=0) stop("Number of imputation networks specified by control$obs.twostage control parameter must divide the nsim control parameter evenly.")
  
  nw <- object$network

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

    sim.m.obs_settings <- simulate(object, monitor=NULL, observational=TRUE, nsim=control$nsim, control=copy_parallel_controls(control,control$obs.simulate), basis=nw, output="stats", ..., return.args="ergm_model")
  }else control$obs.twostage <- FALSE # Ignore two-stage setting if no observational process.

  sim.m_settings <- simulate(object, monitor=NULL, nsim=control$nsim, control=copy_parallel_controls(control,control$simulate), basis=nw, output="stats", ..., return.args="ergm_model")

  message("Constructing GOF model.")
  NVL(GOF) <- if(length(object$formula)==3) object$formula[-2] else object$formula
  gof.m <- ergm_model(GOF, nw=nw, ...)
  nmonitored <- nstats <- nparam(gof.m, canonical=TRUE)
  
  # Indices of monitored elements.
  monitored <- nparam(sim.m_settings$object, canonical=TRUE) + seq_len(nmonitored)

  # The two-stage sample, taken marginally, *is* an unconstrained
  # sample.
  if(!control$obs.twostage){
    message("Simulating unconstrained sample.")
    sim <- do.call(simulate, .update.list(sim.m_settings, list(monitor=gof.m)))
    sim <- sim[,monitored,drop=FALSE]
    SST <- Welford(nrow(sim), colMeans(sim), .col_var(sim))
    if(!save_stats) suppressWarnings(rm(sim))
  }

  # TODO: Make this adaptive: start with a small simulation,
  # increase on fail; or perhaps use a pilot sample.
  if(is.na(object)){

    # Construct a simulate.ergm_state() call list for constrained simulation.
    args <- .update.list(sim.m.obs_settings,
                         list(monitor=gof.m, nsim=control$nsim/control$obs.twostage,
                              return.args="ergm_state"))
    sim.s.obs_settings <- do.call(simulate, args)
    suppressWarnings(rm(sim.m.obs_settings, gof.m))

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
    sim.obs <- do.call(simulate, .update.list(sim.s.obs_settings,
                                              list(nsim=control$nsim, return.args=NULL)))
    suppressWarnings(rm(sim.s.obs_settings))
    sim.obs <- sim.obs[,monitored,drop=FALSE]
  }else{
    sim.obs <- summary(gof.m, object$network, ...)
    suppressWarnings(rm(gof.m))
  }
  message("Collating the simulations.")

  if(save_stats){
    stats <- sim
    stats.obs <- sim.obs
  }

  # Calculate variances for each network and statistic.
  v <- SST$vars
  vo <- if(control$obs.twostage) MV else .col_var(sim.obs)
  # If any statistic for the network has negative variance estimate, stop with an error.
  if(negative_info != "ignore"){
    remain <- v>0 & v-vo<=0
    if(any(remain)){
      msg <- paste0(sum(remain), " network statistics have bad simulations after permitted number of retries. Rerun with higher nsim= control parameter.")
      switch(negative_info,
             error = stop(msg),
             warning = warning(msg),
             message = message(msg))
    }
  }

  m <- SST$means
  mo <- colMeans(sim.obs)

  suppressWarnings(rm(sim, sim.obs))

  message("Summarizing.")
  l <- list()
  l$var <- ifelse(v>0, v, NA)
  l$var.obs <- ifelse(v>0, vo, NA)
  l$observed <- ifelse(v>0, mo, NA)
  l$fitted <- ifelse(v>0, m, NA)
  l$pearson <- ifelse(v>0, (mo-m)/sqrt(v-vo), NA)
  if(save_stats){
    l$stats <- stats
    l$stats.obs <- stats.obs
  }

  structure(l, nw=nw, control=control)
}
