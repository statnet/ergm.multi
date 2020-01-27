#' An experimental function for calculating gofN()-style Pearson residuals for arbitrary statistics.
#'
#' It should probably be moved to `ergm`, perhaps integrated into the `simulate` methods.
#' @export
marg_cond_sim <- function(object, nsim=1, obs.twostage=nsim/2, GOF=NULL, control=control.simulate.ergm(), ...){
  check.control.class(c("simulate.ergm"), "marg_cond_sim")
  control$obs.twostage <- obs.twostage
  control$nsim <- nsim
  if(control$obs.twostage && nsim %% control$obs.twostage !=0) stop("Number of imputation networks specified by control$obs.twostage control parameter must divide the nsim control parameter evenly.")
  
  nw <- object$network

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
  gof.m <- ergm_model(GOF, nw=nw, response = object$response, ...)
  nmonitored <- nstats <- nparam(gof.m, canonical=TRUE)
  
  # Indices of monitored elements.
  monitored <- nparam(sim.m_settings$object, canonical=TRUE) + seq_len(nmonitored)

  # The two-stage sample, taken marginally, *is* an unconstrained
  # sample.
  if(!control$obs.twostage){
    message("Simulating unconstrained sample.")
    sim <- do.call(simulate, .update.list(sim.m_settings, list(monitor=gof.m)))
    sim <- sim[,monitored,drop=FALSE]
  }

  # TODO: Make this adaptive: start with a small simulation,
  # increase on fail; or perhaps use a pilot sample.
  if(!is.null(object$constrained.obs)){

    # Construct a simulate.ergm_state() call list for constrained simulation.
    args <- .update.list(sim.m.obs_settings,
                         list(monitor=gof.m, nsim=control$nsim/control$obs.twostage,
                              do.sim=FALSE))
    sim.s.obs_settings <- do.call(simulate, args)
    rm(sim.m.obs_settings, gof.m); gc()

    if(control$obs.twostage){
      message("Simulating imputed networks.", appendLF=FALSE)
      # Construct a simulate.ergm_state() call list for unconstrained simulation.
      args <- .update.list(sim.m_settings, list(do.sim=FALSE))
      sim.s_settings <- do.call(simulate, args)
      rm(sim.m_settings, args); gc()

      #' @importFrom parallel clusterCall
      sim <- if(!is.null(cl)) unlist(clusterCall(cl, gen_obs_imputation_series, sim.s_settings, sim.s.obs_settings, control, nthreads, monitored), recursive=FALSE)
             else gen_obs_imputation_series(sim.s_settings, sim.s.obs_settings, control, nthreads, monitored)
      message("")
      sim <- do.call(rbind, sim)
      gc()
    }
    message("Simulating constrained sample.")
    sim.obs <- do.call(simulate, .update.list(sim.s.obs_settings,
                                              list(nsim=control$nsim)))
    sim.obs <- sim.obs[,monitored,drop=FALSE]
  }else{
    sim.obs <- summary(gof.m, object$network, response = object$response, ...)
  }
  
  stats <- sim
  stats.obs <- sim.obs
  
  # Calculate variances for each statistic.
  v <- apply(stats, 2, var)
  vo <- if(control$obs.twostage) apply(stats, 2, .mean_var, control$obs.twostage) else apply(stats.obs, 2, var)
  # If any statistic for the network has negative variance estimate, rerun it.
  if(any(v>0 & v-vo<=0)) stop("Some statistics networks have bad simulations. Rerun with higher nsim= control parameter.")

  message("Summarizing.")
  l <- list()
  l$var <- ifelse(v>0, v, NA)
  l$var.obs <- ifelse(v>0, vo, NA)
  l$observed <- ifelse(v>0, colMeans(stats.obs), NA)
  l$fitted <- ifelse(v>0, colMeans(stats), NA)
  l$pearson <- ifelse(v>0, (l$observed-l$fitted)/sqrt(l$var-l$var.obs), NA)
  l$stats <- stats
  l$stats.obs <- stats.obs
  l

  structure(l, nw=nw, control=control)
}
