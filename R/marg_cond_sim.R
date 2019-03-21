#' An experimental function for calculating gofN()-style Pearson residuals for arbitrary statistics.
#'
#' It should probably be moved to `ergm`, perhaps integrated into the `simulate` methods.
#' @export
marg_cond_sim <- function(object, nsim=1, obs.twostage=nsim/2, GOF=NULL, control=control.simulate.ergm(), ...){
  check.control.class(c("simulate.ergm"), "marg_cond_sim")
  if(obs.twostage && nsim %% obs.twostage !=0) stop("Number of imputation networks specified by obs.twostage control parameter must divide the nsim control parameter evenly.")
  
  nw <- object$network

  stats <- stats.obs <- NULL

  cl <- ergm.getCluster(control)
  nthreads <- nthreads(control) # Fix this, so as not to become confused inside a clusterCall().

  message("Constructing simulation model(s).")
  if(!is.null(object$constrained.obs)){

    # First, make sure that the sample size is a multiple of the number of threads.
    if(obs.twostage && obs.twostage%%nthreads!=0){
      obs.twostage.new <- ceiling(obs.twostage/nthreads)*nthreads
      message(sQuote("obs.twostage"), " must be a multiple of number of parallel threads; increasing it to ", obs.twostage.new, ".")
      nsim <- round(obs.twostage.new/obs.twostage*nsim)
      obs.twostage <- obs.twostage.new
    }

    sim.obs_settings <- simulate(object, monitor=NULL, observational=TRUE, nsim=nsim, control=set.control.class("control.simulate.ergm",control), basis=nw, output="stats", response = object$response, ..., do.sim=FALSE)
  }else obs.twostage <- FALSE # Ignore two-stage setting if no observational process.

  sim_settings <- simulate(object, monitor=NULL, nsim=nsim, control=set.control.class("control.simulate.ergm",control), basis=nw, output="stats", response = object$response, ..., do.sim=FALSE)

  message("Constructing GOF model.")
  NVL(GOF) <- if(length(object$formula)==3) object$formula[-2] else object$formula
  gof.m <- ergm_model(GOF, nw=nw, response = object$response, ...)
  nstats <- nparam(gof.m, canonical=TRUE)

  # The two-stage sample, taken marginally, *is* an unconstrained
  # sample.
  if(!obs.twostage){
    message("Simulating unconstrained sample.")
    sim <- do.call(simulate, .update.list(sim_settings, list(monitor=gof.m)))
    sim <- sim[,attr(sim, "monitored"),drop=FALSE]
  }
  
  # TODO: Make this adaptive: start with a small simulation,
  # increase on fail; or perhaps use a pilot sample.
  if(!is.null(object$constrained.obs)){
    sim <-
      if(obs.twostage){
        message("Simulating imputed networks.", appendLF=FALSE)
        sim.net <- sim_settings$basis
        genseries <- function(){
          sim <- list()
          for(i in seq_len(obs.twostage/nthreads)){
            args <- .update.list(sim_settings,
                                 list(
                                   basis=sim.net,
                                   control=.update.list(sim_settings$control,
                                                        list(parallel=0,MCMC.burnin=if(i==1)sim_settings$control$MCMC.burnin else sim_settings$control$MCMC.interval)),
                                   output="pending_update_network", nsim=1))
            sim.net <- do.call(simulate, args)
            args <- .update.list(sim.obs_settings,
                                 list(basis=sim.net, monitor=gof.m, nsim=nsim/obs.twostage,
                                      control=.update.list(sim.obs_settings$control,
                                                           list(parallel=0))))
            sim[[i]] <- do.call(simulate, args)
            sim[[i]] <- sim[[i]][,attr(sim[[i]], "monitored"),drop=FALSE]
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
                                              list(nsim=nsim,
                                                   monitor=gof.m)))
    sim.obs <- sim.obs[,attr(sim.obs,"monitored"),drop=FALSE]
  }else{
    sim.obs <- summary(gof.m, object$network, response = object$response, ...)
  }
    
  stats <- sim
  stats.obs <- sim.obs
  
  # Calculate variances for each statistic.
  v <- apply(stats, 2, var)
  vo <- if(obs.twostage) apply(stats, 2, .mean_var, obs.twostage) else apply(stats.obs, 2, var)
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
