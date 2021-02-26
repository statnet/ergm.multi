#  File R/zzz.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################
#' @import ergm
#' @import network
#' @import statnet.common
#' @import stats
## #' @import rlang

.onAttach <- function(lib, pkg){
  #' @importFrom statnet.common statnetStartupMessage
  sm <- statnetStartupMessage("ergm.multi", c("statnet"), FALSE)
  if(!is.null(sm)){
    packageStartupMessage(sm)
  }
}

.onLoad <- function(lib, pkg){
  # . is used as a placeholder by stantet.common::NVL3().
  utils::globalVariables(".")
  options(ergm.eval.loglik=TRUE)

  .RegisterProposals()
}

.RegisterProposals <- function(){
  ergm_proposal_table("c", "Bernoulli", "|bd&blockdiag", 0, "random", "blockdiag")
  ergm_proposal_table("c", "Bernoulli", "|bd&blockdiag&sparse", 1, "TNT", "blockdiagTNT")
}

#' @useDynLib ergm.multi
.onUnload <- function(libpath){
  library.dynam.unload("ergm.multi",libpath)
}
