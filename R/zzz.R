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

.onAttach <- function(libname, pkgname){
  #' @importFrom statnet.common statnetStartupMessage
  sm <- statnetStartupMessage("ergm.multi", c("statnet"), FALSE)
  if(!is.null(sm)){
    packageStartupMessage(sm)
  }
}

.onLoad <- function(libname, pkgname){
  # . is used as a placeholder by stantet.common::NVL3().
  utils::globalVariables(".")
  options(ergm.eval.loglik=TRUE)

  eval(COLLATE_ALL_MY_CONTROLS_EXPR)

  .RegisterProposals()
  .RegisterKeywords()
}

# TODO: Figure out some automatic way to keep this in sync with statnet.common.
#' @name snctrl
#'
#' @title Statnet Control
#'
#' @description A utility to facilitate argument completion of control lists, reexported from `statnet.common`.
#'
#' @section Currently recognised control parameters:
#' This list is updated as packages are loaded and unloaded.
#'
#' \Sexpr[results=rd,stage=render]{statnet.common::snctrl_names()}
#'
#' @seealso [statnet.common::snctrl()]
#' @docType import
NULL
#' @export
snctrl <- statnet.common::snctrl
## BEGIN boilerplate: should be kept in sync with statnet.common.


eval(UPDATE_MY_SCTRL_EXPR)

.RegisterProposals <- function(){
  ergm_proposal_table("c", "Bernoulli", "|bd|bdmax&blockdiag", 0, "random", "blockdiag")
  ergm_proposal_table("c", "Bernoulli", "|bd|bdmax&blockdiag&sparse", 1, "TNT", "blockdiagTNT")
}

.RegisterKeywords <- function() {
  ergm_keyword(name="layer-aware", short="layer", description="This term operates on multilayer network constructs", popular=TRUE, package="ergm.multi")
}

#' @useDynLib ergm.multi
.onUnload <- function(libpath){
  library.dynam.unload("ergm.multi",libpath)
}
