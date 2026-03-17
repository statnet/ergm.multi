#  File R/zzz.R in package ergm.multi, part of the Statnet suite of packages
#  for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################
#' @import ergm
#' @import network
#' @import statnet.common
#' @import stats
#' @importFrom Rdpack reprompt
#' @importFrom rlang "%@%" "%@%<-"

.onAttach <- function(libname, pkgname){
  #' @importFrom statnet.common statnetStartupMessage
  sm <- statnetStartupMessage("ergm.multi", c("statnet"), FALSE)
  if(!is.null(sm)){
    packageStartupMessage(sm)
  }
}

.onLoad <- function(libname, pkgname){
  # . is used as a placeholder by stantet.common::NVL3().
  utils::globalVariables(c(".", ".fitted", ".against", ".rownames", ".pearson", ".weight", "bids", "bmap"))

  eval(COLLATE_ALL_MY_CONTROLS_EXPR)

  .RegisterKeywords()
  .RegisterCombiners()
}

## BEGIN boilerplate: should be kept in sync with statnet.common.
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
eval(UPDATE_MY_SCTRL_EXPR)
## END boilerplate: should be kept in sync with statnet.common.

.RegisterKeywords <- function() {
  ergm_keyword(name="layer-aware", short="layer", description="operates on multilayer network constructs", popular=TRUE, package="ergm.multi")
}

.RegisterCombiners <- function() {
  ergm.multi_combiner("Layer", c(constructor = "Layer()", construct = "multi-layer network", element = "layer", id = ".LayerID", name = ".LayerName"))
  ergm.multi_combiner("Networks", c(constructor = "Networks()", element = "network", id = ".NetworkID", name = ".NetworkName"))
}

#' @useDynLib ergm.multi
.onUnload <- function(libpath){
  library.dynam.unload("ergm.multi",libpath)
}
