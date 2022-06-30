#  File R/InitWtErgmTerm.multinet.R in package ergm.multi, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2022 Statnet Commons
################################################################################

# Arguments and outputs are identical to the binary version, except for the C routine name.
InitWtErgmTerm..subnets <- function(...){
  # Rename the function to avoid the extra nesting level in the
  # diagnostic messages.
  f <- InitErgmTerm..subnets
  #' @importFrom utils modifyList
  modifyList(f(...), list(name="_wtsubnets"))
}

# Arguments and outputs are identical to the binary version, except for the C routine names.

#' @rdname N-ergmTerm
#' @usage
#' # valued: N(formula, lm=~1, subset=TRUE, weights=1, contrasts=NULL, offset=0, label=NULL)
InitWtErgmTerm.N <- function(...){
  # Rename the function to avoid the extra nesting level in the
  # diagnostic messages.
  f <- InitErgmTerm.N
  term <- f(...)
  term$name <- switch(term$name,
                      MultiNet = "wtMultiNet",
                      MultiNets = "wtMultiNets")
  term
}

#' @rdname ByNetDStats-ergmTerm
#' @usage
#' # valued: ByNetDStats(formula, subset=TRUE)
#' @noRd
InitWtErgmTerm.ByNetDStats <- function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula", "subset"),
                      vartypes = c("formula","formula,logical,numeric,expression,call"),
                      defaultvalues = list(NULL,TRUE),
                      required = c(TRUE,FALSE))
  # Rename the function to avoid the extra nesting level in the
  # diagnostic messages.
  f <- InitErgmTerm.ByNetDStats
  modifyList(f(...), list(name = "wtByNetDStats"))
}
