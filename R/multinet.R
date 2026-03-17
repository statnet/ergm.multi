#  File R/InitErgmTerm.multinet.R in package ergm.multi, part of the Statnet
#  suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################
#' A multinetwork network representation.
#'
#' A function for specifying the LHS of a multi-network
#' (a.k.a. multilevel) ERGM. Typically used in conjunction with the
#' [`N()`][N-ergmTerm] term operator.
#'
#' @param ... network specification, in one of two formats:
#' 
#'   1. An (optionally named) list of networks with same directedness and bipartedness (but possibly different sizes).
#'
#'   1. Several networks as (optionally named) arguments.
#'
#' @return A [`network`] object comprising the provided networks, with multinetwork metadata.
#'
#' @templateVar combiner Networks
#' @template combine_networks_readonly
#'
#' @seealso [`ergmTerm`] for specific terms.
#' @seealso `vignette("Goeyvaerts_reproduction")` for a demonstration.
#' @examples
#'
#' data(samplk)
#'
#' # Method 1: list of networks
#' monks <- Networks(list(samplk1, samplk2))
#' ergm(monks ~ N(~edges))
#'
#' # Method 2: networks as arguments
#' monks <- Networks(samplk1, samplk2)
#' ergm(monks ~ N(~edges))
#' 
#' @export
Networks <- function(...){
  args <- list(...)
  if(all(sapply(args, is, "network"))){
    nwl <- args
  }else if(is.list(args[[1]]) && all(sapply(args[[1]], is, "network"))){
    nwl <- args[[1]]
  }else stop("Unrecognized format for multinetwork specification. See help for information.")

  if(!.same_constraints(nwl, "constraints")) stop("Networks have differing constraint structures. This is not supported at this time.")
  if(!.same_constraints(nwl, "obs.constraints")) stop("Networks have differing observation processes. This is not supported at this time.")
  
  nw <- combine_networks(nwl, ignore.nattr = c(eval(formals(combine_networks)$ignore.nattr), "ergm"))

  nw %n% "ergm" <- combine_ergmlhs(nwl)
  nw <- add_con(nw, blockdiag_tl, nwl[[1]])
  nw %n% ".combiner" <- "Networks"

  nw
}

#' @rdname Networks
#' @description `unNetworks()` extracts the networks into a list.
#'
#' @param object a multinetwork network returned by `Networks()`
#'
#' @export
unNetworks <- function(object) {
  assert_combined_network(object, "Networks", FALSE)
  uncombine_network(object, populate = TRUE)
}


#' An `as_tibble` method for combined networks.
#'
#' A method to obtain a network attribute table from a
#' [`combined_networks`] object, falling back to the
#' [network::as_tibble.network()] if vertex or edge attributes are
#' required.
#'
#' @param x  a [`combined_networks`] (inheriting from [`network::network`]).
#' @param attrnames a list (or a selection index) for attributes to obtain; for combined networks, defaults to all.
#' @param unit whether to obtain edge, vertex, or network attributes.
#' @param ... additional arguments, currently passed to unlist()].
#'
#' @seealso [network::as_tibble.network()]
#' @export
as_tibble.combined_networks<-function(x,attrnames=(match.arg(unit)%in%c("vertices","networks")), ..., unit=c("edges", "vertices", "networks")){
  unit <- match.arg(unit)
  if(unit!="networks") return(NextMethod())

  al <- x %n% ".snattr"

  if(is.logical(attrnames) || is.numeric(attrnames)) attrnames <- na.omit(names(al)[attrnames])
  else intersect(attrnames, names(al))

  out <- al[attrnames] %>% lapply(simplify_simple, ...) %>% as_tibble()

  combiner <- x%n%".combiner"
  if (!is.null(i <- ergm.multi_combiner(combiner)[["id"]])) out[[i]] <- seq_len(nrow(out))
  if (!is.null(i <- ergm.multi_combiner(combiner)[["name"]])) out[[i]] <- names(x%n%".snl")

  out
}
