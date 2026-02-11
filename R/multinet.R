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
  
  nw <- combine_networks(nwl, blockID.vattr=".NetworkID", blockName.vattr=".NetworkName", ignore.nattr = c(eval(formals(combine_networks)$ignore.nattr), "constraints", "obs.constraints", "ergm"), subnet.cache=TRUE)

  nw %n% "ergm" <- combine_ergmlhs(nwl)

  nw %ergmlhs% "constraints" <-
      if(NVL(nwl[[1]] %ergmlhs% "constraints",base_env(~.))==base_env(~.))
        base_env(~blockdiag(".NetworkID", noncontig = "split"))
      else
        append_rhs.formula(nwl[[1]]%ergmlhs%"constraints", list(call("blockdiag", ".NetworkID", noncontig = "split")), TRUE)
  if(!is.null(nwl[[1]]%ergmlhs%"obs.constraints")) nw %ergmlhs% "obs.constraints" <- nwl[[1]]%ergmlhs%"obs.constraints"

  nw
}

#' @rdname Networks
#' @description `unNetworks()` extracts the networks into a list.
#'
#' @param object a multinetwork network returned by `Networks()`
#'
#' @export
unNetworks <- function(object) {
  if (object %n% ".blockID.vattr" != ".NetworkID")
    stop("The specified network is not a sample of disjoint networks at the top level.")

  uncombine_network(object) |> map(ergmlhs_remove_blockdiag, ".NetworkID")
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
#' @template ergmTerm-NetworkIDName
#' @param store.nid whether to include columns with network ID and network name; the columns will be named with the arguments passed to `.NetworkID` and `.NetworkName`.
#' @param ... additional arguments, currently passed to unlist()].
#'
#' @seealso [network::as_tibble.network()]
#' @export
as_tibble.combined_networks<-function(x,attrnames=(match.arg(unit)%in%c("vertices","networks")), ..., unit=c("edges", "vertices", "networks"), .NetworkID=".NetworkID", .NetworkName=".NetworkName", store.nid=FALSE){
  unit <- match.arg(unit)
  if(unit!="networks") return(NextMethod())

  al <-
    if(is(x, "combined_networks") && !is.null(x %n% ".subnetattr")) (x %n% ".subnetattr")[[.NetworkID]]
    else{
      # FIXME: Probably more efficient to use attrnames earlier, in order to save calls to get.network.attribute().
      xl <- if(is.network(x)) subnetwork_templates(x, .NetworkID, .NetworkName) else x
      nattrs <- Reduce(union, lapply(xl, list.network.attributes))
      lapply(nattrs, function(nattr) lapply(xl, get.network.attribute, nattr)) %>% set_names(nattrs)
    }

  if(is.logical(attrnames) || is.numeric(attrnames)) attrnames <- na.omit(names(al)[attrnames])
  else intersect(attrnames, names(al))

  out <- al[attrnames] %>% lapply(simplify_simple, ...) %>% as_tibble()

  if(store.nid){
    out[[.NetworkID]] <- unique(x %v% .NetworkID)
    out[[.NetworkName]] <- unique(x %v% .NetworkName)
  }

  out
}

assert_LHS_Networks <- function(nw, nid, term_trace = TRUE, call = if(term_trace) NULL else rlang::caller_env()){
  if(anyNA(get_combining_attr(nw, nid))){
    msg <- paste0("The LHS of the model is not a multi-network ", sQuote("Networks()"), " construct.")
    if(term_trace) ergm_Init_abort(msg, call=call)
    else abort(msg, call=call)
  }
}
