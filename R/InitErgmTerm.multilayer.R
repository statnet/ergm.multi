#  File R/InitErgmTerm.multilayer.R in package ergm.multi, part of the Statnet
#  suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################

## InitErgmTerm..layer.nets <- function(nw, arglist, ...){
##   a <- check.ErgmTerm(nw, arglist,
##                       varnames = c(),
##                       vartypes = c(),
##                       defaultvalues = list(),
##                       required = c())
##   list(name="_layer_nets", coef.names=c(), inputs=unlist(.block_vertexmap(nw, ".LayerID", TRUE)), dependence=FALSE)
## }


InitErgmTerm..layer.net <- function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("L"),
                      vartypes = c(ERGM_LAYER_SPEC),
                      defaultvalues = list(NULL),
                      required = c(TRUE))

  assert_LHS_Layer(nw)

  nwl <- subnetwork_templates(nw, ".LayerID", ".LayerName", copy.ergmlhs = c())

  L <- ergm_LayerLogic(a$L, nw)
  # Terms on this logical layer will induce dyadic independence if its
  # value depends on more than one other layer value.
  dependence <- length(L%@%"dep") > 1

  if (test_eval.LayerLogic(L, FALSE))
    ergm_Init_stop("Layer specification ", sQuote(deparse1(a$L)),
                   " outputs edges when all input layers are empty.",
                   " This is not supported at this time.", call. = FALSE)

  list(name = "_layer_net", coef.names = c(), dependence = dependence,
       iinputs = c(unlist(.block_vertexmap(nw, ".LayerID", TRUE)),
                   if (is.directed(nw)) sapply(nwl, function(nw) (nw%v%".undirected")[1]),
                   L%@%"C"))
}


#' @templateVar name L
#' @title Evaluation on layers
#' @description Evaluates the terms in `formula` on an observed or
#'   logical layers specified in formula `Ls` and sums the results
#'   elementwise.
#'
#' @usage
#' # binary: L(formula, Ls=~.)
#'
#' @template ergmTerm-formula
#' @templateVar Ls.interp , on which to evaluate `formula`
#' @template ergmTerm-Ls-1
#'
#' @note For the purposes of the terms in the `formula`, nonstandard
#'   network and vertex attributes will be taken from the *first*
#'   layer's network. The subsequent networks' attributes will be
#'   ignored by it. If `Ls` is specified as a [`function`], this
#'   function will be able to access them, however.
#'
#' @template ergmTerm-general
#'
#' @seealso [Layer()]
#' @concept operator
#' @concept layer-aware
InitErgmTerm.L <- function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula", "Ls"),
                      vartypes = c("formula", ERGM_LAYERS_SPEC),
                      defaultvalues = list(NULL, empty_env(~.)),
                      required = c(TRUE, FALSE))

  assert_LHS_Layer(nw)

  nwl <- subnetwork_templates(nw, ".LayerID", ".LayerName", copy.ergmlhs = c())

  Ls <- ergm_LayerLogics(a$Ls, nw)

  w <- rep(1, length(Ls))
  have.LHS <- lengths(Ls) == 3
  w[have.LHS] <- as.numeric(sapply(lapply(Ls[have.LHS], "[[", 2), eval, environment(Ls[[1]])))

  nw1 <- nwl[[1]]
  m <- ergm_model(a$formula, nw1, ..., offset.decorate=FALSE)

  ## FIXME: Is this consistent with extended state API, or do we need to have a different "model" for each layer?
  wm <- wrap.ergm_model(m, nw1, ergm_mk_std_op_namewrap("L", toString(Ls)))
  gs <- wm$emptynwstats
  wm$emptynwstats <- if(!is.null(gs)) gs * sum(w)
  wm$dependence <- wm$dependence || NA # If not determined by the model, set based on the layer logic.

  c(list(name = "OnLayer", iinputs = length(Ls), inputs = w,
         submodel = m, auxiliaries = Ls%@%"aux"),
    wm)
}

#' @templateVar name CMBL
#' @title Conway--Maxwell-Binomial dependence among layers
#' @description Models marginal dependence layers within each dyad by imposing
#'   a Conway--Maxwell-Binomial (CMB) distribution on the number of
#'   layers in each dyad that have a tie.
#'
#'   The term adds one statistic to the model, equalling the sum over
#'   all the dyads in the network of \eqn{\log\{E!(R-E)!/R!\}} , where
#'   \eqn{E} is the number of layers in `Ls` with an edge in that
#'   dyad and \eqn{R} being the total number of layers in `Ls` .
#'
#' @details A positive coefficient induces positive dependence and a negative
#'   one induces negative dependence.
#'
#' @usage
#' # binary: CMBL(Ls=~.)
#' @templateVar Ls.howmany at least two
#' @templateVar Ls.interp .
#' @template ergmTerm-Ls
#'
#' @template ergmTerm-general
#'
#' @concept directed
#' @concept undirected
#' @concept layer-aware
InitErgmTerm.CMBL <- function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("Ls"),
                      vartypes = c(ERGM_LAYERS_SPEC),
                      defaultvalues = list(empty_env(~.)),
                      required = c(FALSE))

  assert_LHS_Layer(nw)

  Ls <- ergm_LayerLogics(a$Ls, nw)

  list(name = "layerCMB", coef.names = paste0("CMBL(", Ls%@%"repr", ")"),
       iinputs = length(Ls), dependence = TRUE, auxiliaries = Ls%@%"aux")
}

################################################################################

#' @templateVar name twostarL
#' @title Multilayer two-star
#' @description This term adds one statistic to the model, equal to the number of
#'   cross-layer two-stars or two-paths in
#'   the network.
#'
#' @usage
#' # binary: twostarL(Ls, type, distinct=TRUE)
#' @templateVar Ls.howmany two
#' @templateVar Ls.interp {} specifying the layers of interest.
#' @template ergmTerm-Ls
#' @param type if the network is directed, an argument to determine which configurations are counted:
#'   \describe{
## #'   \item{`"any"`}{Number of configurations
## #'   \eqn{(i-j), (i-k)}{(i,j), (i,k)}, where
## #'   \eqn{(i-j)}{(i,j)} is in logical layer `Ls[[1]]`
## #'   and \eqn{(i-k)}{(i,k)} is in logical layer `Ls[[2]]`.}
#'   \item{`"out"`}{Number of configurations
#'   \eqn{(i{\rightarrow}j), (i{\rightarrow}k)}{(i,j), (i,k)}, where
#'   \eqn{(i{\rightarrow}j)}{(i,j)} is in logical layer `Ls[[1]]`
#'   and \eqn{(i{\rightarrow}k)}{(i,k)} is in logical layer `Ls[[2]]`.}
#'   \item{`"in"`}{Number of configurations
#'   \eqn{(j{\rightarrow}i), (k{\rightarrow}i)}{(j,i), (k,i)}, where
#'   \eqn{(j{\rightarrow}i)}{(j,i)} is in logical layer `Ls[[1]]`
#'   and \eqn{(k{\rightarrow}i)}{(k,i)} is in logical layer `Ls[[2]]`.}
#'   \item{`"path"`}{Number of configurations
#'   \eqn{(j{\rightarrow}i), (i{\rightarrow}k)}{(j,i), (i,k)}, where
#'   \eqn{(j{\rightarrow}i)}{(j,i)} is in logical layer `Ls[[1]]`
#'   and \eqn{(i{\rightarrow}k)}{(i,k)} is in logical layer `Ls[[2]]`.}
#'   }
#'
#'   This argument is ignored for undirected networks.
#'
#' @param distinct if `TRUE`, \eqn{j} and \eqn{k} above are required to
#'   be distinct. That is, the constituent edges may not be coincident or
#'   reciprocal.
#'
#' @template ergmTerm-general
#'
#' @concept directed
#' @concept undirected
#' @concept layer-aware
InitErgmTerm.twostarL<-function(nw, arglist,  ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("Ls", "type", "distinct"),
                      vartypes = c(ERGM_LAYERS_SPEC, "character", "logical"),
                      defaultvalues = list(NULL, NULL, TRUE),
                      required = c(TRUE, FALSE, FALSE))

  assert_LHS_Layer(nw)

  TYPES <- c("any", "out", "in", "path")
  TYPEREP <- setNames(c("--", "<>", "><", ">>"), TYPES)

  type <- a$type
  if (!is.directed(nw)) type <- "any"
  else if (is.null(type)) ergm_Init_stop(sQuote("type"), " argument is required for directed networks")

  type <- match.arg(tolower(type), TYPES)
  if (type == "any" && is.directed(nw)) ergm_Init_stop(paste0("at this time, ", sQuote('type="any"'), " is only supported for undirected networks"))
  typeID <- match(type, TYPES) - 1L

  Ls <- ergm_LayerLogics(a$Ls, nw) |> rep_len(2) |> ergm_LayerLogics(nw)
  coef.names <- paste0("twostarL(",
                       map_chr(Ls, toString) |> paste0(collapse=TYPEREP[type]),
                       if(a$distinct) ",distinct",
                       ")")

  iinputs <- c(typeID, a$distinct)
  list(name = "twostarL", coef.names = coef.names, iinputs = iinputs,
       auxiliaries = Ls%@%"aux", minval = 0, dependence = TRUE)
}

################################################################################

#' @templateVar name mutualL
#' @title Mutuality
#' @description In binary ERGMs, equal to the number of
#'   pairs of actors \eqn{i} and \eqn{j} for which \eqn{(i{\rightarrow}j)}{(i,j)}
#'   and \eqn{(j{\rightarrow}i)}{(j,i)} both exist.
#'
#' @details This term can only be used with directed networks.
#'
#' @usage
#' # binary: mutualL(same=NULL, diff=FALSE, by=NULL, keep=NULL, Ls=NULL)
#' @param same optional argument. If passed the name of a vertex attribute,
#'   only mutual pairs that match on the attribute are counted. Only one of `same`
#'   or `by` may be used. If both parameters are passed, `same` takes precedent. This
#'   parameter is affected by `diff`.
#' @param diff separate counts for each unique matching value can be obtained by using
#'   `diff=TRUE` with `same`.
#' @param by each node is counted separately for each mutual pair in which it
#'   occurs and the counts are tabulated by unique values of the attribute if
#'   passed the name of a vertex attribute. This means that the sum of the mutual statistics when `by` is used
#'   will equal twice the standard mutual statistic. Only one of `same`
#'   or `by` may be used. If both parameters are passed, `same` takes precedent. This
#'   parameter is not affected by `diff`.
#' @param keep a numerical vector to specify which statistics should be kept whenever the `mutual` term would
#'   ordinarily result in multiple statistics.
#' @templateVar Ls.howmany one or two
#' @templateVar Ls.interp . If given, the statistic will count the number of dyads where a tie in `Ls[[1]]` reciprocates a tie in `Ls[[2]]` and vice versa. (Note that dyad that has mutual ties in `Ls[[1]]` and in `Ls[[2]]` will add 2 to this statistic.) If a formula is given, it is replicated.
#' @template ergmTerm-Ls
#'
#' @template ergmTerm-general
#'
#' @concept directed
#' @concept frequently-used
#' @concept layer-aware
InitErgmTerm.mutualL<-function (nw, arglist, ...) {
  ## Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, directed=TRUE, bipartite=NULL,
                      varnames = c("same", "by", "diff", "keep", "Ls"),
                      vartypes = c("character", "character", "logical", "numeric", ERGM_LAYERS_SPEC),
                      defaultvalues = list(NULL, NULL, FALSE, NULL, NULL),
                      required = c(FALSE, FALSE, FALSE, FALSE, FALSE))

  assert_LHS_Layer(nw)

  ## Process the arguments
  if (!is.null(a$same) || !is.null(a$by)) {
    if (!is.null(a$same)) {
     attrname <- a$same
     if (!is.null(a$by))
       warning("Ignoring 'by' argument to mutual because 'same' exists", call.=FALSE)
    }else{
     attrname <- a$by
    }
    nodecov <- get.node.attr(nw, attrname)
    u <- NVL(a$levels, sort(unique(nodecov)))
    if (!is.null(a$keep)) {
      u <- u[a$keep]
    }
    #   Recode to numeric
    nodecov <- match(nodecov,u,nomatch=length(u)+1)
    # All of the "nomatch" should be given unique IDs so they never match:
    dontmatch <- nodecov==(length(u)+1)
    nodecov[dontmatch] <- length(u) + (1:sum(dontmatch))
    ui <- seq(along=u)
  }

  ### Construct the list to return
  if (!is.null(a$same) || !is.null(a$by)) {
    if (is.null(a$same)) {
      coef.names <- paste("mutual.by", attrname, u, sep=".")
      inputs <- c(ui, nodecov)
    }else{
     if (a$diff) {
      coef.names <- paste("mutual.same", attrname, u, sep=".")
      inputs <- c(ui, nodecov)
     }else{
      coef.names <- paste("mutual", attrname, sep=".")
      inputs <- nodecov
     }
    }
    if (is.null(a$same) && !is.null(a$by)) {
     name <- "mutual_by_attr"
    }else{
     name <- "mutual"
    }
  }else{
     name <- "mutual"
     coef.names <- "mutual"
     inputs <- NULL
  }

  maxval <- network.dyadcount(nw,FALSE)/2

  Ls <- ergm_LayerLogics(a$Ls, nw)
  L1 <- Ls[[1]]
  L2 <- Ls[[2]]
  if(!is.null(L1) || !is.null(L2)){
    NVL(L1) <- L2
    NVL(L2) <- L1
    auxiliaries <- Ls%@%"aux"
    name <- paste(name, "ML", sep="_")
    coef.names <- ergm_mk_std_op_namewrap("L", if (L1 == L2) list(L1) else list(L1, L2))(coef.names)
    maxval <- maxval*2
  }else auxiliaries <- NULL

  list(name=name,                      #name: required
       coef.names = coef.names,        #coef.names: required
       inputs=inputs,
       auxiliaries = auxiliaries,
       minval = 0,
       maxval = maxval)
}


#' @templateVar name hammingL
#' @title Hamming distance between pairs of lairs
#' @description Models marginal dependence of layers within each dyad
#'   by counting their hamming distances.
#'
#'   The term adds one statistic to the model, equalling the sum over
#'   all distinct pairs of specified layers of their Hamming distances.
#'
#' @details A positive coefficient induces negative dependence and a negative
#'   one induces positive dependence.
#'
#' @usage
#' # binary: hammingL(Ls=~.)
#'
#' @templateVar Ls.howmany at least two
#' @templateVar Ls.interp .
#' @template ergmTerm-Ls
#'
#' @template ergmTerm-general
#'
#' @concept directed
#' @concept undirected
#' @concept layer-aware
InitErgmTerm.hammingL <- function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("Ls"),
                      vartypes = c(ERGM_LAYERS_SPEC),
                      defaultvalues = list(empty_env(~.)),
                      required = c(FALSE))

  assert_LHS_Layer(nw)

  Ls <- ergm_LayerLogics(a$Ls, nw)
  if (length(Ls) < 2L) ergm_Init_stop("multiple layers are required")
  deps <- lapply(Ls, attr, "dep")

  affects <- map(seq_len(network.layercount(nw)),
                 function(l) which(map_lgl(deps, function(d) l %in% d)))

  iinputs <- c(0L, cumsum(c(0L, lengths(affects))) + length(affects) + 1L, unlist(affects) - 1L)

  list(name="pairwisedistL", coef.names = paste0('hammingL(', toString(Ls), ')'),
       iinputs = iinputs, dependence = TRUE, auxiliaries = Ls%@%"aux")
}


#' @templateVar name entrainmentL
#' @title Entrainment between pairs of lairs
#' @description Models marginal dependence of layers within each dyad
#'   by counting edges they have in common.
#'
#'   The term adds one statistic to the model, equalling the sum over
#'   all distinct pairs of specified layers of their common edges.
#'
#' @details This is similar to `L(~edges, c(~A & B, ~A & C, ~B & C))`
#'   but more efficient for multiple layers.
#'
#' @usage
#' # binary: entrainmentL(Ls=~.)
#'
#' @templateVar Ls.howmany at least two
#' @templateVar Ls.interp .
#' @template ergmTerm-Ls
#'
#' @template ergmTerm-general
#'
#' @concept directed
#' @concept undirected
#' @concept layer-aware
InitErgmTerm.entrainmentL <- function(...){
  # Rename the function to avoid the extra nesting level in the
  # diagnostic messages.
  f <- InitErgmTerm.hammingL
  term <- f(...)
  term$coef.names <- gsub("^hammingL", "entrainmentL", term$coef.names)
  term$iinputs[1] <- 1L
  term
}
