#  File R/InitErgmTerm.multilayer.R in package ergm.multi, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2022 Statnet Commons
################################################################################
## TODO: LL-Constrained proposals.
## TODO: Check that noncommutative LL operators work as intended.

.same_constraints <- function(nwl, nattr){
  map(nwl, get.network.attribute, nattr) %>% map(NVL, ~.) %>% map(empty_env) %>% all_identical
}

.layer_namemap <- function(nw){
  nwnames <- .peek_vattrv(nw, ".LayerName", missing="NULL")
  if(is.numeric(nwnames)) nwnames <- NULL
  nwids <- .peek_vattrv(nw, ".LayerID")

  o <- structure(nwids, names=nwnames)
  o[!duplicated(o)]
}

.set_layer_namemap <- function(ll, nw){
  namemap <- if(is.network(nw)) .layer_namemap(nw) else nw
  if(is(ll,"formula")) ergm_LayerLogic(ll, namemap)
  else lapply(ll, ergm_LayerLogic, namemap)
}

.lspec_coef.namewrap <- function(Llist, collapse=TRUE){
  reprs <- sapply(seq_along(Llist), function(l){
    name <- names(Llist)[l]
    L <- Llist[[l]]
    s <- NVL3(L, toString(ergm_LayerLogic(.)), "")
    if(NVL(name,"")!="") s <- paste0(name,"=",s)
    s
  })
  if(collapse) ergm_mk_std_op_namewrap("L", reprs) else reprs
}

#' Construct a "view" of a network.
#'
#' Returns a network with edges optionally filtered according to a
#' specified criterion and with edge attributes optionally computed
#' from other edge attributes.
#'
#' @param x a [`network`] object.
#' @param ... a list of attribute or filtering specifications. See
#'   Details.
#' @param .clear whether the edge attributes not set by this call
#'   should be deleted.
#' @param .sep when specifying via a character vector, use this as the
#'   separator for concatenating edge values.
#'
#' @details Attribute specification arguments have the form
#'   `<newattrname> = <expr>`, where `<newattrname>` specifies the
#'   name of the new edge attribute (or attribute to be overwritten)
#'   and `<expr>` can be one of the following:
#' \describe{
#'
#' \item{a function}{The function will be passed two arguments, the
#' edgelist [`tibble`] and the network, and must return a vector of
#' edge attribute values to be set on the edges in the order
#' specified.}
#'
#' \item{a formula}{The expression on the RHS of the formula will be
#' evaluated with names in it referencing the edge attributes. The
#' input network may be referenced as `.nw`. The expression's result
#' is expected to be a vector of edge attribute values to be set on
#' the edges in the order specified.}
#' 
#' \item{a character vector}{If of length one, the edge attribute with
#' that name will simply be copied; if greater than one, the attribute
#' values will be concatenated wtih the `.sep` argument as the
#' separator.}
#'
#' \item{an object enclosed in [I()]}{The object will be used directly
#' to set the edge attribute.}
#' }
#'
#' Filtering arguments are specified the same way as attribute
#' arguments, but they must be unnamed arguments (i.e., must be passed
#' without the `=`) and must return a logical or numeric vector
#' suitable for indexing the edge list. Multiple filtering arguments
#' may be specified, and the edge will be kept if it satisfies
#' *all*. If the conjunction of the edge's original states and the
#' filtering results is ambiguous (i.e., `NA`), it will be set as
#' missing.
#'
#' @return A [`network`] object with modified edges and edge attributes.
#'
#' @examples
#' data(florentine)
#' flo <- flomarriage
#' flo[,,add.edges=TRUE] <- as.matrix(flomarriage) | as.matrix(flobusiness)
#' flo[,, names.eval="m"] <- as.matrix(flomarriage)==1
#' flobusiness[3,5] <- NA
#' flo[,, names.eval="b"] <- as.matrix(flobusiness)==1
#' flo
#' (flob <- network_view(flo, "b"))
#' (flobusiness) # for comparison
#' \dontshow{
#' if(require(testthat, quietly=TRUE))
#' testthat::expect_equivalent(as.matrix(flobusiness),as.matrix(flob))
#' }
#' 
#' (flob <- network_view(flo, ~b&m))
#' (flobusiness & flomarriage) # for comparison
#' \dontshow{
#' if(require(testthat, quietly=TRUE))
#' testthat::expect_equivalent(as.matrix(flobusiness & flomarriage),as.matrix(flob))
#' }
#'
#' as.matrix(flob <- network_view(flo, bm=~b+m), attrname="bm")
#' (as.matrix(flobusiness) + as.matrix(flomarriage)) # for comparison
#' \dontshow{
#' if(require(testthat, quietly=TRUE))
#' testthat::expect_equivalent(as.matrix(flobusiness)+as.matrix(flomarriage),as.matrix(flob, attrname="bm"))
#' }
#'
#' as.matrix(flob <- network_view(flo, ~b, bm=~b+m), attrname="bm")
#' as.matrix(flobusiness)*(1+as.matrix(flomarriage)) # for comparison
#' \dontshow{
#' if(require(testthat, quietly=TRUE))
#' testthat::expect_equivalent(as.matrix(flobusiness)*(1+as.matrix(flomarriage)),as.matrix(flob, attrname="bm"))
#' }
#' 
#' 
#' @export
network_view <- function(x, ..., .clear=FALSE, .sep="."){
  # Handle empty network
  if(network.edgecount(x,na.omit=FALSE)==0) return(x)
  
  exprs <- list(...)
  fes <- exprs[NVL(names(exprs),rep("",length(exprs)))==""]
  oes <- exprs[NVL(names(exprs),rep("",length(exprs)))!=""]

  #' @importFrom rlang abort
  evl <- function(e, el, x){
    switch(class(e),
           `function` = e(.el=el, .nw=x),
           formula = eval(e[[length(e)]], envir=c(list(.nw=x), as.list(el)), enclos=environment(e)),
           character = if(length(e)==1) el[[e]] else do.call(paste, c(as.list(el[e]), sep=.sep)),
           AsIs = e,
           abort("Unsupported specification for network_view."))
  }

  if(length(fes)){
    el <- as_tibble(x, attrnames=list.edge.attributes(x), na.rm=FALSE, store.eid=TRUE)
    keep <- !el$na
    for(e in fes){
      keep <- keep & evl(e, el, x)
    }
    del <- na.omit(el$.eid[!keep])
    nael <- el$.eid[is.na(keep)]
    if(length(del)) delete.edges(x, del)
    if(length(nael)) set.edge.attribute(x, "na", TRUE, nael)
    # This one applies to a rare situation where the edge is missing
    # but the filter says that it should be present.
    add <- el$.eid[el$na & !is.na(keep) & keep]
    if(length(add)) set.edge.attribute(x, "na", FALSE, add)
  }

  for(i in seq_along(oes)){
    el <- as_tibble(x, attrname=list.edge.attributes(x), na.rm=FALSE, store.eid=TRUE)
    
    e <- oes[[i]]
    nm <- names(oes)[[i]]
    
    newval <- evl(e, el, x)
    set.edge.attribute(x, nm, newval, el$.eid)
  }
  
  if(.clear) for(a in setdiff(list.edge.attributes(x), c(nm,"na"))) delete.edge.attribute(x, a)
  x
}


#' Returns a directed version of an undirected binary network
#'
#' @param x a [`network`] object.
#' @param rule a string specifying how the network is to be
#'   constructed.
direct.network <- function(x, rule=c("both", "upper", "lower")){
  rule <- match.arg(rule)

  el <- as.edgelist(x)
  el <- switch(rule,
               both = rbind(el, el[,2:1,drop=FALSE]),
               upper = cbind(pmin(el[,1],el[,2]),pmax(el[,1],el[,2])),
               lower = cbind(pmax(el[,1],el[,2]),pmin(el[,1],el[,2])))
  
  o <- network.initialize(network.size(x), directed=TRUE, bipartite=x%n%"bipartite", loops=has.loops(x), hyper=is.hyper(x), multiple=is.multiplex(x))
  o <- network.edgelist(el, o)
  nvattr.copy.network(o, x)
}

#' A multilayer network representation.
#'
#' A function for specifying the LHS of a multilayer
#' (a.k.a. multiplex, a.k.a. multirelational, a.k.a. multivariate)
#' ERGM in the framework of \insertCite{KrKo20e;textual}{ergm.multi}.
#'
#' @param ... layer specification, in one of three formats:
#' 
#'   1. An (optionally named) list of identically-dimensioned
#'      networks.
#'
#'   1. Several networks as (optionally named) arguments.
#'
#'   1. A single network, a character vector, and several optional
#'      arguments. Then, the layers are values of the named edge
#'      attributes. If the vector has named elements (e.g.,
#'      `c(a="advice", c="collaboration")`), the layers will be
#'      renamed accordingly. The optional arguments `.symmetric` and
#'      `.bipartite` are then interpreted as described below.
#'
#' @param .symmetric If the layer specification is via a single
#'   network with edge attributes and the network is directed, an
#'   optional logical vector to specify which of the layers should be
#'   treated as undirected.
#' @param .bipartite If the layer specification is via a single
#'   network with edge attributes and the network is unipartite, an
#'   optional integer vector to specify which of the layers should be
#'   treated as bipartite and how many `b1` vertices there are.
#' @param .active A [nodal attribute specification][nodal_attributes]
#'   (`? nodal_attributes`) specifying which nodes on each network
#'   *may* have ties, or a list with an element for each network. The
#'   list will be recycled up to the number of layers.
#'
#' @return A network object with layer metadata.
#'
#' @note The resulting network will be the "least common denominator"
#'   network: if not all layers have the same bipartedness, all layers
#'   will appear as unipartite to the statistics, and if any are
#'   directed, all will be. However, [certain operator
#'   terms][ergmTerm], particularly `Symmetrize()` and `S()`, can be
#'   used to construct a bipartite subgraph of a unipartite graph or
#'   change directedness.
#'
#' @section Specifying models for multilayer network:
#' In order to fit a model for multilayer
#' networks, first use [`Layer`] construct an LHS network that
#' [ergm()] will understand as multilayered.
#'
#' Used in the formula directly, most, but not all, \pkg{ergm} terms will
#' sum their statistics over the observed layers.
#'
#' Some terms are *layer-aware*, however. By convention, layer-aware
#' terms have capital `L` appended to them. For example,
#' [`mutualL`][mutualL-ergmTerm] is a layer-aware generalization of
#' [`mutual`][mutual-ergmTerm]. These terms have one or more explicit
#' (usually optional) layer specification arguments. By convention, an
#' argument that requires one layer specification is named \code{L=} and
#' one that requires a list of specifications (constructed by [list()]
#' or [c()]) is named \code{Ls=}; and a specification of the form `~.` is a
#' placeholder for all observed layers.
#'
#' Operator [`L(formula, Ls=...)`][L-ergmTerm] can be used to evaluate
#' arbitrary terms in the `formula` on specified layers.
#'
#' Layer specification documentation follows.
#'
#' ## Layer Logic
#'
#' Each formula's right-hand side describes an observed layer *or* some
#' "logical" layer, whose ties are a function of corresponding ties in
#' observed layers. \insertCite{KrKo20e}{ergm.multi}
#'
#' The observed layers can be referenced either by name or by number (i.e.,
#' order in which they were passed to [`Layer`]). When referencing by
#' number, enclose the number in quotation marks (e.g., "1") or
#' backticks (e.g., \dQuote{`1`}).
#'
#' \link[base:Arithmetic]{Arithmetical}, \link[base:Comparison]{relational},
#' and \link[base:Logic]{logical} operators can be used to combine them. All
#' listed operators are implemented, as well as functions [`abs`],
#' [`round`], and [`sign`]. Standard
#' \link[base:Syntax]{operator precedence} applies, so use of parentheses is
#' recommended to ensure the logical expression is what it looks like.
#'
#' \strong{Important:} For performance reasons, \pkg{ergm.multi}'s
#' Layer Logic implementation uses integer arithmetic. This means, in
#' particular, that `/` will round down instead of returning a
#' fraction (as `%/%` does in \R), and [round()] function without a
#' second argument (which can be negative to round to the nearest 10,
#' 100, etc.) is not meaningful and will be ignored.
#'
#' For example, if LHS is \code{Layer(A=nwA, B=nwB)}, both \code{~`2`} and
#' \code{~B} refer to \code{nwB}, while \code{A&!B} refers to a
#' \dQuote{logical} layer that has ties that are in \code{nwA} but not in
#' \code{nwB}.
#'
#' Transpose function \code{\link{t}} applied to a directed layer will reverse
#' the direction of all relations (transposing the sociomatrix). Unlike the
#' others, it can only be used on an observed layer directly. For example,
#' \code{~t(`1`)&t(`2`)} is valid but \code{~t(`1`&`2`)} is not.
#'
#' At this time, logical expressions that produce complete graphs from empty
#' graph inputs (e.g., \code{A==B} or \code{!A}) are not supported.
#'
#' ## Summing layers
#'
#' Some of the terms that call for a list of layers (i.e., have \code{Ls=}
#' arguments) will sum the statistic over the layers. For example,
#' \code{Layer(nw1,nw2)~L(~edges, c(~`1`,~(`2`&!`1`)))} produces the
#' number of edges in layer 1 plus the number of edges in layer 2 but
#' not in layer 1.
#'
#' For these formulas, one can specify the layer's weight on its left-handside.
#' For example, \code{Layer(nw1,nw2)~L(~edges, c(3~`1`,-1~(`2`&!`1`)))} will
#' produce three times the number of edges in layer 1, minus the number of
#' edges in layer 2 but not in layer 1.
#'
#' @seealso [Help on model specification][ergmTerm] for specific terms.
#'
#' @references \insertAllCited{}
#' 
#' @examples
#'
#' data(florentine)
#'
#' # Method 1: list of networks
#' flo <- Layer(list(m = flomarriage, b = flobusiness))
#' ergm(flo ~ L(~edges, ~m)+L(~edges, ~b))
#'
#' # Method 2: networks as arguments
#' flo <- Layer(m = flomarriage, b = flobusiness)
#' ergm(flo ~ L(~edges, ~m)+L(~edges, ~b))
#'
#' # Method 3: edge attributes (also illustrating renaming):
#' flo <- flomarriage | flobusiness
#' flo[,, names.eval="marriage"] <- as.matrix(flomarriage)
#' flo[,, names.eval="business"] <- as.matrix(flobusiness)
#' flo # edge attributes
#' flo <- Layer(flo, c(m="marriage", b="business"))
#' ergm(flo ~ L(~edges, ~m)+L(~edges, ~b))
#'
#' ### Specifying modes and mixed bipartitedness
#'
#' # Suppose we have a two-mode network with 5 nodes on Mode 1 and 15
#' # on Mode 2, and suppose that we observe two layers, one only among
#' # actors of Mode 1 and the other bipartite between Modes 1 and 2.
#'
#' # Construct the two layers' networks:
#' nw1 <- network.initialize(20, dir=FALSE)
#' nw12 <- network.initialize(20, dir=FALSE, bipartite=5)
#' nw1 %v% "mode" <- rep(1:2,c(5,15))
#'
#' # For testing: the maximal set of edges for each type of network:
#' nw1[1:5,1:5] <- 1
#' nw12[1:5,6:20] <- 1
#'
#' # The .active argument specifies the following:
#' # * nw1's vertices are only active if their mode=1 (i.e., 1-2, 2-1,
#' #   and 2-2 can't have edges).
#' # * nw12's vertices are all active, but the network is bipartite,
#' #   so constraints will be adjusted automatically.
#' lnw <- Layer(nw1, nw12, .active=list(~mode==1, ~TRUE))
#'
#' summary(lnw~
#' edges+ # 5*4/2+5*15 = 10+75 = 85
#' L(~edges,~`1`)+ # 5*4/2 = 10
#' L(~edges,~`2`)+ # 5*15 = 75
#' L(~edges,~(`1`|`2`))+ # This logical layer has contents of both, so also 85.
#' L(~edges,~(`1`&`2`)) # There is no overlap between the two layers, so 0.
#' )
#'
#' # Layer-aware terms can be used:
#'
#' nw1[,] <-0
#' nw1[1,2:3] <- 1
#' nw1[2,3] <- 1
#' nw12[,] <- 0
#' nw12[1,6:7] <- 1
#' nw12[2,6:7] <- 1
#'
#' lnw <- Layer(nw1, nw12, .active=list(~mode==1,~TRUE))
#'
#' summary(lnw~L(~triangles, ~`1`)+ # 1-2-3 triangle.
#'   L(~triangles, ~`1`|`2`)+ # 1-2-3, 1-2-6, 1-2-7 triangles
#'   dgwespL(L.base=~`1`, Ls.path=list(~`2`,~`2`)) # 1-2-6 and 1-2-7 only
#' )
#'
#' # Because the layers are represented as a block-diagonal matrix,
#' # this will only count triangles entirely contained within a single
#' # layer, i.e., 1-2-3:
#' summary(lnw~triangles)
#'
#' # If you need to evaluate bipartite-only statistics on the second
#' # layer, you need to use the S() operator to select the bipartite
#' # view:
#' summary(lnw~L(~S(~b1degree(1:3)+b2degree(1:3),1:5~6:20), ~`2`))
#'
#' @export
Layer <- function(..., .symmetric=NULL, .bipartite=NULL, .active=NULL){
  args <- list(...)
  if(all(sapply(args, is, "network"))){
    nwl <- args
  }else if(is.list(args[[1]]) && all(sapply(args[[1]], is, "network"))){
    nwl <- args[[1]]
  }else if(length(args)>=1 && is(args[[1]], "network") && is(args[[2]], "vector")){
    nw <- args[[1]]
    nwl <-
      lapply(args[[2]], function(eattr){
        network_view(nw, eattr)
      })

    if(is.directed(nw) && !is.null(.symmetric)){
      symm <- as.logical(.symmetric)
      for(i in which(symm)){
        # There is probably a more efficient way to do this, but we need
        # to compute nw1 anyway.
        nw1 <- ergm_symmetrize(nwl[[i]], rule="weak")
        nw2 <- ergm_symmetrize(nwl[[i]], rule="strong")
        if(!identical(as.vector(as.edgelist(nw1)),as.vector(as.edgelist(nw2))))
          stop("Layer specified to be treated as undirected is not symmetric.")
        nwl[[i]] <- nw1
      }
    }

    if(!is.bipartite(nw) && !is.null(.bipartite)){
      bip <- as.integer(.bipartite)
      for(i in which(bip!=0)){
        nw1 <- nwl[[i]]
        nw1 %n% "bipartite" <- bip[i]
        nwl[[i]] <- nw1
      }
    }

    # Set names: if args[[2]] is a named vector, use those where set.
    names(nwl) <- NVL3(names(args[[2]]), ifelse(.=="", args[[2]], .), args[[2]])

  }else stop("Unrecognized format for multilayer specification. See help for information.")

  if(!is.null(.active)){
    if(!is.list(.active)) .active <- list(.active)
    .active <- rep(.active, length.out=length(nwl))
    al <- Map(function(nw, a) ergm_get_vattr(a, nw, accept="logical"),
              nwl, .active)
    if(!all(unlist(al)))
      nwl <- Map(function(nw, a) blacklist_intersect(nw, a, invert=TRUE),
                 nwl, al)
  }

  # nwl may now be a list with networks of heterogeneous bipartitedness.
  bip <- sapply(nwl, `%n%`, "bipartite") %>% sapply(NVL, 0L)
  blockout <- if(all_identical(bip)) rep(FALSE, length(nwl)) else bip

  nwl <- Map(function(nw, b){
    nw %v% ".bipartite" <- b;
    if(b){
      n <- network.size(nw)
      v <- rep(c(TRUE,FALSE), c(b,n-b))
      # Blacklist b1-b1 ties and b2-b2 ties for the purpose of
      # sampling.
      nw <- nw %>% blacklist_intersect(v) %>% blacklist_intersect(!v)
      # Bipartiteness is now enforced by the blacklist.
      nw %n% "bipartite" <- FALSE
    }
    nw
  }, nwl, blockout)

  # nwl may now be a list with networks of heterogeneous directedness.
  
  dir <- sapply(nwl, is.directed)
  symm <- if(all_identical(dir)) rep(FALSE, length(nwl)) else !dir

  nwl <- Map(function(nw, symm) {
    nw %v% ".undirected" <- symm
    if(symm) direct.network(nw,rule="upper") else nw
  }, nwl, symm)

  # nwl is now a list of networks with homogeneous directedness, some
  # networks tagged with vertex attribute .undirected.

  # Perform some checks and imputations for layer names.
  nnames <- names(nwl)
  if(!is.null(nnames) && any(blank<-(nnames==""))){
    warning("Only some of the layers have specified names; they have been imputed with the corresponding layer number.")
    nnames[blank] <- as.character(seq_along(nnames)[blank])
  }
  if(any(
  regexpr('^[0-9]+$',nnames)!=-1 # Names that are integers are potentially problematic,
  & nnames!=seq_along(nnames) # but not if they happen to match layer IDs.
  )) warning("Using numeric layer names is ambiguous.")
  if(anyDuplicated(nnames)) stop("Duplicate layer names.")
  names(nwl) <- nnames

  if(!.same_constraints(nwl, "constraints")) stop("Layers have differing constraint structures. This is not supported at this time.")
  if(!.same_constraints(nwl, "obs.constraints")) stop("Layers have differing observation processes. This is not supported at this time.")
  
  nw <- combine_networks(nwl, blockID.vattr=".LayerID", blockName.vattr=".LayerName", ignore.nattr = c(eval(formals(combine_networks)$ignore.nattr), "constraints", "obs.constraints", "ergm"), subnet.cache=TRUE)

  nw %n% "ergm" <- combine_ergmlhs(nwl)

  nw %ergmlhs% "constraints" <-
      if(NVL(nwl[[1]] %ergmlhs% "constraints",base_env(~.))==base_env(~.))
        base_env(~blockdiag(".LayerID"))
      else
        append_rhs.formula(nwl[[1]] %ergmlhs% "constraints", list(call("blockdiag",".LayerID")), TRUE)
  
  if(any(symm)) nw %ergmlhs% "constraints" <- append_rhs.formula(nw %ergmlhs% "constraints", list(call("upper_tri",".undirected")), TRUE)
  if(any(blockout!=0)||!is.null(.active)) nw %ergmlhs% "constraints" <- append_rhs.formula(nw %ergmlhs% "constraints", list(call("blacklist_block")), TRUE)

  if(!is.null(nwl[[1]]%ergmlhs%"obs.constraints")) nw %ergmlhs% "obs.constraints" <- nwl[[1]] %ergmlhs% "obs.constraints"

  nw
}


## InitErgmTerm..layer.nets <- function(nw, arglist, ...){
##   a <- check.ErgmTerm(nw, arglist,
##                       varnames = c(),
##                       vartypes = c(),
##                       defaultvalues = list(),
##                       required = c())
##   list(name="_layer_nets", coef.names=c(), inputs=unlist(.block_vertexmap(nw, ".LayerID", TRUE)), dependence=FALSE)
## }

.depends_on_layers <- function(commands){
  coms <- commands[-1]
  if(any(coms==0)) coms <- coms[-(which(coms==0)+1)] # Drop all numeric literals (i.e., numbers preceded by 0).
  coms <- coms[coms>=1] # Drop all commands.
  unique(coms)
}

InitErgmTerm..layer.net <- function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("L"),
                      vartypes = c("formula"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))

  
  nwl <- subnetwork_templates(nw,".LayerID",".LayerName")

  ll <- to_ergm_Cdouble(.set_layer_namemap(a$L, nw))
  # Terms on this logical layer will induce dyadic independence if its
  # value depends on more than one other layer value.
  dependence <- length(.depends_on_layers(ll))>1
  
  if(test_eval.LayerLogic(ll, FALSE)) stop("Layer specifications that produce edges on the output layer for empty input layers are not supported at this time.", call.=FALSE)
  
  list(name="_layer_net", coef.names=c(), iinputs=c(unlist(.block_vertexmap(nw, ".LayerID", TRUE)), if(is.directed(nw)) sapply(nwl, function(nw) (nw%v% ".undirected")[1]), ll), dependence=dependence)
}

LL_PREOPMAP <- list(
    # Unary operators
    c(`t` = -21)
  )
LL_POSTOPMAP <- list(
    # Unary operators
    c(`(` = NA,
      `!` = -1,
      `+` = NA,
      `-` = -16,
      `abs` = -17,
      `round` = NA,
      `sign` = -20),
    # Binary operators
    c(`&` = -2,
      `&&` = -2,
      `|` = -3,
      `||` = -3,
      `xor` = -4,
      `==` = -5,
      `!=` = -6,
      `<` = -7,
      `>` = -8,
      `<=` = -9,
      `>=` = -10,
      `+` = -11,
      `-` = -12,
      `*` = -13,
      `/` = -14,
      `%%` = -15,
      `^` = -18,
      `%/%` = -14,
      `round` = -19)
)

LL_IDEMPOTENT <- c("&", "&&", "|", "||")
LL_TAUTOLOGICAL <- c("==", ">=", "<=")
LL_CONTRADICTORY <- c("!=", "xor", "<", ">")

#' Internal representation of Layer Logic
#'
#' @param formula A Layer Logic formula.
#' @param namemap A character vector giving the names of the layers
#'   referenced, or `NULL`.
#'
#' @return A structure with nonce class
#'   `c("ergm_LayerLogic",class(formula))`, comprising the input
#'   `formula` and an attribute `namemap` containing the `namemap`.
#' @keywords internal
#' @export
ergm_LayerLogic <- function(formula, namemap=NULL){
  ## TODO: Check whether we should verify that this is a formula.
  structure(formula, namemap=namemap, class=if(is(formula,"ergm_LayerLogic")) class(formula) else c("ergm_LayerLogic", class(formula)))
}

#' @describeIn ergm_LayerLogic A method to generate coefficient names
#'   associated with the Layer Logic.
#' @param x An `ergm_LayerLogic` object.
#' @param ... Additional arguments, currently unused.
#' @export
toString.ergm_LayerLogic <- function(x, ...){
  fmt <- function(x){
    class(x) <- keep(class(x), `!=`, "ergm_LayerLogic")
    switch(class(x),
           formula = despace(deparse(if(length(x)==2) x[[2]] else x)),
           character = x,
           list = paste0('(',paste(sapply(x,fmt),collapse=","),')'),
           as.character(x))
  }
  fmt(x)
}

# Substitute numeric layer IDs for layer names and simplify some common trivial expressions (e.g., a | a => a).
sub.ergm_LayerLogic <- function(x){
  formula <- x
  namemap <- attr(formula, "namemap")

  lidSub <- function(l){
    switch(class(l),
           logical =,
           numeric = l,
           character =,
           name = {
             l <- as.character(l)
             if(regexpr('^[0-9]+$',l)!=-1) as.name(as.integer(l))
             else if(l %in% names(namemap)) as.name(namemap[l])
             else stop("Unrecognised symbol ",sQuote(l)," in layer logic.", call.=FALSE)
           }
           )
  }

  simplify <- function(call){
    if(is.call(call)){
      op <- call[[1]]

      call[-1] <- lapply(call[-1], simplify)

      if(as.character(op) %in% LL_IDEMPOTENT && all_identical(as.list(call)[-1])) call[[2]]
      else if(as.character(op) %in% LL_TAUTOLOGICAL && all_identical(as.list(call)[-1])) TRUE
      else if(as.character(op) %in% LL_CONTRADICTORY && length(call)%%2==0 && all_identical(as.list(call)[-1])) FALSE
      else call
    }else{
      lidSub(call)
    }
  }

  replace(formula, length(formula), list(simplify(ult(formula))))
}

#' @describeIn ergm_LayerLogic A method to encode and serialize the
#'   Layer Logic into a postfix program understood by the C code.
#' @export
to_ergm_Cdouble.ergm_LayerLogic <- function(x, ...){
  formula <- sub.ergm_LayerLogic(x)

  lidMap <- function(l){
    switch(class(l),
           logical =,
           numeric = c(0,l),
           character =,
           name = if(regexpr('^[0-9]+$',l)!=-1) as.integer(as.character(l))
                  else stop("Unrecognised symbol ",sQuote(l)," in layer logic.", call.=FALSE))
  }

  postfix <- function(call, coml=c()){
    if(is.call(call)){
      op <- call[[1]]
      if(as.character(op) %in% unlist(lapply(LL_PREOPMAP, names))){
        coml <- c(coml, LL_PREOPMAP[[length(call)-1]][[as.character(op)]])
        postop <- FALSE
      }else postop <- TRUE
      for(i in seq_along(call[-1])+1){
        coml <- c(coml, postfix(call[[i]]))
      }
      if(postop) coml <- c(coml, LL_POSTOPMAP[[length(call)-1]][as.character(op)])
    }else{
      coml <- c(coml, lidMap(call))
    }
    coml[!is.na(coml)]
  }
  
  com <- postfix(ult(formula))
  c(sum(com!=0 & !com%in%unlist(LL_PREOPMAP)), com)
}

test_eval.LayerLogic <- function(commands, lv, lvr = lv){
  coms <- as.integer(commands[-1])
  lv <- rep(lv, length.out=max(coms))
  stack <- integer(0)
  if(sum(coms!=0 & !coms%in%unlist(LL_PREOPMAP))!=commands[1]) stop("Layer specification command vector specifies incorrect number of commands.", call.=FALSE)
  for(i in 1:commands[1]){
    com <- coms[1]
    if(com==0){
      coms <- coms[-1]
      com <- coms[1]
      stack <- c(com, stack)
    }else if(com==-1){
      x0 <- stack[1]; stack <- stack[-1]
      stack <- c(!x0, stack)
    }else if(com==-2){
      x0 <- stack[1]; stack <- stack[-1]
      y0 <- stack[1]; stack <- stack[-1]
      stack <- c(x0 && y0, stack)
    }else if(com==-3){
      x0 <- stack[1]; stack <- stack[-1]
      y0 <- stack[1]; stack <- stack[-1]
      stack <- c(x0 || y0, stack)
    }else if(com==-4){
      x0 <- stack[1]; stack <- stack[-1]
      y0 <- stack[1]; stack <- stack[-1]
      stack <- c(xor(x0, y0), stack)
    }else if(com==-5){
      x0 <- stack[1]; stack <- stack[-1]
      y0 <- stack[1]; stack <- stack[-1]
      stack <- c(x0 == y0, stack)
    }else if(com==-6){
      x0 <- stack[1]; stack <- stack[-1]
      y0 <- stack[1]; stack <- stack[-1]
      stack <- c(x0 != y0, stack)
    }else if(com==-7){
      x0 <- stack[1]; stack <- stack[-1]
      y0 <- stack[1]; stack <- stack[-1]
      stack <- c(x0 < y0, stack)
    }else if(com==-8){
      x0 <- stack[1]; stack <- stack[-1]
      y0 <- stack[1]; stack <- stack[-1]
      stack <- c(x0 > y0, stack)
    }else if(com==-9){
      x0 <- stack[1]; stack <- stack[-1]
      y0 <- stack[1]; stack <- stack[-1]
      stack <- c(x0 <= y0, stack)
    }else if(com==-10){
      x0 <- stack[1]; stack <- stack[-1]
      y0 <- stack[1]; stack <- stack[-1]
      stack <- c(x0 >= y0, stack)
    }else if(com==-11){
      x0 <- stack[1]; stack <- stack[-1]
      y0 <- stack[1]; stack <- stack[-1]
      stack <- c(x0 + y0, stack)
    }else if(com==-12){
      x0 <- stack[1]; stack <- stack[-1]
      y0 <- stack[1]; stack <- stack[-1]
      stack <- c(x0 - y0, stack)
    }else if(com==-13){
      x0 <- stack[1]; stack <- stack[-1]
      y0 <- stack[1]; stack <- stack[-1]
      stack <- c(x0 * y0, stack)
    }else if(com==-14){
      x0 <- stack[1]; stack <- stack[-1]
      y0 <- stack[1]; stack <- stack[-1]
      stack <- c(x0 %/% y0, stack)
    }else if(com==-15){
      x0 <- stack[1]; stack <- stack[-1]
      y0 <- stack[1]; stack <- stack[-1]
      stack <- c(x0 %% y0, stack)
    }else if(com==-16){
      x0 <- stack[1]; stack <- stack[-1]
      stack <- c(-x0, stack)
    }else if(com==-17){
      x0 <- stack[1]; stack <- stack[-1]
      stack <- c(abs(x0), stack)
    }else if(com==-18){
      x0 <- stack[1]; stack <- stack[-1]
      y0 <- stack[1]; stack <- stack[-1]
      stack <- c(as.integer(x0 ^ y0), stack)
    }else if(com==-19){
      x0 <- stack[1]; stack <- stack[-1]
      y0 <- stack[1]; stack <- stack[-1]
      stack <- c(as.integer(round(x0, y0)), stack)
    }else if(com==-20){
      x0 <- stack[1]; stack <- stack[-1]
      stack <- c(sign(x0), stack)
    }else if(com==-21){ 
      coms <- coms[-1]
      x0 <- coms[1]
      stack <- c(lvr[x0], stack)
    }else{
      stack <- c(lv[com], stack)
    }
    coms <- coms[-1]
  }
  if(length(stack)!=1) stop("Invalid layer specification command sequence.", call.=FALSE)
  stack
}

.all_layers_terms <- function(n, LHS=NULL){
  if(is.null(LHS))
    lapply(seq_len(n), function(i) empty_env(as.formula(substitute(~x,list(x=as.name(i))))))
  else
    lapply(seq_len(n), function(i) empty_env(as.formula(substitute(lhs~x,list(lhs=LHS, x=as.name(i))))))
}

.mk_.layer.net_auxform <- function(ll){
  trmcalls <- .layers_expand_dot(ll) %>%
    lapply(function(ltrm){
      if(length(ltrm)==3) ltrm[[2]] <- NULL
      ltrm
    }) %>%
    lapply(sub.ergm_LayerLogic)
  # Get the formula as a list of term calls.o
  trmcalls <- lapply(trmcalls, function(ltrm) call(".layer.net", ltrm))
  # TODO: Double-check that base_env here doesn't break anything.
  base_env(append_rhs.formula(~.,trmcalls)[-2])
}

.layers_expand_dot <- function(ll){
  if(is(ll, "formula")) ll <- list(ll)
  nl <- length(attr(ll[[1]], "namemap"))
  # Replace . with all layers.
  do.call(c, lapply(ll, function(f) if(f[[length(f)]]=='.') .all_layers_terms(nl, LHS = if(length(f)==3) f[[2]]) else list(as.formula(f))))
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
#' @template ergmTerm-general
#'
#' @concept operator
#' @concept layer-aware
InitErgmTerm.L <- function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula", "Ls"),
                      vartypes = c("formula", "formula,list"),
                      defaultvalues = list(NULL, empty_env(~.)),
                      required = c(TRUE, FALSE))

  nwl <- subnetwork_templates(nw,".LayerID",".LayerName")

  Ls <- .set_layer_namemap(a$Ls, nw)
  if(is(Ls, "formula")) Ls <- list(Ls)
  Ls.dotexp <- .layers_expand_dot(Ls)

  auxiliaries <- .mk_.layer.net_auxform(Ls)
  nltrms <- length(list_rhs.formula(auxiliaries))

  w <- rep(1,nltrms)
  have.LHS <- sapply(Ls.dotexp, length)==3
  w[have.LHS] <- as.numeric(sapply(lapply(Ls.dotexp[have.LHS], "[[", 2), eval,environment(Ls[[1]])))
  
  nw1 <- nwl[[1]]
  m <- ergm_model(a$formula, nw1, ..., offset.decorate=FALSE)

  ## FIXME: Is this consistent with extended state API, or do we need to have a different "model" for each layer?
  wm <- wrap.ergm_model(m, nw1, function(x) .lspec_coef.namewrap(list(a$Ls))(x))
  gs <- wm$emptynwstats
  wm$emptynwstats <- if(!is.null(gs)) gs*nltrms
  wm$dependence <- wm$dependence || !is.dyad.independent(auxiliaries, basis=nw)

  c(list(name="OnLayer", iinputs=nltrms, inputs=w, submodel=m, auxiliaries = auxiliaries),
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
                      vartypes = c("formula,list"),
                      defaultvalues = list(empty_env(~.)),
                      required = c(FALSE))

  Ls <- .set_layer_namemap(a$Ls, nw)
  auxiliaries <- .mk_.layer.net_auxform(Ls)
  nltrms <- length(list_rhs.formula(auxiliaries))

  list(name="layerCMB", coef.names = paste0('CMBL(',despace(deparse(Ls)),')'), iinputs=nltrms, dependence=TRUE, auxiliaries = auxiliaries)
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
#' @param type determines what is counted:
#'   1) *"any"* Number of configurations
#'   \eqn{(i-j), (i-k)}{(i,j), (i,k)} , where
#'   \eqn{(i-j)}{(i,j)} is in logical layer `Ls[[1]]`
#'   and \eqn{(i-k)}{(i,k)} is in logical layer `Ls[[2]]` .
#'   2) *"out"* Number of configurations
#'   \eqn{(i{\rightarrow}j), (i{\rightarrow}k)}{(i,j), (i,k)}, where
#'   \eqn{(i{\rightarrow}j)}{(i,j)} is in logical layer `Ls[[1]]`
#'   and \eqn{(i{\rightarrow}k)}{(i,k)} is in logical layer `Ls[[2]]`.
#'   3) *"in"* Number of configurations
#'   \eqn{(j{\rightarrow}i), (k{\rightarrow}i)}{(j,i), (k,i)}, where
#'   \eqn{(j{\rightarrow}i)}{(j,i)} is in logical layer `Ls[[1]]`
#'   and \eqn{(k{\rightarrow}i)}{(k,i)} is in logical layer `Ls[[2]]`.
#'   4) *"path"* Number of configurations
#'   \eqn{(j{\rightarrow}i), (i{\rightarrow}k)}{(j,i), (i,k)}, where
#'   \eqn{(j{\rightarrow}i)}{(j,i)} is in logical layer `Ls[[1]]`
#'   and \eqn{(i{\rightarrow}k)}{(i,k)} is in logical layer `Ls[[2]]`.
#'
#'   At this time, `"any"` is only supported for undirected networks, and if the network is undirected, `type` is ignored and `"any"` is assumed.
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
                      vartypes = c("formula,list", "character", "logical"),
                      defaultvalues = list(NULL, NULL, TRUE),
                      required = c(TRUE, TRUE, FALSE))
  TYPES <- c("any", "out", "in", "path")
  TYPEREP <- setNames(c("--", "<>", "><", ">>"), TYPES)
  type <- match.arg(tolower(a$type), TYPES)
  if(!is.directed(nw)) type <- "any"
  else if(type == "any") ergm_Init_abort(paste0("at this time, ", sQuote('type="any"'), " is only supported for undirected networks"))
  typeID <- match(type, TYPES) - 1L

  Ls <- .set_layer_namemap(a$Ls, nw)
  if(is(Ls, "formula")) Ls <- list(Ls)
  Ls <- rep(Ls, length.out=2)

  auxiliaries <- .mk_.layer.net_auxform(Ls)
  reprs <- .lspec_coef.namewrap(Ls, collapse=FALSE)
  coef.names <- paste0("twostarL(",
                       paste0(reprs, collapse=TYPEREP[type]),
                       if(a$distinct) ",distinct",
                       ")")

  iinputs <- c(typeID, a$distinct)
  list(name="twostarL", coef.names=coef.names, iinputs=iinputs, auxiliaries=auxiliaries, minval=0, dependence=TRUE)
}
