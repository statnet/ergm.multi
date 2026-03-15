#  File R/combine.networks.R in package ergm.multi, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################

#' A single block-diagonal network created by combining multiple networks
#'
#' Given a list of compatible networks, the [combine_networks()] returns a single
#' block-diagonal network, preserving attributes that can be
#' preserved.
#'
#' @param nwl a list of [`network::network`]s to be combined; they
#'   must have similar fundamental properties: directedness and
#'   bipartedness, though their sizes (and the size of each bipartite
#'   group) can vary.
#'
#' @param ignore.nattr,ignore.vattr,ignore.eattr network, vertex, and
#'   edge attributes not to be processed as described below.
#'
#' @param blockID.vattr name of the vertex attribute into which to store
#'   the index of the network to which that vertex originally belonged.
#'
#' @param blockName.vattr if not `NULL`, the name of the vertex
#'   attribute into which to store the name of the network to which
#'   that vertex originally belonged.
#'
#' @param detect.edgecov if `TRUE`, combine network attributes that
#'   look like dyadic covariate ([`ergm::edgecov`][ergm::edgecov-ergmTerm]) matrices into a
#'   block-diagonal matrix.
#'
#' @param keep.unshared.attr whether to keep those network, vertex,
#'   and edge attributes not shared by all networks in the list; if
#'   \code{TRUE}, positions corresponding to networks lacking the
#'   attribute are replaced with \code{NA}, \code{NULL}, or some other
#'   placeholder; incompatible with \code{detect.edgecov==TRUE}.
#'
#' @param subnet.cache whether to save the input network list as an
#'   attribute of the combined network, so that if the network is
#'   resplit using on the same attribute (e.g. using
#'   [uncombine_network()]), an expensive call to [split.network()]
#'   can be avoided, at the cost of storage.
#'
#' @return an object of class `combined_networks` inheriting from
#'   [`network::network`], implemented as a list with two elements:
#'   `$nw` (the constituent networks) and `$gal` (network-level
#'   metadata), analogous to the [`networkLite::networkLite`] design.
#'   The combined object presents a block-diagonal network view (or its
#'   bipartite equivalent) of the input networks. In particular,
#'
#' * the returned network's size is the sum of the input networks';
#'
#' * its basic properties (directedness and bipartednes) are the same;
#'
#' * the input networks' sociomatrices (both edge presence and edge
#'   attributes) are the blocks when the object is converted to a
#'   matrix or a [`networkLite::networkLite`];
#'
#' * vertex attributes are presented as concatenated across constituent
#'   networks (with bipartite reordering if applicable);
#'
#' * network attributes are stored in a list in `$gal`; if
#'   `detect.edgecov==TRUE`, those network attributes that have the
#'   same dimension as the sociomatrices of the constituent networks
#'   are combined into a single block-diagonal matrix and stored there.
#'
#' Two virtual vertex attributes, specified by `blockID.vattr` and
#' (optionally) `blockName.vattr`, give the index and name of the
#' constituent network from which each vertex came. The name is
#' determined as follows:
#'
#' 1. If `nwl` is a named list, names from the list are used.
#'
#' 2. If not 1, but the network has an attribute `title`, it is used.
#'
#' 3. Otherwise, a numerical index is used.
#'
#' If `blockID.vattr` already exists on the constituent networks, the
#' index is *prepended* to the attribute.
#'
#' The values of `blockID.vattr` and `blockName.vattr` are stored in
#' network attributes `".blockID.vattr"` and `".blockName.vattr"`.
#'
#' @aliases combined_networks
#' @examples
#'
#' data(samplk)
#'
#' o1 <- combine_networks(list(samplk1, samplk2, samplk3))
#' image(as.matrix(o1))
#' head(get.vertex.attribute(o1, ".NetworkID"))
#' o2 <- combine_networks(list(o1, o1))
#' image(as.matrix(o2))
#' head(get.vertex.attribute(o2, ".NetworkID", unlist=FALSE))
#'
#' data(florentine)
#' f1 <- combine_networks(list(business=flobusiness, marriage=flomarriage),
#'                        blockName.vattr=".NetworkName")
#' image(as.matrix(f1))
#' head(get.vertex.attribute(f1, ".NetworkID"))
#' head(get.vertex.attribute(f1, ".NetworkName"))
#' @export
combine_networks <- function(nwl, ignore.nattr=c("mnext"), ignore.vattr=c(), ignore.eattr=c(), blockID.vattr=".NetworkID", blockName.vattr=NULL, detect.edgecov=FALSE, keep.unshared.attr=FALSE, subnet.cache=FALSE){
  if(any(sapply(nwl, is.bipartite)))
    .combine_networks.bipartite(nwl=nwl, ignore.nattr=ignore.nattr, ignore.vattr=ignore.vattr, ignore.eattr=ignore.eattr, blockID.vattr=blockID.vattr, blockName.vattr=blockName.vattr, detect.edgecov=detect.edgecov, keep.unshared.attr=keep.unshared.attr)
  else
    .combine_networks.unipartite(nwl=nwl, ignore.nattr=ignore.nattr, ignore.vattr=ignore.vattr, ignore.eattr=ignore.eattr, blockID.vattr=blockID.vattr, blockName.vattr=blockName.vattr, detect.edgecov=detect.edgecov, keep.unshared.attr=keep.unshared.attr)
}


.combine_networks.unipartite <- function(nwl, ignore.nattr=c("mnext"), ignore.vattr=c(), ignore.eattr=c(), blockID.vattr=".NetworkID", blockName.vattr=NULL, detect.edgecov=FALSE, keep.unshared.attr=FALSE){
  if(any(diff(sapply(nwl, is.directed)))) stop("All networks must have the same directedness.")
  if(keep.unshared.attr && detect.edgecov) stop("Detection of edge covariates is not compatible with retaining unshared attributes.")
  attrset <- if(keep.unshared.attr) union else intersect

  ns <- sapply(nwl, network.size)
  blks <- c(0, cumsum(ns))

  # Build .subnetattr metadata
  sna <- list()
  sna[[blockID.vattr]] <- list()

  # gal starts with required network properties
  gal <- list(
    n = sum(ns),
    directed = is.directed(nwl[[1]]),
    bipartite = FALSE,
    loops = FALSE,
    hyper = FALSE,
    multiple = FALSE,
    ".blockID.vattr" = blockID.vattr,
    ".blockName.vattr" = blockName.vattr
  )

  for(a in setdiff(Reduce(attrset,lapply(nwl, list.network.attributes)),
                          ignore.nattr)){ # I.e., iterate through common attributes.
    sna[[blockID.vattr]][[a]] <- vl <- lapply(nwl, get.network.attribute, a, unlist=FALSE)

    # Here, try to autodetect covariate matrices and combine them.
    if(detect.edgecov
       && all(sapply(vl, is.matrix))
       && all(sapply(vl, nrow)==ns)
       && all(sapply(vl, ncol)==ns)
       && all_identical(sapply(vl, mode))){

      # A logical vector that extracts off-diagonal element of the ns*ns matrix.


      offdiags <- unlist(lapply(ns, function(n) c(diag(1,n)==0)))
      # It doesn't matter what the "filler" elements are, as long as
      # adding them doesn't add another category and it's not NA. So,
      # what this does is as follows: grab the off-diagonal elements
      # of each covariate matrix, concatenate them into one vector,
      # remove the NAs, and take 0 (if it's present) or the minimum
      # value. (0 as filler can be helpful for sparse matrix
      # representations.)
      dummyvals <- na.omit(unlist(lapply(vl, "c"))[offdiags])
      dummyval <- if(0 %in% dummyvals) 0 else min(dummyvals)
      m <- matrix(dummyval, sum(ns), sum(ns))
      mode(m) <- mode(vl[[1]])

      for(b in seq_along(vl)){
        inds <- blks[b]+seq_len(ns[b])
        m[inds, inds] <- vl[[b]]
      }

      gal[[a]] <- m
    }
  }
  sna[[blockID.vattr]]$..class <- lapply(nwl, class)
  gal[[".subnetattr"]] <- sna

  # Determine block names
  bn <- if(!is.null(names(nwl))) names(nwl)
        else if("title" %in% Reduce(intersect, lapply(nwl, list.network.attributes)))
          sapply(nwl, `%n%`, "title")
        else as.character(seq_along(ns))
  gal[[".blocknames"]] <- bn

  structure(
    list(nw = nwl, gal = gal),
    class = c("combined_networks", "network")
  )
}


.combine_networks.bipartite <- function(nwl, ignore.nattr=c("mnext"), ignore.vattr=c(), ignore.eattr=c(), blockID.vattr=".NetworkID", blockName.vattr=NULL, detect.edgecov=FALSE, keep.unshared.attr=FALSE){
  if(!all(sapply(nwl, is.bipartite))) stop("This function operates only on bipartite networks.")
  if(any(sapply(nwl, is.directed))) stop("Bipartite directed networks are not supported at this time.")
  if(keep.unshared.attr && detect.edgecov) stop("Detection of edge covariates is not compatible with retaining unshared attributes.")
  attrset <- if(keep.unshared.attr) union else intersect

  ns <- sapply(nwl, network.size)
  es <- sapply(nwl, b1.size)
  eblks <- c(0, cumsum(es))
  bip <- eblks[length(eblks)]
  ablks <- cumsum(c(bip, ns-es))

  # Build .subnetattr metadata
  sna <- list()
  sna[[blockID.vattr]] <- list()

  # gal starts with required network properties
  gal <- list(
    n = sum(ns),
    directed = FALSE,
    bipartite = bip,
    loops = FALSE,
    hyper = FALSE,
    multiple = FALSE,
    ".blockID.vattr" = blockID.vattr,
    ".blockName.vattr" = blockName.vattr
  )

  for(a in setdiff(Reduce(attrset,lapply(nwl, list.network.attributes)),
                          ignore.nattr)){ # I.e., iterate through common attributes.
    sna[[blockID.vattr]][[a]] <- vl <- lapply(nwl, get.network.attribute, a, unlist=FALSE)

    # Here, try to autodetect covariate matrices and combine them.
    if(detect.edgecov
       && all(sapply(vl, is.matrix))
       && all(sapply(vl, nrow)==es)
       && all(sapply(vl, ncol)==ns-es)
       && all_identical(sapply(vl, mode))){

      # It doesn't matter what the "filler" elements are, as long as
      # adding them doesn't add another category and it's not NA. So,
      # what this does is as follows: grab the elements of each
      # covariate matrix, concatenate them into one vector, remove the
      # NAs, and take 0 (if it's present) or the minimum value. (0 as
      # filler can be helpful for sparse matrix representations.)
      dummyvals <- na.omit(unlist(lapply(vl, "c")))
      dummyval <- if(0 %in% dummyvals) 0 else min(dummyvals)
      m <- matrix(dummyval, sum(es), sum(ns-es))
      mode(m) <- mode(vl[[1]])

      for(b in seq_along(vl)){
        einds <- eblks[b]+seq_len(es[b])
        ainds <- ablks[b]+seq_len(ns[b]-es[b])
        m[einds, ainds-sum(es)] <- vl[[b]]
      }

      gal[[a]] <- m
    }
  }
  sna[[blockID.vattr]]$..class <- lapply(nwl, class)
  gal[[".subnetattr"]] <- sna

  # Determine block names
  bn <- if(!is.null(names(nwl))) names(nwl)
        else if("title" %in% Reduce(intersect, lapply(nwl, list.network.attributes)))
          sapply(nwl, `%n%`, "title")
        else as.character(seq_along(ns))
  gal[[".blocknames"]] <- bn

  structure(
    list(nw = nwl, gal = gal),
    class = c("combined_networks", "network")
  )
}

#' In a combined network, obtain the top-level attribute vector on which it was combined.
#'
#' @param x a [combine_networks()] network
#' @param attrname a string
#' @param missing what to return if the attribute is not found
#'
#' @keywords internal
#' @export
get_combining_attr <- function(x, attrname, missing = c("NA", "NULL")) {
  missing <- match.arg(missing)
  if (missing == "NULL"
      && ! attrname %in% list.vertex.attributes(x)) return(NULL)

  av <- get.vertex.attribute(x, attrname, unlist = FALSE)
  sapply(av, "[", 1)
}

# Used only for the legacy (networkLite-based) uncombine path.
.pop_vattrv <- function(nw, vattr){
  av <- get.vertex.attribute(nw, vattr, unlist=FALSE)
  a <- sapply(av, "[", 1)
  rest <- lapply(av, "[", -1)

  if(all(lengths(rest)==0)) delete.vertex.attribute(nw, vattr)
  else set.vertex.attribute(nw, vattr, rest)

  list(nw = nw, vattr = a)
}


#' A [split()] method for [`network::network`] objects.
#'
#' Split a network into subnetworks on a factor.
#'
#' @param x a [`network::network`] object.
#'
#' @param f,drop,sep,lex.order see [split()]; note that `f` must have length equal to `network.size(x)`.
#'
#' @param ... additional arguments, currently unused.
#'
#' @return A [`network.list`] containing the networks. These networks
#'   will inherit all vertex and edge attributes, as well as relevant
#'   network attributes.
#'
#' @seealso [network::get.inducedSubgraph()]
#' @export
split.network <- function(x, f, drop = FALSE, sep = ".", lex.order = FALSE, ...)
{
  ### NOTE: This is taken from the split.default() implementation, but is trivial.
  if(!missing(...))
    .NotYetUsed(deparse1(...), error = FALSE)
  if(is.list(f))
    f <- interaction(f, drop = drop, sep = sep, lex.order = lex.order)
  else if (!is.factor(f))
    f <- as.factor(f)
  else if (drop)
    f <- factor(f)
  ### END Taken from split.default().

  o <- lapply(levels(f), function(l) network::get.inducedSubgraph(x, which(f==l)))
  class(o) <- c("network.list", class(o))
  o
}

#' Split up a network into a list of subgraphs
#'
#' Given a network created by [combine_networks()], [uncombine_network()] returns a list of networks,
#' preserving attributes that can be preserved.
#'
#' @param nw a [`network::network`] created by [combine_networks()].
#'
#' @param split.vattr name of the vertex attribute on which to split,
#'   defaulting to the value of the `".blockID.vattr"` network
#'   attribute.
#'
## #' @param detect.edgecov if `TRUE`, split up network attributes that
## #'   look like dyadic covariate ([`ergm::edgecov`][ergm::edgecov-ergmTerm]) matrices.
#'
#' @param names.vattr optional name of the vertex attribute to use as
#'   network names in the output list, defaulting to the value of the
#'   `".blockName.vattr"` network attribute.
#'
#' @param use.subnet.cache whether to use subnet cache if available;
#'   this is only safe to do if the network is *not* used for its
#'   edges but only for its vertex and network attributes.
#'
#' @return a list of [`network::network`]s containing subgraphs on `split.vattr`. In particular,
#'
#' * their basic properties (directedness and bipartednes) are the same as those of the input network;
#'
#' * vertex attributes are split;
#'
#' * edge attributes are assigned to their respective edges in
#'   the returned networks.
#'
#' If `split.vattr` is a vector, only the first element is used and it's "popped".
#'
#' @seealso [split.network()]
#' @examples
#'
#' data(samplk)
#'
#' o1 <- combine_networks(list(samplk1, samplk2, samplk3))
#' image(as.matrix(o1))
#'
#' ol <- uncombine_network(o1)
#'
#' @export
uncombine_network <- function(nw, split.vattr=nw %n% ".blockID.vattr", names.vattr=nw %n% ".blockName.vattr", use.subnet.cache=FALSE){
  if(!is(nw, "combined_networks")) stop("Specified network was not constructed by ", sQuote("combine_networks()"), ".")

  if(!is.null(nw$nw)){
    # New list-based combined_networks: constituent networks are stored in nw$nw.
    # They are already complete objects (with edges and attributes), so just
    # return them directly -- no splitting or attribute restoration needed.
    nwl <- nw$nw
    if(!is.null(names.vattr)) names(nwl) <- nw$gal[[".blocknames"]]
  } else {
    # Legacy / post-simulation networkLite-based combined_networks:
    # split the block-diagonal by copying edges using known block offsets.
    tmp <- .pop_vattrv(nw, split.vattr); nw <- tmp$nw; f <- tmp$vattr
    if(!is.null(names.vattr)){ tmp <- .pop_vattrv(nw, names.vattr); nw <- tmp$nw; nwnames <- tmp$vattr }

    sna <- (nw %n% ".subnetattr")[[split.vattr]]

    # Use subnet cache (empty templates) + copy edges from block-diagonal.
    templates <- if(!is.null(nw %n% ".subnetcache") && names(nw %n% ".subnetcache") == split.vattr)
      (nw %n% ".subnetcache")[[split.vattr]]
    else
      NULL

    nwl <- .uncombine_network_from_blockdiag(nw, f, sna, templates)

    if(!is.null(names.vattr)) names(nwl) <- unique(nwnames)
  }

  class(nwl) <- c("network.list", class(nwl))
  nwl
}

# Helper: reconstruct constituent networks from a block-diagonal network/networkLite.
# Uses known block offsets to copy edges, avoiding get.inducedSubgraph().
.uncombine_network_from_blockdiag <- function(nw, f, sna, templates=NULL){
  n <- network.size(nw)
  bip <- b1.size(nw)
  levels <- unique(f)
  nb <- length(levels)

  # Determine per-block vertex ranges in the combined ordering.
  ns <- sna$n
  if(NVL(bip, 0)){
    es <- sna$bipartite
    eblks <- c(0, cumsum(es))
    ablks <- cumsum(c(bip, ns - es))
    # For bipartite, vertices for block i are at:
    #   b1: eblks[i]+1 .. eblks[i+1]
    #   b2: ablks[i]+1 .. ablks[i+1]  (ablks indexed starting from 1 = first b2 block)
    ranges <- lapply(seq_len(nb), function(i)
      c(eblks[i] + seq_len(es[i]),
        ablks[i] + seq_len(ns[i] - es[i])))
    offsets <- lapply(seq_len(nb), function(i)
      c(rep(eblks[i], es[i]),
        rep(ablks[i] - es[i], ns[i] - es[i])))
  } else {
    blks <- c(0, cumsum(ns))
    ranges <- lapply(seq_len(nb), function(i) blks[i] + seq_len(ns[i]))
    offsets <- lapply(seq_len(nb), function(i) rep(blks[i], ns[i]))
  }

  # Get the full edge list of the block-diagonal network.
  el_full <- as_tibble(nw, attrnames = TRUE, unit = "edge", na.rm = FALSE)
  eattr_names <- names(el_full)[-(1:2)]

  nwl <- vector("list", nb)
  for(i in seq_len(nb)){
    r <- ranges[[i]]
    off <- offsets[[i]]

    # Build (or reuse) a template for this constituent.
    if(!is.null(templates) && i <= length(templates)){
      nwi <- templates[[i]]
    } else {
      constructor <- if(is(nw, "networkLite")) networkLite::networkLite else network.initialize
      if(NVL(bip, 0)){
        nwi <- constructor(ns[i], directed = FALSE, bipartite = es[i])
      } else {
        nwi <- constructor(ns[i], directed = is.directed(nw))
      }
      # Copy vertex attributes
      for(va in setdiff(list.vertex.attributes(nw),
                        c(list.network.attributes(nw), ".NetworkID", ".NetworkName",
                          ".LayerID", ".LayerName"))){
        vals <- get.vertex.attribute(nw, va, unlist = FALSE)
        if(!is.null(vals)) nwi <- set.vertex.attribute(nwi, va, vals[r])
      }
    }

    # Copy edges that belong to this block.
    if(nrow(el_full) > 0){
      tail_in <- el_full[[1]] %in% r
      head_in <- el_full[[2]] %in% r
      el_i <- el_full[tail_in & head_in, , drop = FALSE]
      if(nrow(el_i) > 0){
        local_tails <- el_i[[1]] - off[match(el_i[[1]], r)]
        local_heads <- el_i[[2]] - off[match(el_i[[2]], r)]
        if(length(eattr_names)){
          enames <- rep(list(as.list(eattr_names)), nrow(el_i))
          evals  <- transpose(el_i[, eattr_names, drop = FALSE])
          nwi <- add.edges(nwi, local_tails, local_heads,
                           names.eval = enames, vals.eval = evals)
        } else {
          nwi <- add.edges(nwi, local_tails, local_heads)
        }
      }
    }

    # Restore network attributes from .subnetattr.
    class(nwi) <- sna$..class[[i]]
    for(nattr in c(".subnetattr", ".blockID.vattr", ".blockName.vattr", ".blocknames"))
      if(nattr %in% list.network.attributes(nwi))
        nwi <- delete.network.attribute(nwi, nattr)
    for(nattr in setdiff(names(sna), c("..class", "n", "directed", "bipartite", "loops")))
      nwi %n% nattr <- sna[[nattr]][[i]]

    nwl[[i]] <- nwi
  }
  nwl
}
#' Dynamic registry of functions that combine networks.
#'
#' This is used primarily by developers to provide informative error
#' messages.
#'
#' @param blockID the vertex attribute used as the ID of the block.
#'
#' @param combiner a character vector of length 1 or 2, giving the
#'   top-level (i.e. user-visible) function that combines on it and
#'   optionally what is being combined.
#'
#' @keywords internal
#' @export
ergm.multi_combiner <- local({
  cache <- list()

  function(blockID, combiner=NULL){
    if(missing(blockID)) cache
    else if(is.null(combiner)) cache[[blockID]]
    else cache[[blockID]] <<- c(combiner, tolower(combiner))
  }
})

#' Obtain empty networks representing constituents of a combined network
#'
#' This utility uncombines a [combine_networks()] network using subnetwork cache (which contains only empty networks). It is used primarily by initialisation functions.
#'
#' @param nw,split.vattr,names.vattr see [uncombine_network()].
#'
#' @param copy.ergmlhs a character vector of [`%ergmlhs%`] settings that are to be copied into the constituent networks.
#'
#' @return A list of [`network`]s.
#'
#' @keywords internal
#' @export
subnetwork_templates <- function(nw, split.vattr=nw%n%".blockID.vattr", names.vattr=nw%n%".blockName.vattr", copy.ergmlhs=c("response")){
  if(NVL3(nw.split.vattr<-nw%n%".blockID.vattr", split.vattr != ., FALSE)){
    ergm_Init_stop(
      "The LHS was (at the top level) created by ",
      NVL3(ergm.multi_combiner(nw.split.vattr)[1],
           sQuote(.),
           paste0("an unknown function that uses attribute ", sQuote(nw.split.vattr))),
      " but the term is trying to extract its ",
      NVL3(ergm.multi_combiner(split.vattr)[2],
           paste0(., " (created by ", sQuote(ergm.multi_combiner(split.vattr)[1]), ")."),
           paste0("subgraphs defined by ", split.vattr)),
      " Nesting of terms must match the nesting of constructors; this may change in the future.", immediate.=TRUE)
  }

  if(!is.null(nw$nw)){
    # New list-based: return constituent networks with edges cleared.
    nwl <- lapply(nw$nw, .clear_edges)
    if(!is.null(names.vattr)) names(nwl) <- nw$gal[[".blocknames"]]
  } else {
    # Legacy / networkLite-based: use subnet cache if available.
    nwl <- uncombine_network(nw, split.vattr = split.vattr, names.vattr = names.vattr, use.subnet.cache = TRUE)
  }

  if (length(copy.ergmlhs)) nwl <- map(nwl, function(nw1) {
    for (name in copy.ergmlhs) nw1%ergmlhs%name <- nw%ergmlhs%name
    nw1
  })
  else nwl
}

#' Calculate a vector that maps the combined (block-diagonal) LHS network Vertex indices within-layer/within-network Vertex and a Vertex to layer/network lookup table.
#'
#' @param nw combined network.
#' @param by vertex attribute on which to split blocks.
#' @param same_dim whether all blocks must have the same dimensions (usually `FALSE` for multinetwork and `TRUE` for multilayer objects).
#'
#' @return A list with the following elements:
#'
#' \item{nb}{Number of blocks.}
#'
#' \item{bids}{A vector equal in length to the size of the combined
#' network containing the 1-based IDs of the block to which each
#' vertex belongs.}
#'
#' \item{bmap}{A vector equal in length to the size of the combined
#' network containing the 1-based IDs of the each
#' vertex's within-block ID.}
#'
#' @noRd
.block_vertexmap <- function(nw, by=nw %n% ".blockID.vattr", same_dim=FALSE){
  a <- get_combining_attr(nw, by)
  n <- length(a)
  bip <- b1.size(nw)
  if(NVL(bip,0)){
    ea <- a[seq_len(bip)]
    aa <- a[bip+seq_len(n-bip)]
    el <- rle(ea)$lengths
    al <- rle(aa)$lengths
    if(same_dim) if(!all_identical(el) || !all_identical(al)) stop("Layers must be networks of the same dimensions.", call.=FALSE)

    eoff <- rep(cumsum(c(0,el[-length(el)])), el)
    aoff <- rep(cumsum(c(0,al[-length(al)]))-el, al) + sum(el)

    o <- list(nb = length(el), bids = a, bmap = seq_len(n) - c(eoff,aoff), ns = rbind(el,al))
  }else{
    l <- rle(a)$lengths
    if(same_dim) if(!all_identical(l)) stop("Layers must be networks of the same size.", call.=FALSE)
    off <- rep(cumsum(c(0,l[-length(l)])), l)
    o <- list(nb = length(l), bids = a, bmap = seq_len(n) - off, ns = l)
  }
  o
}

#' @noRd
#' @importFrom networkLite as.networkLite
#' @export
as.networkLite.combined_networks <- function(x, ...){
  if(!is.null(x$nw)){
    # New list-based: build block-diagonal networkLite from constituent networks.
    nwl <- x$nw
    ns  <- sapply(nwl, network.size)
    bip <- b1.size(x)
    blockID.vattr   <- x$gal[[".blockID.vattr"]]
    blockName.vattr <- x$gal[[".blockName.vattr"]]

    if(NVL(bip, 0)){
      es     <- sapply(nwl, b1.size)
      eblks  <- c(0, cumsum(es))
      ablks  <- cumsum(c(bip, ns - es))

      out <- networkLite::networkLite(sum(ns), directed = FALSE, bipartite = bip)

      # Vertex attributes (bipartite order: b1 of all networks, then b2).
      vattrnames <- setdiff(Reduce(union, lapply(nwl, list.vertex.attributes)),
                            c(blockID.vattr, blockName.vattr))
      for(va in vattrnames){
        vl <- lapply(nwl, get.vertex.attribute, va, unlist = FALSE)
        v <- vector("list", sum(ns))
        for(b in seq_along(nwl)){
          einds <- eblks[b] + seq_len(es[b])
          ainds <- ablks[b] + seq_len(ns[b] - es[b])
          v[einds] <- vl[[b]][seq_len(es[b])]
          v[ainds] <- vl[[b]][es[b] + seq_len(ns[b] - es[b])]
        }
        out <- set.vertex.attribute(out, va, v)
      }

      # Block-ID vertex attribute.
      bid <- rep(rep(seq_along(nwl), 2), c(es, ns - es))
      bname <- if(!is.null(blockName.vattr)){
        bn <- x$gal[[".blocknames"]]
        rep(rep(bn, 2), c(es, ns - es))
      }

      # Edge list with remapped indices.
      el <- map(seq_along(nwl), function(b) {
        df <- as_tibble(nwl[[b]], attrname = TRUE, unit = "edge", na.rm = FALSE)
        remap <- function(i) i + ifelse(i <= es[b], eblks[b], -es[b] + ablks[b])
        nms  <- rep(list(as.list(names(df)[-(1:2)])), nrow(df))
        vals <- transpose(df[, -(1:2)])
        list(tails = remap(df[[1]]), heads = remap(df[[2]]), names = nms, vals = vals)
      }) %>% transpose()
    } else {
      blks <- c(0, cumsum(ns))

      out <- networkLite::networkLite(sum(ns), directed = is.directed(nwl[[1]]))

      # Vertex attributes.
      vattrnames <- setdiff(Reduce(union, lapply(nwl, list.vertex.attributes)),
                            c(blockID.vattr, blockName.vattr))
      for(va in vattrnames){
        out <- set.vertex.attribute(out, va,
                                    do.call(c, lapply(nwl, get.vertex.attribute, va, unlist = FALSE)))
      }

      # Block-ID vertex attribute.
      bid <- rep(seq_along(nwl), ns)
      bname <- if(!is.null(blockName.vattr)){
        bn <- x$gal[[".blocknames"]]
        rep(bn, ns)
      }

      # Edge list with remapped indices.
      el <- map(seq_along(nwl), function(b) {
        df <- as_tibble(nwl[[b]], attrname = TRUE, unit = "edge", na.rm = FALSE)
        nms  <- rep(list(as.list(names(df)[-(1:2)])), nrow(df))
        vals <- transpose(df[, -(1:2)])
        list(tails = df[[1]] + blks[b], heads = df[[2]] + blks[b], names = nms, vals = vals)
      }) %>% transpose()
    }

    out <- add.edges(out, unlist(el$tails), unlist(el$heads),
                     names.eval = do.call(c, el$names), vals.eval = do.call(c, el$vals))

    # Add block-ID and (optionally) block-name vertex attributes.
    out <- set.vertex.attribute(out, blockID.vattr, .as_nested_bid(x, bid))
    if(!is.null(blockName.vattr))
      out <- set.vertex.attribute(out, blockName.vattr, .as_nested_bname(x, bname))

    # Copy network-level attributes from gal (excluding internals).
    internal_gal <- c("n", "directed", "bipartite", "loops", "hyper", "multiple",
                       ".blockID.vattr", ".blockName.vattr", ".subnetattr", ".blocknames")
    for(a in setdiff(names(x$gal), internal_gal))
      out <- set.network.attribute(out, a, x$gal[[a]])

    # Copy combined_networks metadata.
    out <- set.network.attribute(out, ".blockID.vattr",   blockID.vattr)
    out <- set.network.attribute(out, ".blockName.vattr", blockName.vattr)
    out <- set.network.attribute(out, ".subnetattr",      x$gal[[".subnetattr"]])

    # Store empty-network subnet cache for efficient post-simulation uncombine.
    # Convert templates to networkLite to ensure they have proper vertex
    # attributes (including nested block IDs) for multi-level uncombining.
    snc <- list()
    snc[[blockID.vattr]] <- lapply(nwl, function(nw1) as.networkLite(.clear_edges(nw1)))
    out <- set.network.attribute(out, ".subnetcache", snc)

    extra_classes <- class(x)[seq_len(which(class(x) == "combined_networks")[1])]
    structure(out, class = c(extra_classes, class(out)))
  } else {
    # Already a networkLite-based combined_networks (e.g., post-simulation).
    extra_classes <- class(x)[seq_len(which(class(x) == "combined_networks")[1])]
    class(x) <- class(x)[-seq_len(which(class(x) == "combined_networks")[1])]
    x <- as.networkLite(x)
    structure(x, class = c(extra_classes, class(x)))
  }
}

# Helpers for building nested block-ID / block-name vertex attributes.
.as_nested_bid <- function(x, top_bid){
  blockID.vattr <- x$gal[[".blockID.vattr"]]
  # Check whether constituent networks already have a blockID.vattr (nested case).
  any_nested <- any(sapply(x$nw, function(nw1) blockID.vattr %in% list.vertex.attributes(nw1)))
  if(!any_nested) return(top_bid)

  bip <- b1.size(x)
  ns  <- sapply(x$nw, network.size)
  es  <- if(NVL(bip,0)) sapply(x$nw, b1.size) else ns

  # Build the combined nested-bid list in block-diagonal vertex order.
  if(NVL(bip,0)){
    inner <- lapply(x$nw, function(nw1){
      v <- get.vertex.attribute(nw1, blockID.vattr, unlist = FALSE)
      v1 <- v[seq_len(b1.size(nw1))]
      v2 <- v[b1.size(nw1) + seq_len(network.size(nw1) - b1.size(nw1))]
      list(b1 = v1, b2 = v2)
    })
    inner_flat <- c(
      do.call(c, lapply(inner, `[[`, "b1")),
      do.call(c, lapply(inner, `[[`, "b2"))
    )
  } else {
    inner_flat <- do.call(c, lapply(x$nw, function(nw1)
      get.vertex.attribute(nw1, blockID.vattr, unlist = FALSE)))
  }

  mapply(function(outer, inner) c(outer, inner), top_bid, inner_flat,
         SIMPLIFY = FALSE)
}

.as_nested_bname <- function(x, top_bname){
  blockName.vattr <- x$gal[[".blockName.vattr"]]
  if(is.null(blockName.vattr)) return(top_bname)
  any_nested <- any(sapply(x$nw, function(nw1) blockName.vattr %in% list.vertex.attributes(nw1)))
  if(!any_nested) return(top_bname)

  bip <- b1.size(x)
  if(NVL(bip,0)){
    inner <- lapply(x$nw, function(nw1){
      v <- get.vertex.attribute(nw1, blockName.vattr, unlist = FALSE)
      list(b1 = v[seq_len(b1.size(nw1))],
           b2 = v[b1.size(nw1) + seq_len(network.size(nw1) - b1.size(nw1))])
    })
    inner_flat <- c(do.call(c, lapply(inner, `[[`, "b1")),
                    do.call(c, lapply(inner, `[[`, "b2")))
  } else {
    inner_flat <- do.call(c, lapply(x$nw, function(nw1)
      get.vertex.attribute(nw1, blockName.vattr, unlist = FALSE)))
  }
  mapply(function(outer, inner) c(outer, inner), top_bname, inner_flat,
         SIMPLIFY = FALSE)
}
