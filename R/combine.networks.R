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
#' Given a list of compatible networks, [combine_networks()] returns a single
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
#' @param keep.unshared.attr whether to keep those network, vertex,
#'   and edge attributes not shared by all networks in the list; if
#'   \code{TRUE}, positions corresponding to networks lacking the
#'   attribute are replaced with \code{NA}, \code{NULL}, or some other
#'   placeholder; incompatible with \code{detect.edgecov==TRUE}.
#'
#' @return an object of class `combined_networks` inheriting from
#'   [`network::network`] with a block-diagonal structure (or its
#'   bipartite equivalent) comprising the networks passed in `nwl`. In
#'   particular,
#'
#' * the returned network's size is the sum of the input networks';
#'
#' * its basic properties (directedness and bipartednes) are the same;
#'
#' * the input networks' sociomatrices (both edge presence and edge
#'   attributes) are the blocks in the sociomatrix of the returned
#'   network;
#'
#' * vertex attributes are concatenated;
#'
#' * edge attributes are assigned to their respective edges in
#'   the returned network;
#'
#' * network attributes are stored in a list network attribute `".snattr"`.
#'
#' * a new network attribute, ".snid" identifies the original
#'   networks; if combined networks are combined again, ".snid"
#'   becomes a space-separated string of IDs from outside in.
#'
#' * a  new network  attribute, ".snl",  containing networks  in `nwl`
#'   (preserving their names, if any) but with all edges removed.
#'
#' 1. If `nwl` is a named list, names from the list are used.
#'
#' 2. If not 1, but the network has an attribute `title`, it is used.
#'
#' 3. Otherwise, a numerical index is used.
#'
#' @aliases combined_networks
#' @examples
#'
#' data(samplk)
#'
#' o1 <- combine_networks(list(samplk1, samplk2, samplk3))
#' image(as.matrix(o1))
#' head(get.vertex.attribute(o1, ".snid"))
#' o2 <- combine_networks(list(o1, o1))
#' image(as.matrix(o2))
#' head(get.vertex.attribute(o2, ".snid", unlist=FALSE))
#'
#' data(florentine)
#' f1 <- combine_networks(list(business=flobusiness, marriage=flomarriage))
#' image(as.matrix(f1))
#' head(get.vertex.attribute(f1, ".snid"))
#' head(names(get.network.attribute(f1, ".snl")))
#' @export
combine_networks <- function(nwl, ignore.nattr=c("mnext"), ignore.vattr=c(), ignore.eattr=c(), keep.unshared.attr=FALSE){
  if (!all_identical(map_lgl(nwl, is.directed))) stop("All networks must have the same directedness.")

  bm <- mk_block_maps(nwl)

  attrset <- if(keep.unshared.attr) union else intersect

  constructor <- if(is(nwl[[1]], "networkLite")) networkLite::networkLite else network.initialize
  out <-  constructor(sum(bm$ns), directed = is.directed(nwl[[1]]), bipartite = bm$bip)

  # Concatenate common network attributes.
  nattrs <- setdiff(reduce(map(nwl, list.network.attributes), attrset), ignore.nattr)
  out %n% ".snattr" <- map(nattrs, \(a) map(nwl, get.network.attribute, a, unlist = FALSE)) |> setNames(nattrs)

  # Concatenate vertex attributes.

  for (a in setdiff(reduce(map(nwl, list.vertex.attributes), attrset),
                    ignore.vattr)) { # I.e., iterate through common attributes.
    val <- map(nwl, get.vertex.attribute, a, unlist = FALSE)

    for (b in seq_along(val))
      out <- set.vertex.attribute(out, a, val[[b]], bm$v2V(b, seq_along(val[[b]])))
  }

  # Add ties and attributes
  el <- map(seq_along(nwl), function(b) {
    df <- as_tibble(nwl[[b]], attrname = TRUE, unit = "edge", na.rm = FALSE)
    names <- rep(list(as.list(names(df)[-(1:2)])), nrow(df))
    vals <- transpose(df[,-(1:2)])
    list(tails = bm$v2V(b, df[[1]]), heads = bm$v2V(b, df[[2]]), names = names, vals = vals)
  }) |> transpose()

  out <- add.edges(out, unlist(el$tails), unlist(el$heads), names.eval=do.call(c, el$names), vals.eval=do.call(c, el$vals))

  # Finally, add a vertex attribute specifying the blocks
  b <- bm$V2b(seq_len(network.size(out)))
  if(".snid" %in% list.vertex.attributes(out)){ # .snid already exists
    b <- map2(b, # each element of b
              get.vertex.attribute(out, ".snid"), # with the corresponding element of out %v% ".snid".
              paste)
  }
  out <- set.vertex.attribute(out, ".snid", b)

  out %n% ".snl" <- lapply(nwl, `[<-.network`, value = 0)

  if (is.null(names(nwl))) {
    titles <- map(nwl, `%n%`, "title") |> replace(is.null, "") |> unlist()
    names(nwl) <- titles
  }


  out %n% ".bm" <- bm

  structure(out, class = c("combined_networks", class(out)))
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

#' @rdname combine_networks
#' @description [uncombine_network()] undoes [combine_networks()], returning a list
#' of networks, optionally mapping the edges of the combined network
#' onto them.
#'
#' @param nw a [`network::network`] created by [combine_networks()].
#' @param populate logical; whether to map the edges from the combined
#'   network to the subnetworks, or whether to return them empty.
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
uncombine_network <- function(nw, populate = TRUE) {
  if(!is(nw, "combined_networks")) stop("Specified network was not constructed by ", sQuote("combine_networks()"), ".")
  nwl <- nw%n%".snl"
  bm <- nw%n%".bm"

  if (populate) {
    df <- as_tibble(nw, attrname = TRUE, unit = "edge", na.rm = FALSE)
    names <- rep(list(as.list(names(df)[-(1:2)])), nrow(df))
    vals <- transpose(df[,-(1:2)])
    el <- list(tails = bm$V2v(df[[1]]), heads = bm$V2v(df[[2]]), names = names, vals = vals) |> transpose()
    blk <- bm$V2b(df[[1]])
    ell <- split(el, factor(blk, levels = seq_along(nwl)))
    for (b in seq_along(nwl)) {
      bel <- transpose(ell[[b]])
      if (length(bel))
        nwl[[b]] <- add.edges(nwl[[b]], unlist(bel$tails), unlist(bel$heads),
                              names.eval = bel$names, vals.eval = bel$vals)
    }
  }

  nwl
}

#' Dynamic registry of functions that combine networks.
#'
#' This is used primarily by developers to provide informative error
#' messages.
#'
#' @param combiner the combiner name (typically the first element of
#'   `nw %n% ".combiner"`)
#'
#' @param description a named character vector or list with elements
#'   `constructor`, `element`, `elements`, `construct`, `id`, and
#'   `name`, containing, respectively, the top-level
#'   (i.e. user-visible) function that combines on it, what is being
#'   combined (singular and plural), and what columns, if any, should
#'   be used to identify its index and name in the attributes
#'   table. `elements` and `construct` will be extrapolated from the
#'   others if not given.
#'
#' @keywords internal
#' @export
ergm.multi_combiner <- local({
  cache <- list()

  function(combiner, description = NULL) {
    if (missing(combiner)) cache
    else if (is.null(description)) cache[[combiner]]
    else cache[[combiner]] <<- {
      description <- as.list(description)
      NVL(description[["elements"]]) <- paste0(description$element, "s")
      NVL(description[["construct"]]) <- paste0("multi-", description$element, " construct")
      description
    }
  }
})

#' Obtain empty networks representing constituents of a combined network
#'
#' This utility uncombines a [combine_networks()] network using subnetwork cache (which contains only empty networks). It is used primarily by initialisation functions.
#'
#' @param nw see [uncombine_network()].
#' @param copy.ergmlhs a character vector of [`%ergmlhs%`] settings
#'   that are to be copied into the constituent networks.
#' @param subset indices of the constituent networks to return.
#'
#' @return A list of [`network`]s.
#'
#' @keywords internal
#' @export
subnetwork_templates <- function(nw, copy.ergmlhs = c("response"), subset = TRUE) {
  nwl <- uncombine_network(nw, populate = FALSE)[subset]

  if (length(copy.ergmlhs)) nwl <- map(nwl, function(nw1) {
    for (name in copy.ergmlhs) nw1%ergmlhs%name <- nw%ergmlhs%name
    nw1
  })
  nwl
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
.block_vertexmap <- function(nw) {
  with(nw %n% ".bm",
       c(list(nb = length(ns), bids = bids, bmap = bmap),
         if (bip) list(ns = rbind(nes, nas))
         else list(ns = ns)))
}

#' @importFrom utils head
mk_block_maps <- function(nwl) {
  ns <- map_int(nwl, network.size) |> unname()
  n <- sum(ns)
  nes <- map_int(nwl, \(nw) b1.size(nw) %||% 0L) |> unname()

  info <- if (bip <- sum(nes)) {
    nas <- ns - nes

    eoff <- cumsum(c(0L, head(nes, -1L)))
    aoff <- cumsum(c(0L, head(nas, -1L))) - nes + sum(nes)

    list(bids = c(rep(seq_along(nes), nes),
                  rep(seq_along(nas), nas)),
         bmap = seq_len(n) - c(rep(eoff, nes),
                               rep(aoff, nas)),
         nes = nes,
         nas = nas,
         eoff = eoff,
         aoff = aoff)
  } else {
    off <- cumsum(c(0L, head(ns, -1L)))

    list(bids = rep(seq_along(ns), ns),
         bmap = seq_len(n) - rep(off, ns),
         off = off)
  }
  info <- list2env(c(info, list(n = n, bip = if (bip) bip else FALSE, ns = ns)),
                   parent = baseenv())

  c(as.list(info),
    local(list(
      V2b = function(v) bids[v],
      V2v = function(v) bmap[v],
      v2V = if (bip) function(b, v) ifelse(v <= nes[b], eoff[b] + v, aoff[b] + v)
            else function(b, v) off[b] + v), info))
}

#' @noRd
#' @importFrom networkLite as.networkLite
#' @export
as.networkLite.combined_networks <- function(x, ...){
  if(!is.null(x %n% ".snl")){
    snc <- x %n% ".snl"
    snc[[1]] <- lapply(snc[[1]], as.networkLite)
    x %n% ".snl" <- snc
  }
  extra_classes <- class(x)[seq_len(which(class(x) == "combined_networks")[1])]
  class(x) <- class(x)[-seq_len(which(class(x) == "combined_networks")[1])]
  x <- as.networkLite(x)
  structure(x, class = c(extra_classes, class(x)))
}

#' Confirm that this is a combined network.
#'
#' @param nw a [`network`]
#' @param combiners a character vector listing valid combiners supported by the function
#' @param term whether this is a term that should use [ergm_Init_stop()] or some other function, which should use [rlang::abort()]
#' @param call for `term_trace == FALSE` where to report the error message from
#'
#' @keywords internal
#' @export
assert_combined_network <- function(nw, combiners, term = TRUE, call = rlang::caller_env()) {
  if (!any(nw%n%".combiner" %in% combiners)) {
    info <- map(combiners, ergm.multi_combiner)
    valid <- paste.and(paste0(map_chr(info, "construct"), " constructed by ", sQuote(map_chr(info, "constructor"))),
                       con = "nor")
    msg <- paste0("This network is not a ", valid, " (nor compatible) at its top level.")
    if (term) ergm_Init_stop(msg)
    else abort(msg, call = call)
    FALSE
  } else TRUE
}
