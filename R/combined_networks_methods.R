#  File R/combined_networks_methods.R in package ergm.multi, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################

# S3 methods that make the new list-based combined_networks behave like a
# network object.  The internal structure is:
#
#   list(
#     nw  = list of constituent network objects,
#     gal = named list of network-level attributes (analogous to networkLite$gal)
#   )
#   class = c("combined_networks", "network")
#
# Most methods either delegate to x$gal (network-level attributes) or
# compute/aggregate across x$nw (vertex-level attributes, size, etc.).

# ---------------------------------------------------------------------------
# Network attribute methods
# ---------------------------------------------------------------------------

#' @export
get.network.attribute.combined_networks <- function(x, attrname, unlist = FALSE, ...) {
  out <- x$gal[[attrname]]
  if(unlist) out <- base::unlist(out)
  out
}

#' @export
set.network.attribute.combined_networks <- function(x, attrname, value, ...) {
  x$gal[[attrname]] <- value
  x
}

#' @export
list.network.attributes.combined_networks <- function(x, ...) {
  sort(unique(names(x$gal)))
}

#' @export
delete.network.attribute.combined_networks <- function(x, attrname, ...) {
  for(a in attrname) x$gal[[a]] <- NULL
  x
}

# ---------------------------------------------------------------------------
# Network-level property accessors (delegating to gal)
# ---------------------------------------------------------------------------

#' @export
network.size.combined_networks <- function(x, ...) x$gal[["n"]]

#' @export
is.directed.combined_networks <- function(x) x$gal[["directed"]]

# has.loops, has.hyper, is.multiplex are not S3 generics from network, but
# they call get.network.attribute internally, which we have overridden.
# Similarly for is.bipartite:
#   is.bipartite(x) uses x %n% "bipartite", which dispatches via our
#   get.network.attribute method.

# ---------------------------------------------------------------------------
# Vertex attribute methods
# ---------------------------------------------------------------------------

#' @export
list.vertex.attributes.combined_networks <- function(x, ...) {
  blockID.vattr   <- x$gal[[".blockID.vattr"]]
  blockName.vattr <- x$gal[[".blockName.vattr"]]

  base_attrs <- sort(unique(Reduce(union, lapply(x$nw, list.vertex.attributes))))
  extra <- c(blockID.vattr, blockName.vattr)
  sort(unique(c(base_attrs, extra[!is.null(extra)])))
}

#' @export
get.vertex.attribute.combined_networks <- function(x, attrname,
                                                    unlist = TRUE,
                                                    null.na = TRUE,
                                                    na.omit = FALSE,
                                                    ...) {
  blockID.vattr   <- x$gal[[".blockID.vattr"]]
  blockName.vattr <- x$gal[[".blockName.vattr"]]
  bip <- b1.size(x)

  if(!is.null(blockID.vattr) && attrname == blockID.vattr) {
    return(.get_blockID_vattr(x, unlist = unlist))
  }
  if(!is.null(blockName.vattr) && attrname == blockName.vattr) {
    return(.get_blockName_vattr(x, unlist = unlist))
  }

  # Concatenate the attribute from each constituent network.
  # Bipartite ordering: all b1 vertices first, then all b2 vertices.
  if(NVL(bip, 0)){
    nwl <- x$nw
    ns  <- sapply(nwl, network.size)
    es  <- sapply(nwl, b1.size)
    b1_vals <- lapply(nwl, function(nw1) {
      v <- get.vertex.attribute(nw1, attrname, unlist = FALSE, null.na = null.na)
      v[seq_len(b1.size(nw1))]
    })
    b2_vals <- lapply(nwl, function(nw1) {
      v <- get.vertex.attribute(nw1, attrname, unlist = FALSE, null.na = null.na)
      v[b1.size(nw1) + seq_len(network.size(nw1) - b1.size(nw1))]
    })
    vals_list <- c(do.call(c, b1_vals), do.call(c, b2_vals))
  } else {
    vals_list <- do.call(c,
      lapply(x$nw, get.vertex.attribute, attrname, unlist = FALSE, null.na = null.na))
  }

  if(na.omit) vals_list <- vals_list[!sapply(vals_list, is.na)]

  if(unlist) base::unlist(vals_list, recursive = FALSE) else vals_list
}

# Compute the block-ID vertex attribute for the combined network.
# For nested combined networks the IDs are prepended to the inner IDs.
.get_blockID_vattr <- function(x, unlist = TRUE) {
  blockID.vattr <- x$gal[[".blockID.vattr"]]
  nwl <- x$nw
  ns  <- sapply(nwl, network.size)
  bip <- b1.size(x)

  if(NVL(bip, 0)){
    es <- sapply(nwl, b1.size)
    top_bid <- rep(rep(seq_along(nwl), 2), c(es, ns - es))
  } else {
    top_bid <- rep(seq_along(nwl), ns)
  }

  # Check for nested block IDs in constituent networks.
  has_nested <- sapply(nwl, function(nw1)
    blockID.vattr %in% list.vertex.attributes(nw1))

  if(!any(has_nested)){
    if(unlist) return(top_bid)
    return(as.list(top_bid))
  }

  # Build combined nested list in block-diagonal vertex order.
  if(NVL(bip, 0)){
    inner_parts <- lapply(nwl, function(nw1){
      if(blockID.vattr %in% list.vertex.attributes(nw1)){
        v <- get.vertex.attribute(nw1, blockID.vattr, unlist = FALSE)
        list(b1 = v[seq_len(b1.size(nw1))],
             b2 = v[b1.size(nw1) + seq_len(network.size(nw1) - b1.size(nw1))])
      } else {
        n1 <- network.size(nw1); e1 <- b1.size(nw1)
        list(b1 = as.list(rep(NA, e1)), b2 = as.list(rep(NA, n1 - e1)))
      }
    })
    inner_flat <- c(
      do.call(c, lapply(inner_parts, `[[`, "b1")),
      do.call(c, lapply(inner_parts, `[[`, "b2"))
    )
  } else {
    inner_flat <- do.call(c, lapply(nwl, function(nw1){
      if(blockID.vattr %in% list.vertex.attributes(nw1))
        get.vertex.attribute(nw1, blockID.vattr, unlist = FALSE)
      else
        as.list(rep(NA, network.size(nw1)))
    }))
  }

  result <- mapply(function(outer, inner)
    if(is.na(inner[[1]])) list(outer) else c(list(outer), inner),
    top_bid, inner_flat, SIMPLIFY = FALSE)

  if(unlist) sapply(result, `[`, 1) else result
}

# Compute the block-name vertex attribute for the combined network.
.get_blockName_vattr <- function(x, unlist = TRUE) {
  blockName.vattr <- x$gal[[".blockName.vattr"]]
  bn  <- x$gal[[".blocknames"]]
  nwl <- x$nw
  ns  <- sapply(nwl, network.size)
  bip <- b1.size(x)

  if(NVL(bip, 0)){
    es <- sapply(nwl, b1.size)
    top_bname <- rep(rep(bn, 2), c(es, ns - es))
  } else {
    top_bname <- rep(bn, ns)
  }

  # Check for nested block names in constituent networks.
  has_nested <- sapply(nwl, function(nw1)
    blockName.vattr %in% list.vertex.attributes(nw1))

  if(!any(has_nested)){
    if(unlist) return(top_bname)
    return(as.list(top_bname))
  }

  if(NVL(bip, 0)){
    inner_parts <- lapply(nwl, function(nw1){
      if(blockName.vattr %in% list.vertex.attributes(nw1)){
        v <- get.vertex.attribute(nw1, blockName.vattr, unlist = FALSE)
        list(b1 = v[seq_len(b1.size(nw1))],
             b2 = v[b1.size(nw1) + seq_len(network.size(nw1) - b1.size(nw1))])
      } else {
        n1 <- network.size(nw1); e1 <- b1.size(nw1)
        list(b1 = as.list(rep(NA, e1)), b2 = as.list(rep(NA, n1 - e1)))
      }
    })
    inner_flat <- c(
      do.call(c, lapply(inner_parts, `[[`, "b1")),
      do.call(c, lapply(inner_parts, `[[`, "b2"))
    )
  } else {
    inner_flat <- do.call(c, lapply(nwl, function(nw1){
      if(blockName.vattr %in% list.vertex.attributes(nw1))
        get.vertex.attribute(nw1, blockName.vattr, unlist = FALSE)
      else
        as.list(rep(NA, network.size(nw1)))
    }))
  }

  result <- mapply(function(outer, inner)
    if(is.na(inner[[1]])) list(outer) else c(list(outer), inner),
    top_bname, inner_flat, SIMPLIFY = FALSE)

  if(unlist) sapply(result, `[`, 1) else result
}

#' @export
set.vertex.attribute.combined_networks <- function(x, attrname, value,
                                                    v = seq_len(network.size(x)),
                                                    ...) {
  # Distribute values across constituent networks using block offsets.
  # Not typically called on the combined network (attributes live in
  # the constituent networks), but provided for completeness.
  blockID.vattr   <- x$gal[[".blockID.vattr"]]
  blockName.vattr <- x$gal[[".blockName.vattr"]]
  if(identical(attrname, blockID.vattr) || identical(attrname, blockName.vattr))
    return(x)  # virtual attributes; ignore

  bip <- b1.size(x)
  nwl <- x$nw
  ns  <- sapply(nwl, network.size)

  if(NVL(bip, 0)){
    es    <- sapply(nwl, b1.size)
    eblks <- c(0, cumsum(es))
    ablks <- cumsum(c(bip, ns - es))
    n <- sum(ns)
    for(i in seq_along(nwl)){
      einds <- eblks[i] + seq_len(es[i])
      ainds <- ablks[i] + seq_len(ns[i] - es[i])
      combined_inds <- c(einds, ainds)
      nwl[[i]] <- set.vertex.attribute(nwl[[i]], attrname, value[combined_inds])
    }
  } else {
    blks <- c(0, cumsum(ns))
    for(i in seq_along(nwl)){
      inds <- blks[i] + seq_len(ns[i])
      nwl[[i]] <- set.vertex.attribute(nwl[[i]], attrname, value[inds])
    }
  }
  x$nw <- nwl
  x
}

#' @export
delete.vertex.attribute.combined_networks <- function(x, attrname, ...) {
  x$nw <- lapply(x$nw, function(nw1) {
    for(a in attrname)
      if(a %in% list.vertex.attributes(nw1))
        nw1 <- delete.vertex.attribute(nw1, a)
    nw1
  })
  x
}

# ---------------------------------------------------------------------------
# Edge methods
# ---------------------------------------------------------------------------

#' @export
network.edgecount.combined_networks <- function(x, na.omit = TRUE, ...) {
  sum(sapply(x$nw, network.edgecount, na.omit = na.omit))
}

#' @export
as.matrix.combined_networks <- function(x, matrix.type = "adjacency", ...) {
  # Build block-diagonal matrix from constituent networks.
  n   <- network.size(x)
  bip <- b1.size(x)
  nwl <- x$nw
  ns  <- sapply(nwl, network.size)

  if(NVL(bip, 0)){
    es    <- sapply(nwl, b1.size)
    eblks <- c(0, cumsum(es))
    ablks <- cumsum(c(bip, ns - es))
    m <- matrix(0, bip, n - bip)
    for(b in seq_along(nwl)){
      einds <- eblks[b] + seq_len(es[b])
      ainds <- ablks[b] + seq_len(ns[b] - es[b]) - bip
      m[einds, ainds] <- as.matrix(nwl[[b]], matrix.type = matrix.type, ...)
    }
  } else {
    blks <- c(0, cumsum(ns))
    m <- matrix(0, n, n)
    for(b in seq_along(nwl)){
      inds <- blks[b] + seq_len(ns[b])
      m[inds, inds] <- as.matrix(nwl[[b]], matrix.type = matrix.type, ...)
    }
  }
  m
}

# ---------------------------------------------------------------------------
# Printing / summary
# ---------------------------------------------------------------------------

#' @describeIn combine_networks A wrapper to print constituent network
#'   information.
#' @param x,object a combined network.
#' @param ... additional arguments to methods.
#' @export
print.combined_networks <- function(x, ...) {
  .print_combined_networks_info(x)
  invisible(x)
}

#' @describeIn combine_networks A summary method for combined networks.
#' @export
summary.combined_networks <- function(object, ...) {
  .print_combined_networks_info(object)
  invisible(object)
}
