#  File R/block_blacklist.R in package ergm.multi, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2022 Statnet Commons
################################################################################
#' Generate a list of (at least somewhat) globally unique identifiers
#' (vectors of integers).
#' @param n number of identifiers to generate.
#' @param minbytes minimum length of an identifier in bytes; will be
#'   rounded up to the nearest multiple of 4.
#' @param runif_only if `FALSE`, also mix the time stamp and the
#'   process ID into the value, so that even the same random seed will
#'   (probably) produce a different result.
#' @return
#' A list of lists of integers.
#' @noRd
guid <- function(n, minbytes=8, runif_only=FALSE){
  ints <- ceiling(minbytes/4)
  # R really needs an API for generating 4 random bytes.

  x <- runif(ints*n)
  if(!runif_only){
    # After generating the random numbers, also shift them by PID and current time so that they are unique even if the user uses set.seed().
    x <- (x + (Sys.getpid()/.Machine$integer.max)) %% 1
    x <- (x + (as.integer(Sys.time())/.Machine$integer.max)) %% 1
  }
  x <- as.integer((x-0.5)*2*.Machine$integer.max)
  split(x, rep(seq_len(n),each=ints))
}

#' Add or append to a vertex attribute a set of globally unique
#' identifiers.
#' @param nw a [`network`] object.
#' @param blocks a factor vector equal in length to the network size of `nw`
#'   whose values identify the (not necessarily contiguous) blocks.
#' @param attrname vertex attribute name to use.
#' @return A `list()` with two elements: \describe{
#' 
#' \item{nw}{A network with `attrname` augmented with the new IDs.}
#' \item{id}{A list of IDs generated in the same order as `levels(factor(blocks))`.}
#' }
#' @noRd
add_block_guid <- function(nw, blocks, attrname=".ubid"){
  if(!is.factor(blocks)){
    warning("Recommend specifying blocks as a factor, so that there would be no ambiguity in levels.")
    blocks <- factor(blocks)
  }
  id <- unname(guid(nlevels(blocks)))
  new_id <- id[as.integer(blocks)]
  old_id <- get.vertex.attribute(nw, attrname, null.na=FALSE, na.omit=FALSE, unlist=FALSE)
  list(nw=set.vertex.attribute(nw, attrname, mapply(function(new,old)c(list(new),old,recursive=FALSE), new_id, old_id, SIMPLIFY=FALSE, USE.NAMES=FALSE)),
       id=id)
}

#' Identify which vertices have a particular globally unique ID (among
#' others).
#' @param v a list of lists of `guid`s.
#' @param guid a one or more `guid`, either as an integer vector or a
#'   list containing one or more integer vectors.
#' @return If `guid` is an integer vector or a list whose sole element
#'   is such a vector, returns a logical vector of length
#'   `network.size(nw)` indicating whether each vertex carries that
#'   ID. If `guid` is a list of several integer vectors, returns a
#'   corresponding matrix, with vertices in columns and corresponding
#'   `guid` presences in rows.
#' @noRd
in_block <- function(v, guid){
  if(!is.list(guid)) guid <- list(guid)
  sapply(v, `%in%`, x=guid)
}

#' Blacklist some blocks from toggling based on where certain vertex
#' values are TRUE.
#' @param nw a [`network`].
#' @param tails,heads logical or numeric vectors selecting which tails
#'   or heads of dyads should be blacklisted or whitelisted.
#' @param block_vattr,blacklist_nattr vertex and network attributes
#'   into which to put block IDs and the blacklists.
#' @param invert whether everything *but* the specified intersection
#'   should be blacklisted.
#' @return A modified [`network`] object with `block_vattr` augmented
#'   with `guid`s identifying the block they belong to, and
#'   `blacklist_nattr` network attribute containing a list of `guid`
#'   pairs identifying the blacklisted blocks.
#' @noRd
blacklist_intersect <- function(nw, tails, heads=tails, block_vattr=".ubid", blacklist_nattr=".block_blacklist", invert=FALSE){
  if(!is.logical(tails)) tails <- unwhich(tails, network.size(nw))
  if(!is.logical(heads)) heads <- unwhich(heads, network.size(nw))
  vtype <- list(tails,heads) %>% transpose %>% match(list(list(FALSE,FALSE),list(FALSE,TRUE),list(TRUE,FALSE),list(TRUE,TRUE)))
  uvtype <- sort(unique(vtype))
  tmp <- add_block_guid(nw, factor(vtype,levels=1:4), attrname=block_vattr)
  nw <- tmp$nw
  id <- tmp$id
  rm(tmp)

  l <- expand.grid(tail=c(3L,4L),head=c(2L,4L)) %>% unname %>% transpose
  if(invert){
    al <- expand.grid(tail=uvtype,head=uvtype) %>% unname %>% transpose
    l <- setdiff(al,l)
  }
  l <- lapply(l, lapply, function(i)id[[i]]) # Remap vtypes to guids.
  nw %n% blacklist_nattr <- c(nw %n% blacklist_nattr, l, recursive=FALSE)
  nw
}

flatten_guid_pairs <- function(guidll){
  if(length(guidll))
    guidll %>% unlist %>% unname %>%
      split(.,(seq_along(.)-1L)%/%4L) %>% unname %>%
      lapply(split, c(1L,1L,2L,2L)) %>% lapply(unname)
}

get_all_bl_subnetattr <- function(l){
  lapply(l, function(subl)
    c(subl$.block_blacklist, lapply(subl$.subnetattr, get_all_bl_subnetattr), recursive=FALSE)
    ) %>% flatten_guid_pairs
}

#' Obtain and concatenate unique block blacklists from the specified
#' networks, which may be combined networks.
#' @noRd
get_all_bl <- function(nw, blacklist_nattr=".block_blacklist"){
  c(nw%n%blacklist_nattr, get_all_bl_subnetattr(nw%n%".subnetattr")) %>% unique
}

