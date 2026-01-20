#  File R/InitErgmConstraint.R in package ergm.multi, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################

##########################################################################################
# Each of the <InitErgmConstraint.X> functions accepts an existing constraint list, 'conlist',
# and to this adds an empty constraint list for term X; if any arguments are passed besides
# 'conlist", execution will halt.
#
# --PARAMETERS--
#   conlist: a list, presumably of constraints for other terms
#
# --RETURNED--
#   conlist: updated to include the initialized empty constraint list for term X
#
##########################################################################################

#' @templateVar name upper_tri
#' @title Only dyads in the upper-triangle of the sociomatrix may be
#'   toggled
#' @description For a directed network, only dyads \eqn{(i,j)} for
#'   which \eqn{i < j} may be toggled. Optional argument `attr`
#'   controls which subgraphs are thus restricted.
#'
#' @usage
#' # upper_tri(attr = NULL)
#' @template ergmTerm-attr
#'
#' @template ergmConstraint-general
#'
#' @concept dyad-independent
#' @concept directed
#' @import rle
InitErgmConstraint.upper_tri<-function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist, directed = TRUE,
                      varnames = c("attr"),
                      vartypes = c(ERGM_VATTR_SPEC),
                      defaultvalues = list(NULL),
                      required = c(FALSE))
  n <- network.size(nw)
  storage.mode(n) <- "integer"
  list(attr=a$attr,
       free_dyads = {
         restrict <- if(is.null(a$attr)) rep(TRUE, n) else ergm_get_vattr(a$attr, nw, accept = "logical")
         # The pattern is TRUE,...,TRUE,FALSE,...,FALSE for those
         # columns i where restrict[i]==TRUE, and it's just
         # TRUE,...,TRUE,TRUE,...,TRUE where restrict[i]==FALSE.
         d <- do.call(c, lapply(seq_len(n), function(i) rep(c(rle(TRUE),rle(!restrict[i])), c(i-1, n-i+1),scale="run")))
         rlebdm(compress(d), n)
       },
       dependence = FALSE
       )
}

#' @templateVar name blacklist_block
#' @title Blacklist blocks of dyads from toggling in a way that
#'   propagates through combined networks
#' @description This constraint is primarily for internal use,
#'   alongside the block blacklist API implemented in
#'   `block_blacklist.R`, which assigns and interprets vertex and
#'   network attributes allowing certain blocks of dyads to be
#'   blacklisted in a way that will survive even network manipulations
#'   such as subgraph extraction and combining into block-diagonal
#'   using [combine_networks()].
#'
#' @usage
#' # blacklist_block(block_vattr = ".ubid",
#' #                 blacklist_nattr = ".block_blacklist")
#' @param block_vattr a character string giving the name of vertex
#'   attribute containing globally unique vertex block membership IDs.
#' @param blacklist_nattr a character string giving the name of a
#'   network attribute containing the list of pairs of vertex IDs that
#'   are to be blacklisted.
#'
#' @template ergmConstraint-general
#'
#' @concept dyad-independent
#' @concept directed
#' @concept undirected
#' @keywords internal
InitErgmConstraint.blacklist_block<-function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("block_vattr", "blacklist_nattr"),
                      vartypes = c("character", "character"),
                      defaultvalues = list(".ubid", ".block_blacklist"),
                      required = c(FALSE, FALSE))
  n <- as.integer(network.size(nw))
  v <- get.vertex.attribute(nw, a$block_vattr, unlist=FALSE, na.omit=FALSE, null.na=FALSE)
  bl <- get_all_bl(nw, a$blacklist_nattr)

  list(
    block_vattr=a$block_vattr,
    blacklist_nattr=a$blacklist_nattr,
    free_dyads=if(length(bl)) (lapply(bl, function(b){
      btail <- in_block(v, b[[1]])
      bhead <- in_block(v, b[[2]])

      col0 <- rep(rle(FALSE),n)
      col1 <- rle(btail)

      # TODO: Optimize.
      do.call(c, lapply(seq_len(n), function(i) if(bhead[i]) col1 else col0)) %>% compress %>% `!`
    }) %>% reduce(`&`) %>% rlebdm(n)),
    dependence=FALSE
  )
}
