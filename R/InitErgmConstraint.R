#  File R/InitErgmConstraint.R in package ergm.multi, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2024 Statnet Commons
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

#### This constraint has been provisionally moved to 'ergm'. It may be
#### moved back here once 'ergm.multi' is on CRAN.
## #' @templateVar name blockdiag
## #' @title Block-diagonal structure constraint
## #' @description Force a block-diagonal structure (and its bipartite analogue) on
## #'   the network. Only dyads \eqn{(i,j)} for which
## #'   `attr(i)==attr(j)` can have edges.
## #'
## #'   Note that the current implementation requires that blocks be
## #'   contiguous for unipartite graphs, and for bipartite
## #'   graphs, they must be contiguous within a partition and must have
## #'   the same ordering in both partitions. (They do not, however,
## #'   require that all blocks be represented in both partitions, but
## #'   those that overlap must have the same order.)
## #'
## #'   If multiple block-diagonal constraints are given, or if
## #'   `attr` is a vector with multiple attribute names, blocks
## #'   will be constructed on all attributes matching.
## #'
## #' @usage
## #' # blockdiag(attr)
## #' @template ergmTerm-attr
## #'
## #' @template ergmConstraint-general
## #'
## #' @concept dyad-independent
## #' @concept directed
## #' @concept undirected
## #' @import rle
## InitErgmConstraint.blockdiag<-function(lhs.nw, attr=NULL, ...){
##   if(length(list(...)))
##     stop(paste("Block diagonal constraint takes one argument at this time."), call.=FALSE)
##   list(attr=attr,
##        free_dyads = {
##          n <- network.size(lhs.nw)
##          storage.mode(n) <- "integer"
##          a <- c(ergm_get_vattr(attr, lhs.nw)) # Strip attributes, which confuse rle().
##          if(NVL(lhs.nw%n%"bipartite",0)){
##            bip <- lhs.nw %n% "bipartite"
##            ea <- a[seq_len(bip)]
##            aa <- a[bip+seq_len(n-bip)]
##            if(length(rle(ea)$lengths)!=length(unique(rle(ea)$values)) || length(rle(aa)$lengths)!=length(unique(rle(aa)$values))) stop("Current implementation of block-diagonal sampling requires that the blocks of the egos and the alters be contiguous. See ", sQuote("ergmConstraint?blockdiag"), " for more information.")

##            tmp <- .double.rle(ea, aa)
##            el <- tmp$lengths1
##            al <- tmp$lengths2

##            o <- rlebdm(c(rep(rle(FALSE), bip*n, scale="run"),
##                          do.call(c,rep(
##                                      mapply(function(blen,bend){rep(rle(c(FALSE,TRUE,FALSE)), c(bend-blen, blen, n-bend), scale="run")},
##                                             el, cumsum(el), SIMPLIFY=FALSE),
##                                      al)
##                                  )), n)
##            # Future-proofing: in case it's bipartite directed, add
##            # both thte blocks and their transposes. (If undirected,
##            # it'll get filtered out by the .attributes constraints.)
##            ot <- rlebdm(c(do.call(c,rep(
##                                       mapply(function(blen,bend){rep(rle(c(FALSE,TRUE,FALSE)), c(bip+bend-blen, blen, n-bip-bend), scale="run")},
##                                              al, cumsum(al), SIMPLIFY=FALSE),
##                                       el)
##                                   ),
##                           rep(rle(FALSE), (n-bip)*n, scale="run")), n)
##            compress(o | ot)
##          }else{
##            a <- rle(a)
##            rlebdm(compress(do.call(c,rep(
##                                        mapply(function(blen,bend){rep(rle(c(FALSE,TRUE,FALSE)), c(bend-blen, blen, n-bend), scale="run")},
##                                               a$lengths, cumsum(a$lengths), SIMPLIFY=FALSE),
##                                        a$lengths)
##                                    )), n)
##          }
##        },
##        dependence = FALSE)
## }

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
