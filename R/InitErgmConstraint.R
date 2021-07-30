#  File R/InitErgmConstraint.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################
#============================================================================
# This file contains the following 12 functions for initializing empty
# constraint lists (each prependend with "InitErgmConstraint")
#         <edges>                   <odegreedist>
#         <degrees>=<nodedegrees>   <bd>
#         <degreesTetrad>           <idegrees>
#         <degreesHexad>            <odegrees>
#         <degreedist>              <hamming>
#         <idegreedist>            <observed>
#============================================================================

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

#' @name blockdiag-ergmConstraint
#' @title Block-diagonal structure constraint
#' @description Block-diagonal structure constraint
#' @details Force a block-diagonal structure (and its bipartite analogue) on
#'   the network. Only dyads \eqn{(i,j)} for which
#'   `attrname(i)==attrname(j)` can have edges.
#'
#'   Note that the current implementation requires that blocks be
#'   contiguous for unipartite graphs, and for bipartite
#'   graphs, they must be contiguous within a partition and must have
#'   the same ordering in both partitions. (They do not, however,
#'   require that all blocks be represented in both partitions, but
#'   those that overlap must have the same order.)
#'
#'   If multiple block-diagonal constraints are given, or if
#'   `attrname` is a vector with multiple attribute names, blocks
#'   will be constructed on all attributes matching.
#'
#' @usage
#' # blockdiag(attrname)
#' @template ergmTerm-attr
#'
#' @template ergmConstraint-general
#'
#' @import rle
#'
#' @concept dyad-independent
InitErgmConstraint.blockdiag<-function(lhs.nw, attr=NULL, ...){
  if(length(list(...)))
    stop(paste("Block diagonal constraint takes one argument at this time."), call.=FALSE)
  list(attr=attr,
       free_dyads = {
         n <- network.size(lhs.nw)
         storage.mode(n) <- "integer"
         a <- c(ergm_get_vattr(attr, lhs.nw)) # Strip attributes, which confuse rle().
         if(NVL(lhs.nw%n%"bipartite",0)){
           bip <- lhs.nw %n% "bipartite"
           ea <- a[seq_len(bip)]
           aa <- a[bip+seq_len(n-bip)]
           if(length(rle(ea)$lengths)!=length(unique(rle(ea)$values)) || length(rle(aa)$lengths)!=length(unique(rle(aa)$values))) stop("Current implementation of block-diagonal sampling requires that the blocks of the egos and the alters be contiguous. See help('ergm-constraints') for more information.")
           
           tmp <- .double.rle(ea, aa)
           el <- tmp$lengths1
           al <- tmp$lengths2
           
           o <- rlebdm(c(rep(rle(FALSE), bip*n, scale="run"),
                         do.call(c,rep(
                                     mapply(function(blen,bend){rep(rle(c(FALSE,TRUE,FALSE)), c(bend-blen, blen, n-bend), scale="run")},
                                            el, cumsum(el), SIMPLIFY=FALSE),
                                     al)
                                 )), n)
           # Future-proofing: in case it's bipartite directed, add
           # both thte blocks and their transposes. (If undirected,
           # it'll get filtered out by the .attributes constraints.)
           ot <- rlebdm(c(do.call(c,rep(
                                      mapply(function(blen,bend){rep(rle(c(FALSE,TRUE,FALSE)), c(bip+bend-blen, blen, n-bip-bend), scale="run")},
                                             al, cumsum(al), SIMPLIFY=FALSE),
                                      el)
                                  ),
                          rep(rle(FALSE), (n-bip)*n, scale="run")), n)
           compress(o | ot)
         }else{
           a <- rle(a)
           rlebdm(compress(do.call(c,rep(
                                       mapply(function(blen,bend){rep(rle(c(FALSE,TRUE,FALSE)), c(bend-blen, blen, n-bend), scale="run")},
                                              a$lengths, cumsum(a$lengths), SIMPLIFY=FALSE),
                                       a$lengths)
                                   )), n)
         }
       },
       dependence = FALSE)
}


InitErgmConstraint.upper_tri<-function(lhs.nw, attrname=NULL, ...){
  if(length(list(...)))
    stop(paste("Upper-triangular constraint takes one argument at this time."), call.=FALSE)
  if(!is.directed(lhs.nw)) stop("Upper-triangular constraint can be used only with directed networks.")
  n <- network.size(lhs.nw)
  storage.mode(n) <- "integer"
  list(attrname=attrname,
       free_dyads = {
         restrict <- if(is.null(attrname)) rep(TRUE, n) else as.logical(lhs.nw %v% attrname)
         # The pattern is TRUE,...,TRUE,FALSE,...,FALSE for those
         # columns i where restrict[i]==TRUE, and it's just
         # TRUE,...,TRUE,TRUE,...,TRUE where restrict[i]==FALSE.
         d <- do.call(c, lapply(seq_len(n), function(i) rep(c(rle(TRUE),rle(!restrict[i])), c(i-1, n-i+1),scale="run")))
         rlebdm(compress(d), n)
       },
       dependence = FALSE
       )
}

InitErgmConstraint.blacklist_block<-function(lhs.nw, block_vattr=".ubid", blacklist_nattr=".block_blacklist"){
  n <- as.integer(network.size(lhs.nw))
  v <- get.vertex.attribute(lhs.nw,block_vattr,unlist=FALSE,na.omit=FALSE,null.na=FALSE)
  bl <- get_all_bl(lhs.nw, blacklist_nattr)

  list(
    block_vattr=block_vattr,
    blacklist_nattr=blacklist_nattr,
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
