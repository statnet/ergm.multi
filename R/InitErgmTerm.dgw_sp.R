#  File R/InitErgmTerm.dgw_sp.R in package ergm.multi, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2024 Statnet Commons
################################################################################

.spcache.auxL <- function(type, Ls.path, L.in_order){
  type <- toupper(type)
  base_env(as.formula(as.call(list(as.name('~'), as.call(list(as.name('.spcache.netL'),type=if(type=='ITP')'OTP' else type,Ls.path=Ls.path,L.in_order=L.in_order))))))
}

.sp.handle_layers <- function(nw, a, type, has_base, cache.sp=FALSE){
  assert_LHS_Layer(nw)

  out <- list()

  namemap <- .layer_namemap(nw)

  if(is(a$Ls.path,"formula")) a$Ls.path <- list(a$Ls.path)
  if(length(a$Ls.path) == 1) a$Ls.path <- rep(a$Ls.path, 2)
  L.path1 <- NVL3(a$Ls.path[[1]], .set_layer_namemap(., namemap))
  L.path2 <- NVL3(a$Ls.path[[2]], .set_layer_namemap(., namemap))
  L.base <- NVL3(a$L.base, .set_layer_namemap(., namemap))

  if(is.null(L.path1) && is.null(L.path2) && is.null(L.base)) return(out)

  if(type=="RTP") stop("Layer-aware shared partner terms do not support reciprocated two-paths at this time.",call.=FALSE)

  NVL(L.path1) <- NVL(L.path2, L.base)
  NVL(L.path2) <- NVL(L.path1, L.base)
  if(has_base) NVL(L.base) <- NVL(L.path1, L.path2)
  
  layer0 <- .set_layer_namemap(
    if(has_base)
      as.formula(call("~",call("|",call("|", L.path1[[2]], L.path2[[2]]),L.base[[2]])))
    else
      as.formula(call("~",call("|", L.path1[[2]], L.path2[[2]]))),
    namemap)

  out$auxiliaries <- .mk_.layer.net_auxform(layer0)
  aux1 <- .mk_.layer.net_auxform(L.path1)
  out$auxiliaries[[2]] <- call("+", out$auxiliaries[[2]], aux1[[2]])
  aux2 <- .mk_.layer.net_auxform(L.path2)
  out$auxiliaries[[2]] <- call("+", out$auxiliaries[[2]], aux2[[2]])
  if(has_base){
    aux3 <- .mk_.layer.net_auxform(L.base)
    out$auxiliaries[[2]] <- call("+", out$auxiliaries[[2]], aux3[[2]])
  }
  if(cache.sp){
    aux4 <- .spcache.auxL(type, c(L.path1, L.path2), a$L.in_order)
    out$auxiliaries[[2]] <- call("+", out$auxiliaries[[2]], aux4[[2]])
  }
  
  out$any_order <- if(type=="UTP" || (type%in%c("OSP","ISP") && !has_base)) TRUE else !a$L.in_order
  out$coef.namewrap <- .lspec_coef.namewrap(list(pth=c(L.path1,if(L.path2!=L.path1)L.path2),bse=if(has_base) L.base,inord=a$L.in_order))
  out$name_suffix <- "_ML"
  out$nw1 <- subnetwork_templates(nw, ".LayerID", ".LayerName")[[1]] # Needed for emptynwstats.

  out
}

no_layer_err <- function(instead){
  ergm_Init_abort(paste("No layer specification found. Use", sQuote(paste0(instead, "()")), "instead."))
}

wrap_ergm_sp_call <- function(ergm_name, nw, a, has_base, d0 = FALSE, cache.sp = TRUE, ...) {
  # Construct and call the ergm term with only the arguments it
  # supports.
  not_ergm <- c("L.base", "Ls.path", "L.in_order")
  totransfer <- setdiff(names(a)[!attr(a, "missing")[names(a)]], not_ergm)
  cl <- as.call(c(ergm_name, a[totransfer]))
  trm <- call.ErgmTerm(cl, nw = nw, cache.sp = cache.sp, ...)

  # If null, return null.
  if (is.null(trm)) return(NULL)

  # Infer the type and construct layer info.
  type <- c("UTP", "OTP", "ITP", "RTP", "OSP", "ISP")[trm$iinputs[1] + 1]
  linfo <- .sp.handle_layers(nw, a, type, has_base, cache.sp)

  if(length(linfo)) nw <- linfo$nw1

  .emptynwstats <-
    if (d0 && any(a$d == 0)) {
      n <-
        if (is.bipartite(nw))
          switch(type,
                 UTP = network.size(nw),
                 OSP = b1.size(nw),
                 ISP = b2.size(nw))
        else network.size(nw)
      replace(dbl_along(a$d), a$d == 0,
              choose(n, 2L) * (is.directed(nw) + 1L))
    }

  # Replace the parts that are different for the layered term.
  if (length(linfo))
    within(trm, {
      name <- paste0(name, linfo$name_suffix)
      pkgname <- "ergm.multi"
      coef.names <- linfo$coef.namewrap(coef.names)
      auxiliaries <- linfo$auxiliaries
      iinputs <- c(linfo$any_order, iinputs)
      emptynwstats <- .emptynwstats
    })
  else no_layer_err(ergm_name)
}


################################################################################
#Term to count ESP statistics, where the shared partners may be any of
#several distinct types.
#
#Type codes are as follows (where (i,j) is the focal edge):
#
#  UTP - Undirected two-path (undirected graphs only)
#  OTP - Outgoing two-path (i->k->j)
#  ITP - Incoming two-path (i<-k<-j)
#  RTP - Reciprocated two-path (i<->k<->j)
#  OSP - Outgoing shared partner (i->k<-j)
#  ISP - Incoming shared partner (i<-k->j)
#
#Only one type may be specified per esp term.  The default, OTP, retains
#the original behavior of esp/gwesp.  In the undirected case, the UTP
#routine is used (since it is safe for undirected graphs), irrespective of
#the user's selection.  UTP cannot be chosen otherwise, since it won't work.
#

#' @templateVar name despL
#' @title Edgewise shared partners on layers
#' @description This term adds one network statistic to the model for each element in `d` where the \eqn{i} th such statistic equals the number of edges in the network with exactly `d[i]` shared partners. For a directed network, multiple shared partner definitions are possible.
#'
#' @usage
#' # binary: despL(d, type="OTP", L.base=NULL, Ls.path=NULL, L.in_order=FALSE)
#' @param d a vector of distinct integers
#' @template ergmTerm-sp-type
#' @template ergmTerm-L-base
#' @template ergmTerm-Ls-path
#'
#' @template ergmTerm-sp-types
#' @template ergmTerm-cache-sp
#' @template ergmTerm-general
#'
#' @concept directed
#' @concept undirected
#' @concept layer-aware
InitErgmTerm.despL<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("d","type","L.base","Ls.path","L.in_order"),
                      vartypes = c("numeric","character","formula","formula,list","logical"),
                      defaultvalues = list(NULL,"OTP",NULL,NULL,FALSE),
                      required = c(TRUE, FALSE,FALSE,FALSE,FALSE))

  wrap_ergm_sp_call("esp", nw, a, TRUE, ...)
}

#' @templateVar name despL
#' @template ergmTerm-rdname
#' @aliases espL-ergmTerm
#' @description `espL` and `despL` are aliases for consistency with \pkg{ergm}.
#' @usage
#' # binary: espL(d, type="OTP", L.base=NULL, Ls.path=NULL, L.in_order=FALSE)
InitErgmTerm.espL <- InitErgmTerm.despL

################################################################################
#Geometrically weighted edgewise shared partner term, where shared partners
#are defined by one of various specified types.
#
#Type codes are as follows (where (i,j) is the focal edge):
#
#  UTP - Undirected two-path (undirected graphs only)
#  OTP - Outgoing two-path (i->k->j)
#  ITP - Incoming two-path (i<-k<-j)
#  RTP - Reciprocated two-path (i<->k<->j)
#  OSP - Outgoing shared partner (i->k<-j)
#  ISP - Incoming shared partner (i<-k->j)
#
#Only one type may be specified per esp term.  The default, OTP, retains
#the original behavior of esp/gwesp.  In the undirected case, UTP is
#always used (since it is directedness-safe), and the user's input is
#overridden.  UTP cannot be chosen otherwise, since it won't work.
#

#' @templateVar name dgwespL
#' @title Geometrically weighted edgewise shared partner distribution on layers
#' @description This term adds a statistic equal to the geometrically weighted edgewise (not dyadwise) shared partner distribution with decay parameter. For a directed network, multiple shared partner definitions are possible.
#'
#' @usage
#' # binary: dgwespL(decay, fixed=FALSE, cutoff=30, type="OTP", L.base=NULL,
#' #                 Ls.path=NULL, L.in_order=FALSE)
#' @templateVar multiplicand shared partner or selected directed analogue count
#' @template ergmTerm-gw-decay-fixed
#' @templateVar underlying ESP
#' @template ergmTerm-gw-cutoff
#' @template ergmTerm-sp-type
#' @template ergmTerm-L-base
#' @template ergmTerm-Ls-path
#'
#' @template ergmTerm-sp-types
#' @template ergmTerm-cache-sp
#' @template ergmTerm-general
#'
#' @concept directed
#' @concept undirected
#' @concept layer-aware
InitErgmTerm.dgwespL<-function(nw, arglist, gw.cutoff=30, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("decay","fixed","cutoff","type", "alpha","L.base","Ls.path","L.in_order"),
                      vartypes = c("numeric","logical","numeric","character", "numeric","formula","formula,list","logical"),
                      defaultvalues = list(NULL, FALSE, gw.cutoff,"OTP", NULL,NULL,NULL,FALSE),
                      required = c(FALSE, FALSE, FALSE, FALSE, FALSE,FALSE,FALSE,FALSE))

  wrap_ergm_sp_call("gwesp", nw, a, TRUE, ...)
}

#' @templateVar name dgwespL
#' @template ergmTerm-rdname
#' @aliases gwespL-ergmTerm
#' @description `gdwespL` and `dgwespL` are aliases for consistency with \pkg{ergm}.
#' @usage
#' # binary: gwespL(decay, fixed=FALSE, cutoff=30, type="OTP", L.base=NULL,
#' #                Ls.path=NULL, L.in_order=FALSE)
InitErgmTerm.gwespL <- InitErgmTerm.dgwespL


#Term to count DSP statistics, where the shared partners may be any of
#several distinct types.
#
#Type codes are as follows (where (i,j) is the focal edge):
#
#  UTP - Undirected two-path (undirected graphs only)
#  OTP - Outgoing two-path (i->k->j)
#  ITP - Incoming two-path (i<-k<-j)
#  RTP - Reciprocated two-path (i<->k<->j)
#  OSP - Outgoing shared partner (i->k<-j)
#  ISP - Incoming shared partner (i<-k->j)
#
#Only one type may be specified per dsp term.  The default, OTP, retains
#the original behavior of dsp/gwdsp.  In the undirected case, the UTP
#routine is used (since it is safe for undirected graphs), irrespective of
#the user's selection.  UTP cannot be chosen otherwise, since it won't work.
#

#' @templateVar name ddspL
#' @title Dyadwise shared partners on layers
#' @description This term adds one network statistic to the model for each element in `d` where the \eqn{i} th such statistic equals the number of dyads in the network with exactly `d[i]` shared partners. For a directed network, multiple shared partner definitions are possible.
#'
#' @usage
#' # binary: ddspL(d, type="OTP", Ls.path=NULL, L.in_order=FALSE)
#'
#' @template ergmTerm-Ls-path
#' @template ergmTerm-sp-types
#' @template ergmTerm-cache-sp
#' @template ergmTerm-general
#'
#' @concept directed
#' @concept undirected
#' @concept layer-aware
InitErgmTerm.ddspL<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("d","type","Ls.path","L.in_order"),
                      vartypes = c("numeric","character","formula,list","logical"),
                      defaultvalues = list(NULL,"OTP",NULL,FALSE),
                      required = c(TRUE, FALSE,FALSE,FALSE))

  wrap_ergm_sp_call("dsp", nw, a, FALSE, TRUE, ...)
}

#' @templateVar name ddspL
#' @template ergmTerm-rdname
#' @aliases dspL-ergmTerm
#' @description `dspL` and `ddspL` are aliases for consistency with \pkg{ergm}.
#' @usage
#' # binary: dspL(d, type="OTP", Ls.path=NULL, L.in_order=FALSE)
InitErgmTerm.dspL <- InitErgmTerm.ddspL

################################################################################

#' @templateVar name dgwdspL
#' @title Geometrically weighted dyadwise shared partner distribution on layers
#' @description This term adds one network statistic to the model equal to the geometrically weighted dyadwise shared partner distribution with decay parameter. Note that the GWDSP statistic is equal to the sum of GWNSP plus GWESP. For a directed network, multiple shared partner definitions are possible.
#'
#' @usage
#' # binary: dgwdspL(decay, fixed=FALSE, cutoff=30, type="OTP",
#' #                 Ls.path=NULL, L.in_order=FALSE)
#' @templateVar multiplicand shared partner or selected directed analogue count
#' @template ergmTerm-gw-decay-fixed
#' @templateVar underlying DSP
#' @template ergmTerm-gw-cutoff
#' @template ergmTerm-sp-type
#' @template ergmTerm-Ls-path
#'
#' @template ergmTerm-sp-types
#' @template ergmTerm-cache-sp
#' @template ergmTerm-general
#'
#' @concept directed
#' @concept undirected
#' @concept layer-aware
InitErgmTerm.dgwdspL<-function(nw, arglist, cache.sp=TRUE, gw.cutoff=30, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("decay","fixed","cutoff","type", "alpha","Ls.path","L.in_order"),
                      vartypes = c("numeric","logical","numeric","character", "numeric","formula,list","logical"),
                      defaultvalues = list(NULL, FALSE, gw.cutoff,"OTP", NULL,NULL,FALSE),
                      required = c(FALSE, FALSE, FALSE, FALSE, FALSE,FALSE,FALSE))

  wrap_ergm_sp_call("gwdsp", nw, a, FALSE, ...)
}

#' @templateVar name dgwdspL
#' @template ergmTerm-rdname
#' @aliases gwdspL-ergmTerm 
#' @description `gdwdspL` and `dgwdspL` are aliases for consistency with \pkg{ergm}.
#' @usage
#' # binary: gwdspL(decay, fixed=FALSE, cutoff=30, type="OTP",
#' #                Ls.path=NULL, L.in_order=FALSE)
InitErgmTerm.gwdspL <- InitErgmTerm.dgwdspL

#Term to count NSP statistics, where the shared partners may be any of
#several distinct types.
#
#Type codes are as follows (where (i,j) is the focal edge):
#
#  UTP - Undirected two-path (undirected graphs only)
#  OTP - Outgoing two-path (i->k->j)
#  ITP - Incoming two-path (i<-k<-j)
#  RTP - Reciprocated two-path (i<->k<->j)
#  OSP - Outgoing shared partner (i->k<-j)
#  ISP - Incoming shared partner (i<-k->j)
#
#Only one type may be specified per NSP term.  The default, OTP, retains
#the original behavior of nsp/gwnsp.  In the undirected case, the UTP
#routine is used (since it is safe for undirected graphs), irrespective of
#the user's selection.  UTP cannot be chosen otherwise, since it won't work.
#

#' @templateVar name dnspL
#' @title Non-edgewise shared partners and paths on layers
#' @description This term adds one network statistic to the model for each element in `d` where the \eqn{i} th such statistic equals the number of non-edges in the network with exactly `d[i]` shared partners. For a directed network, multiple shared partner definitions are possible.
#'
#' @usage
#' # binary: dnspL(d, type="OTP", L.base=NULL, Ls.path=NULL, L.in_order=FALSE)
#' @param d a vector of distinct integers
#' @template ergmTerm-sp-type
#' @template ergmTerm-L-base
#' @template ergmTerm-Ls-path
#'
#' @template ergmTerm-sp-types
#' @template ergmTerm-cache-sp
#' @template ergmTerm-general
#'
#' @concept directed
#' @concept undirected
#' @concept layer-aware
InitErgmTerm.dnspL<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("d","type","L.base","Ls.path","L.in_order"),
                      vartypes = c("numeric","character","formula","formula,list","logical"),
                      defaultvalues = list(NULL,"OTP",NULL,NULL,FALSE),
                      required = c(TRUE, FALSE,FALSE,FALSE,FALSE))

  wrap_ergm_sp_call("nsp", nw, a, TRUE, TRUE, ...)
}

#' @templateVar name dnspL
#' @template ergmTerm-rdname
#' @aliases nspL-ergmTerm
#' @description `nspL` and `dnspL` are aliases for consistency with \pkg{ergm}.
#' @usage
#' # binary: nspL(d, type="OTP", L.base=NULL, Ls.path=NULL, L.in_order=FALSE)
InitErgmTerm.nspL <- InitErgmTerm.dnspL

################################################################################

#' @templateVar name dgwnspL
#' @title Geometrically weighted non-edgewise shared partner distribution on layers
#' @description This term is just like [`gwespL`][gwespL-ergmTerm] and [`gwdspL`][gwdspL-ergmTerm] except it adds a statistic equal to the geometrically weighted nonedgewise (that is, over dyads that do not have an edge) shared partner distribution with decay parameter. For a directed network, multiple shared partner definitions are possible.
#'
#' @usage
#' # binary: dgwnspL(decay, fixed=FALSE, cutoff=30, type="OTP", L.base=NULL,
#' #                 Ls.path=NULL, L.in_order=FALSE)
#' @templateVar multiplicand shared partner or selected directed analogue count
#' @template ergmTerm-gw-decay-fixed
#' @templateVar underlying NSP
#' @template ergmTerm-gw-cutoff
#' @template ergmTerm-sp-type
#' @template ergmTerm-L-base
#' @template ergmTerm-Ls-path
#'
#' @template ergmTerm-sp-types
#' @template ergmTerm-cache-sp
#' @template ergmTerm-general
#'
#' @concept directed
#' @concept undirected
#' @concept layer-aware
InitErgmTerm.dgwnspL<-function(nw, arglist, gw.cutoff=30, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("decay","fixed","cutoff","type", "alpha","L.base","Ls.path","L.in_order"),
                      vartypes = c("numeric","logical","numeric","character", "numeric","formula","formula,list","logical"),
                      defaultvalues = list(NULL, FALSE, gw.cutoff,"OTP", NULL,NULL,NULL,FALSE),
                      required = c(FALSE, FALSE, FALSE, FALSE, FALSE,FALSE,FALSE,FALSE))

  wrap_ergm_sp_call("gwnsp", nw, a, TRUE, ...)
}

#' @templateVar name dgwnspL
#' @template ergmTerm-rdname
#' @aliases gwnspL-ergmTerm
#' @description `gdwnspL` and `dgwnspL` are aliases for consistency with \pkg{ergm}.
#' @usage
#' # binary: gwnspL(decay, fixed=FALSE, cutoff=30, type="OTP", L.base=NULL,
#' #                Ls.path=NULL, L.in_order=FALSE)
InitErgmTerm.gwnspL <- InitErgmTerm.dgwnspL


################################################################################

#' @templateVar name b1dspL
#' @title Dyadwise shared partners for dyads in the first bipartition on layers
#' @description This term adds one
#'   network statistic to the model for each element in `d`; the \eqn{i}th
#'   such statistic equals the number of dyads in the first bipartition with exactly
#'   `d[i]` shared partners. (Those shared partners, of course, must be members
#'   of the second bipartition.) This term can only be used with bipartite networks.
#'
#' @usage
#' # binary: b1dsp(d, Ls.path=NULL)
#'
#' @param d a vector of distinct integers.
#' @template ergmTerm-Ls-path
#'
#' @template ergmTerm-cache-sp
#' @template ergmTerm-general
#'
#' @concept bipartite
#' @concept undirected
#' @concept layer-aware
InitErgmTerm.b1dspL <- function(nw, arglist, cache.sp=TRUE, ...){
  a <- check.ErgmTerm(nw, arglist, bipartite = TRUE,
                      varnames = c("d", "Ls.path"),
                      vartypes = c("numeric", "formula,list"),
                      defaultvalues = list(NULL, NULL),
                      required = c(TRUE, FALSE))

  wrap_ergm_sp_call("b1dsp", nw, a, FALSE, TRUE, ...)
}

################################################################################

#' @templateVar name gwb1dsp
#' @title Geometrically weighted dyadwise shared partner distribution for dyads in the first bipartition on layers
#' @description This term adds one network statistic to the model equal to the geometrically
#'   weighted dyadwise shared partner distribution for dyads in the first bipartition with decay parameter
#'   `decay` parameter, which should be non-negative. This term can only be used with bipartite networks.
#'
#' @usage
#' # binary: gwb1dsp(decay=0, fixed=FALSE, cutoff=30, Ls.path=NULL)
#'
#' @templateVar multiplicand shared partner counts
#' @template ergmTerm-gw-decay-fixed
#' @templateVar underlying b1dsp
#' @template ergmTerm-gw-cutoff
#' @template ergmTerm-Ls-path
#'
#' @template ergmTerm-cache-sp
#' @template ergmTerm-general
#'
#' @concept bipartite
#' @concept undirected
#' @concept curved
#' @concept layer-aware
InitErgmTerm.gwb1dspL<-function(nw, arglist, cache.sp=TRUE, gw.cutoff=30, ...) {
  a <- check.ErgmTerm(nw, arglist, bipartite = TRUE,
                      varnames = c("decay", "fixed", "cutoff", "alpha", "Ls.path"),
                      vartypes = c("numeric", "logical", "numeric", "numeric", "formula,list"),
                      defaultvalues = list(NULL, FALSE, gw.cutoff, NULL, NULL),
                      required = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE))

  wrap_ergm_sp_call("gwb1dsp", nw, a, FALSE, ...)
}

################################################################################

#' @templateVar name b2dspL
#' @title Dyadwise shared partners for dyads in the second bipartition on layers
#' @description This term adds one network statistic to the model for each element in `d` ; the \eqn{i} th
#'   such statistic equals the number of dyads in the second bipartition with exactly
#'   `d[i]` shared partners. (Those shared partners, of course, must be members
#'   of the first bipartition.) This term can only be used with bipartite networks.
#'
#' @usage
#' # binary: b2dsp(d, Ls.path=NULL)
#'
#' @param d a vector of distinct integers
#' @template ergmTerm-Ls-path
#'
#' @template ergmTerm-cache-sp
#' @template ergmTerm-general
#'
#' @concept bipartite
#' @concept undirected
#' @concept layer-aware
InitErgmTerm.b2dspL <- function(nw, arglist, cache.sp=TRUE, ...){
  a <- check.ErgmTerm(nw, arglist, bipartite = TRUE,
                      varnames = c("d", "Ls.path"),
                      vartypes = c("numeric", "formula,list"),
                      defaultvalues = list(NULL, NULL),
                      required = c(TRUE, FALSE))

  wrap_ergm_sp_call("b2dsp", nw, a, FALSE, TRUE, ...)
}

################################################################################

#' @templateVar name gwb2dspL
#' @title Geometrically weighted dyadwise shared partner distribution for dyads in the second bipartition on layers
#' @description This term adds one network statistic to the model equal to the geometrically
#'   weighted dyadwise shared partner distribution for dyads in the second bipartition with decay parameter
#'   `decay` parameter, which should be non-negative. This term can only be used with bipartite networks.
#'
#' @usage
#' # binary: gwb2dsp(decay=0, fixed=FALSE, cutoff=30, Ls.path=NULL)
#'
#' @templateVar multiplicand shared partner counts
#' @template ergmTerm-gw-decay-fixed
#' @templateVar underlying b2dsp
#' @template ergmTerm-gw-cutoff
#' @template ergmTerm-Ls-path
#'
#' @template ergmTerm-cache-sp
#' @template ergmTerm-general
#'
#' @concept bipartite
#' @concept undirected
#' @concept curved
#' @concept layer-aware
InitErgmTerm.gwb2dsp<-function(nw, arglist, cache.sp=TRUE, gw.cutoff=30, ...) {
  a <- check.ErgmTerm(nw, arglist, bipartite = TRUE,
                      varnames = c("decay", "fixed", "cutoff", "alpha", "Ls.path"),
                      vartypes = c("numeric", "logical", "numeric", "numeric", "formula,list"),
                      defaultvalues = list(NULL, FALSE, gw.cutoff, NULL, NULL),
                      required = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE))

  wrap_ergm_sp_call("gwb2dsp", nw, a, FALSE, ...)
}
