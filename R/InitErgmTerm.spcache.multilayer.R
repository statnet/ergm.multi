#  File R/InitErgmTerm.spcache.multilayer.R in package ergm.multi, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################

InitErgmTerm..spcache.netL<-function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("type", "Ls.path", "L.in_order"),
                      vartypes = c("character", "formula,list", "logical"),
                      defaultvalues = list(NULL,NULL,NULL),
                      required = c(TRUE, TRUE, TRUE))

  assert_LHS_Layer(nw)

  type <- match.arg(tolower(a$type), c("otp","osp","isp","utp")) # ITP not included, because it's just OTP with direction reversed.

  if (is.directed(nw) == (type == "utp")
      && !(is.bipartite(nw) && type %in% c("osp", "isp")))
    stop("Type UTP may only be used with undirected networks, OSP and ISP with bipartite or directed, and the rest only with directed.")

  dname <- paste0("_",type,"_wtnet")
  linfo <- .sp.handle_layers(nw, a, type, FALSE)

  
  list(name=paste0(dname,linfo$name_suffix), auxiliaries=linfo$auxiliaries, iinputs=c(linfo$any_order),
       coef.names=c(), dependence=TRUE)
}
