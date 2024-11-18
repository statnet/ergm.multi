#  File R/InitErgmTerm.dgw_sp.R in package ergm.multi, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2024 Statnet Commons
################################################################################

#  ------------------------------------------------------------------ 
#   Description of the input and output parameters of the  
#   InitErgmTerm.xxx function, where xxx is the name of your term
#  ------------------------------------------------------------------ 
#
#  INPUTS:
#  Each InitErgmTerm function takes three arguments:
#	  		nw: The network of interest
#      arglist: The list of arguments passed to the term xxx
#         ... : There may be other arguments passed by 
#               ergm_model, so each InitErgmTerm function 
#               must include the ... argument
#  These inputs are automatically supplied by ergm_model.
#
#  OUTPUTS:
#  Each InitErgmTerm function should return a list.  
#     REQUIRED LIST ITEMS:
#          name: This names the C changestats function for term xxx, 
#                but does so by excluding the d_ prefix. The 
#                changestats function is named d_xxxy and 'name' is
#                consequently "xxxy". For example, the b1starmix
#                term has 2 changestats functions based on
#                whether the homophily argument is set. These are
#                d_b1starmix and d_b1starmixhomophily. The 'name' 
#                returned by InitErgmTerm.b1starmix is then one of 
#                "b1starmix" or "b1starmixhomophily" as appropriate.
#    coef.names: Vector of names for the coefficients (parameters)
#                as they will be reported in the output.
#       pkgname: This names the package containing the C changestats
#                function d_[name]. The default is "ergm", which means
#                that if you have code that exists as part of the 
#                (say) "ergm.userterms" package, you MUST specify 
#                pkgname="ergm.userterms"
#
#    OPTIONAL LIST ITEMS:
#        inputs: Vector of (double-precision numeric) inputs that the 
#                changestat function called d_'name' may require.
#                The default is NULL; no inputs are required.  But it
#                MUST be a vector!  Thus, if some of the inputs are,  
#                say, matrices, they must be "flattened" to vectors; if 
#                some are categorical character-valued variables, they
#                must be converted to numbers. Optionally, the inputs 
#                vector may have an attribute named "ParamsBeforeCov",
#                which is the number of input parameters preceding the 
#                covariate vector in 'inputs'.  This is necessary for 
#                compatibility with some of the existing d_xxx changestats 
#                functions in ergm, but is not necessary in general.
#    dependence: Logical variable telling whether addition of this term to
#                the model makes the model into a dyadic dependence model.
#                If none of the terms sets dependence==TRUE, then the model
#                is assumed to be a dyadic independence model, which means
#                that the pseudolikelihood estimate coincides with the
#                maximum likelihood estimate.  The default value is TRUE.
#  emptynwstats: Vector of values (if nonzero) for the statistics evaluated
#                on the empty network.  If all are zero for this term, this
#                argument may be omitted.  For example, the degree0 term 
#                would require 'emptynwstats' since degree0 = number of 
#                nodes for the empty network.
#        params: For curved exponential family model terms only, a list of 
#                (numeric) initial values for the parameters of  
#                curved exponential family model terms. Each item in the  
#                list should be named with the corresponding parameter name 
#                (one or more of these will probably coincide with the 
#                 coef.names).  For example, the gwesp term returns 
#                params=list(gwesp=NULL,gwesp.decay=decay), where decay
#                was specified as an argument to the gwesp term. 
#           map: For curved exponential family model terms only, a function 
#                giving the map from the canonical parameters, theta,
#                associated with the statistics for this term, to eta, 
#                the corresponding curved parameters.  The length of eta 
#                is the same as the length of the 'params' list above.
#                The function takes two arguments:  theta and length(eta).
#      gradient: For curved exponential family model terms only, a function 
#                giving the gradient of the 'map'. If theta has length p 
#                and eta has length q, then gradient should return a
#                p by q matrix. This function takes two arguments:  theta 
#                and length(eta).
#


#  ------------------------------------------------------------------------- 
#   Description of the input parameters to the d_xxxy changestats function, 
#   where xxxy corresponds to the 'name' returned by InitErgmTerm.xxx.
#  -------------------------------------------------------------------------- 
#
#  INPUTS:
#  Each d_xxxy function takes five arguments:
#	    ntoggles: the number of toggles as described in 
#                 "ergm.userterms: A template package"
#          heads: a pointer to the array of the head nodes of the 
#                 proposed edges to be toggled
#          tails: a pointer to the array of the tail nodes of the
#                 proposed edges to be toggled
#            mtp: a pointer to the model, which includes the following:
#                 dstats      : a pointer to the array of changestats,
#                               macro-ed as CHANGE_STAT
#                 nstats      : the length of 'dstats', macro-ed as
#                               N_CHANGE_STATS
#                 inputparams : a pointer to the vector of input 
#                               parameters. This is supplied by the
#                               'inputs' returned by InitErgmTerm.xxx
#                               and is macro-ed as INPUT_PARAM
#                 ninputparams: the length of 'inputparams', macro-ed
#                               as N_INPUT_PARAMS
#            nwp: a pointer to the network.  This includes several 
#                 components and several macros exist for accessing
#                 these. See the changestat.h file for a list of these
#                 components and their macros. 
#  These inputs are automatically supplied to the d_xxxy function by the 
#  network_stats_wrapper function 

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
  ergm_Init_abort(paste("No layer specification found. Use", sQuote(instead), "instead."))
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
InitErgmTerm.despL<-function(nw, arglist, cache.sp=TRUE, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("d","type","L.base","Ls.path","L.in_order"),
                      vartypes = c("numeric","character","formula","formula,list","logical"),
                      defaultvalues = list(NULL,"OTP",NULL,NULL,FALSE),
                      required = c(TRUE, FALSE,FALSE,FALSE,FALSE))
  d<-a$d
  ld<-length(d)
  if(ld==0){return(NULL)}
  
  type<-toupper(a$type[1])
  type.vec<-c("OTP","ITP","RTP","OSP","ISP")
  if(!(type%in%type.vec))
    stop("Illegal type code for esp; valid types are:",paste(type.vec, collapse=","))
  dname<-"esp"
  if(is.directed(nw)){
    conam <- paste("esp",type,sep=".")
    typecode<-which(type==type.vec)
    dname <- "desp"
  }else{
    dname <- "desp"
    conam<-"esp"
    type<-"UTP"
    typecode<-0
  }

  linfo <- .sp.handle_layers(nw, a, type, TRUE, cache.sp)
  
  if(length(linfo)) list(name=paste0(dname,linfo$name_suffix), coef.names=linfo$coef.namewrap(paste(conam,d,sep="")), auxiliaries=linfo$auxiliaries, iinputs=c(linfo$any_order,typecode,d), minval=0)
  else no_layer_err("desp()")
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
InitErgmTerm.dgwespL<-function(nw, arglist, cache.sp=TRUE, gw.cutoff=30, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("decay","fixed","cutoff","type", "alpha","L.base","Ls.path","L.in_order"),
                      vartypes = c("numeric","logical","numeric","character", "numeric","formula","formula,list","logical"),
                      defaultvalues = list(NULL, FALSE, gw.cutoff,"OTP", NULL,NULL,NULL,FALSE),
                      required = c(FALSE, FALSE, FALSE, FALSE, FALSE,FALSE,FALSE,FALSE))
  if(!is.null(a$alpha)){
    stop("For consistency with gw*degree terms, in all gw*sp and dgw*sp terms the argument ", sQuote("alpha"), " has been renamed to " ,sQuote("decay"), ".", call.=FALSE)
  }
  
  decay<-a$decay;fixed<-a$fixed
  cutoff<-a$cutoff
  decay=decay[1] # Not sure why anyone would enter a vector here, but...
  type<-toupper(a$type[1])
  type.vec<-c("OTP","ITP","RTP","OSP","ISP")
  if(!(type%in%type.vec))
    stop("Illegal type code for gwesp; valid types are:",paste(type.vec, collapse=","))
  dname<-"desp"
  if(!is.directed(nw)){  
    type <- "UTP"
    typecode<-0
    basenam<-paste("gwesp",sep=".")
  }else{
    typecode<-which(type==type.vec)
    basenam<-paste("gwesp",type,sep=".")
  }
  
  linfo <- .sp.handle_layers(nw, a, type, TRUE, cache.sp)
  
  if(!fixed){ # This is a curved exponential family model
    if(!is.null(a$decay)) warning("In term 'dgwesp': decay parameter 'decay' passed with 'fixed=FALSE'. 'decay' will be ignored. To specify an initial value for 'decay', use the 'init' control parameter.", call.=FALSE)

    maxesp <- min(cutoff,network.size(nw)-2)
    d <- 1:maxesp
    ld<-length(d)
    if(ld==0){return(NULL)}
    params<-list(gwesp=NULL,gwesp.decay=decay)
    names(params)<-c(basenam,paste(basenam,"decay",sep="."))

    if(length(linfo)) c(list(name=paste0(dname,linfo$name_suffix),
                             coef.names=linfo$coef.namewrap(if(is.directed(nw)) paste("esp.",type,"#",d,sep="") else paste("esp#",d,sep="")),auxiliaries=linfo$auxiliaries,
                             iinputs=c(linfo$any_order,typecode,d), params=params), GWDECAY)
    else no_layer_err("dgwesp()")
  }else{
    if(is.null(a$decay)) stop("Term 'dgwesp' with 'fixed=TRUE' requires a decay parameter 'decay'.", call.=FALSE)

    dname<-"dgwesp"
    maxesp <- min(cutoff,network.size(nw)-2)
    if(is.directed(nw))
      coef.names <- paste(paste("gwesp",type,"fixed.",sep="."),decay, sep="")
    else
      coef.names <- paste("gwesp.fixed.",decay,sep="")

    if(length(linfo)) list(name=paste0(dname,linfo$name_suffix), coef.names=linfo$coef.namewrap(coef.names), inputs=decay, iinputs=c(linfo$any_order,typecode,maxesp), auxiliaries=linfo$auxiliaries)
    else no_layer_err("dgwesp()")
  }
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
InitErgmTerm.ddspL<-function(nw, arglist, cache.sp=TRUE, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("d","type","Ls.path","L.in_order"),
                      vartypes = c("numeric","character","formula,list","logical"),
                      defaultvalues = list(NULL,"OTP",NULL,FALSE),
                      required = c(TRUE, FALSE,FALSE,FALSE))
  d<-a$d
  ld<-length(d)
  if(ld==0){return(NULL)}
  
  type<-toupper(a$type[1])
  type.vec<-c("OTP","ITP","RTP","OSP","ISP")
  if(!(type%in%type.vec))
    stop("Illegal type code for sp; valid types are:",paste(type.vec, collapse=","))
  dname<-"ddsp"
  if(is.directed(nw)){
    conam <- paste("dsp",type,sep=".")
    typecode<-which(type==type.vec)
    dname <- "ddsp"
  }else{
    conam <- paste("dsp",sep=".")
    type<-"UTP"
    typecode<-0
  }

  linfo <- .sp.handle_layers(nw, a, type, FALSE, cache.sp)
  if(length(linfo)) nw <- linfo$nw1
  
  if (any(d==0)) {
    emptynwstats <- rep(0, length(d))
    if(is.bipartite(nw)){
      nb1 <- get.network.attribute(nw, "bipartite")
      nb2 <- network.size(nw) - nb1
      emptynwstats[d==0] <- nb1*(nb1-1)/2 + nb2*(nb2-1)/2
    }else{
      emptynwstats[d==0] <- network.dyadcount(nw,FALSE)
    }
  }else{
    emptynwstats <- NULL
  }
  
  if(length(linfo)) list(name=paste0(dname,linfo$name_suffix), coef.names=linfo$coef.namewrap(paste0(conam,d)), auxiliaries=linfo$auxiliaries, iinputs=c(linfo$any_order,typecode,d), minval=0, emptynwstats=emptynwstats)
  else no_layer_err("ddsp()")
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
  if(!is.null(a$alpha)){
    stop("For consistency with gw*degree terms, in all gw*sp and dgw*sp terms the argument ", sQuote("alpha"), " has been renamed to " ,sQuote("decay"), ".", call.=FALSE)
  }
  
  decay<-a$decay;fixed<-a$fixed
  cutoff<-a$cutoff
  decay=decay[1] # Not sure why anyone would enter a vector here, but...
  
  type<-toupper(a$type[1])
  type.vec<-c("OTP","ITP","RTP","OSP","ISP")
  if(!(type%in%type.vec))
    stop("Illegal type code; valid types are:",paste(type.vec, collapse=","))
  dname<-"ddsp"
  
  if(!is.directed(nw)){  
    type <- "UTP"
    basenam<-"gwdsp"
    typecode<-0
  }else{
    typecode<-which(type==type.vec)
    basenam<-paste("gwdsp",type,sep=".")
  }

  linfo <- .sp.handle_layers(nw, a, type, FALSE, cache.sp)
  
  if(!fixed){ # This is a curved exponential family model
    if(!is.null(a$decay)) warning("In term 'dgwdsp': decay parameter 'decay' passed with 'fixed=FALSE'. 'decay' will be ignored. To specify an initial value for 'decay', use the 'init' control parameter.", call.=FALSE)

    #   d <- 1:(network.size(nw)-1)
    maxesp <- min(cutoff,network.size(nw)-2)
    d <- 1:maxesp
    ld<-length(d)
    if(ld==0){return(NULL)}
    
    params<-list(gwdsp=NULL,gwdsp.decay=decay)
    names(params)<-c(basenam,paste(basenam,"decay",sep="."))
    
    if(length(linfo)) c(list(name=paste0(dname,linfo$name_suffix),
                             coef.names=linfo$coef.namewrap(if(is.directed(nw)) paste("dsp.",type,"#",d,sep="") else paste("dsp#",d,sep="")),
                             iinputs=c(linfo$any_order,typecode,d), params=params,
                             auxiliaries = linfo$auxiliaries), GWDECAY)
    else no_layer_err("dgwdsp()")
  }else{
    if(is.null(a$decay)) stop("Term 'dgwdsp' with 'fixed=TRUE' requires a decay parameter 'decay'.", call.=FALSE)

    dname<-"dgwdsp"
    maxesp <- min(cutoff,network.size(nw)-2)
    if (is.directed(nw)) 
      coef.names <- paste("gwdsp",type,"fixed",decay,sep=".")
    else
      coef.names <- paste("gwdsp.fixed",decay,sep=".")
    
    if(length(linfo)) list(name=paste0(dname,linfo$name_suffix), coef.names=linfo$coef.namewrap(coef.names), inputs=decay, iinputs=c(linfo$any_order,typecode,maxesp), auxiliaries=linfo$auxiliaries)
    else no_layer_err("dgwdspL()")
  }
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
InitErgmTerm.dnspL<-function(nw, arglist, cache.sp=TRUE, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("d","type","L.base","Ls.path","L.in_order"),
                      vartypes = c("numeric","character","formula","formula,list","logical"),
                      defaultvalues = list(NULL,"OTP",NULL,NULL,FALSE),
                      required = c(TRUE, FALSE,FALSE,FALSE,FALSE))
  d<-a$d
  ld<-length(d)
  if(ld==0){return(NULL)}
  
  type<-toupper(a$type[1])
  type.vec<-c("OTP","ITP","RTP","OSP","ISP")
  if(!(type%in%type.vec))
    stop("Illegal type code for sp; valid types are:",paste(type.vec, collapse=","))
  dname<-"dnsp"
  if(is.directed(nw)){
    conam <- paste("nsp",type,sep=".")
    typecode<-which(type==type.vec)
  }else{
    conam<-"nsp"
    type<-"UTP"
    typecode<-0
  }

  linfo <- .sp.handle_layers(nw, a, type, TRUE, cache.sp)
  if(length(linfo)) nw <- linfo$nw1

  if (any(d==0)) {
    emptynwstats <- rep(0, length(d))
    if(is.bipartite(nw)){
      nb1 <- get.network.attribute(nw, "bipartite")
      nb2 <- network.size(nw) - nb1
      emptynwstats[d==0] <- nb1*(nb1-1)/2 + nb2*(nb2-1)/2
    }else{
      emptynwstats[d==0] <- network.dyadcount(nw,FALSE)
    }
  }else{
    emptynwstats <- NULL
  }
  if(length(linfo)) list(name=paste0(dname,linfo$name_suffix), coef.names=linfo$coef.namewrap(paste0(conam,d)), auxiliaries=linfo$auxiliaries, iinputs=c(linfo$any_order,typecode,d), minval=0, emptynwstats=emptynwstats)
  else no_layer_err("dnspL()")
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
InitErgmTerm.dgwnspL<-function(nw, arglist, cache.sp=TRUE, gw.cutoff=30, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("decay","fixed","cutoff","type", "alpha","L.base","Ls.path","L.in_order"),
                      vartypes = c("numeric","logical","numeric","character", "numeric","formula","formula,list","logical"),
                      defaultvalues = list(NULL, FALSE, gw.cutoff,"OTP", NULL,NULL,NULL,FALSE),
                      required = c(FALSE, FALSE, FALSE, FALSE, FALSE,FALSE,FALSE,FALSE))
  if(!is.null(a$alpha)){
    stop("For consistency with gw*degree terms, in all gw*sp and dgw*sp terms the argument ", sQuote("alpha"), " has been renamed to " ,sQuote("decay"), ".", call.=FALSE)
  }
  
  decay<-a$decay;fixed<-a$fixed
  cutoff<-a$cutoff
  decay=decay[1] # Not sure why anyone would enter a vector here, but...
  
  type<-toupper(a$type[1])
  type.vec<-c("OTP","ITP","RTP","OSP","ISP")
  if(!(type%in%type.vec))
    stop("Illegal type code; valid types are:",paste(type.vec, collapse=","))
  dname<-"dnsp"
  
  if(!is.directed(nw)){  
    type <- "UTP"
    basenam<-"gwdsp"
    typecode<-0
  }else{
    typecode<-which(type==type.vec)
    basenam<-paste("gwnsp",type,sep=".")
  }
  
  linfo <- .sp.handle_layers(nw, a, type, TRUE, cache.sp)

  if(!fixed){ # This is a curved exponential family model
    if(!is.null(a$decay)) warning("In term 'dgwnsp': decay parameter 'decay' passed with 'fixed=FALSE'. 'decay' will be ignored. To specify an initial value for 'decay', use the 'init' control parameter.", call.=FALSE)

    #   d <- 1:(network.size(nw)-1)
    maxesp <- min(cutoff,network.size(nw)-2)
    d <- 1:maxesp
    ld<-length(d)
    if(ld==0){return(NULL)}
    
    params<-list(gwnsp=NULL,gwnsp.decay=decay)
    names(params)<-c(basenam,paste(basenam,"decay",sep="."))
    
    if(length(linfo)) c(list(name=paste0(dname,linfo$name_suffix),
                             coef.names=linfo$coef.namewrap(if(is.directed(nw)) paste("nsp.",type,"#",d,sep="") else paste("nsp#",d,sep="")),
                             iinputs=c(linfo$any_order,typecode,d), params=params,
                             auxiliaries = linfo$auxiliaries), GWDECAY)
    else no_layer_err("dgwnsp()")
  }else{
    if(is.null(a$decay)) stop("Term 'dgwnsp' with 'fixed=TRUE' requires a decay parameter 'decay'.", call.=FALSE)

    dname<-"dgwnsp"
    maxesp <- min(cutoff,network.size(nw)-2)
    if (is.directed(nw)) 
      coef.names <- paste("gwnsp",type,"fixed",decay,sep=".")
    else
      coef.names <- paste("gwnsp.fixed",decay,sep=".")
    
    if(length(linfo)) list(name=paste0(dname,linfo$name_suffix), coef.names=linfo$coef.namewrap(coef.names), inputs=decay, iinputs=c(linfo$any_order,typecode,maxesp), auxiliaries=linfo$auxiliaries)
    else no_layer_err("dgwnspL()")
  }
}

#' @templateVar name dgwnspL
#' @template ergmTerm-rdname
#' @aliases gwnspL-ergmTerm
#' @description `gdwnspL` and `dgwnspL` are aliases for consistency with \pkg{ergm}.
#' @usage
#' # binary: gwnspL(decay, fixed=FALSE, cutoff=30, type="OTP", L.base=NULL,
#' #                Ls.path=NULL, L.in_order=FALSE)
InitErgmTerm.gwnspL <- InitErgmTerm.dgwnspL
