#  File R/InitErgmTerm.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################

#NOTE: a number of undocumented terms have been removed from this file
# the terms still exist on the experimental_terms svn branch


#===========================================================================
# This file contains the following 74 new, easier-to-write ergm-term
# initialization functions (each prepended with "InitErgmTerm"):
#   A:   <absdiff>          <absdiffcat>      <altkstar>
#        <asymmetric>       <adegcor>
#   B:   <b1concurrent>     <b1degree>        <b1factor>
#        <b1star>           <b1starmix>       <b1twostar>
#        <b2concurrent>     <b2degree>        <b2factor>         
#        <b2star>           <b2starmix>       <b2twostar>
#        <balance>
#   C:   <concurrent>       <cycle>           <ctriple>=<ctriad> 
#   D:   <degree>           <degreepopularity><density>         <dsp>
#        <dyadcov>          <degcrossprod>    <degcor>
#   E:   <edgecov>          <edges>           <esp>
#   G:   <gwb1degree>       <gwb2degree>      <gwdegree>
#        <gwdsp>            <gwesp>           <gwidegree>
#        <gwnsp>            <gwodegree>
#   H:   <hamming>          <hammingmix>
#   I:   <idegree>          <intransitive>    <idegreepopularity> 
#        <isolates>         <istar>
#   K:   <kstar>
#   L:   <localtriangle>
#   M:   <m2star>           <meandeg>         <mutual>
#   N:   <nearsimmelian>    <nodefactor>      <nodecov>=<nodemain> 
#        <nodeicov>         <nodeifactor>     <nodematch>=<match>
#        <nodemix>          <nodeocov>        <nodeofactor>       
#        <nsp>
#   O:   <odegree>          <opentriad>       <ostar>
#        <odegreepopularity>  
#   P:   <pdegcor>
#   R:   <receiver>         <rdegcor>
#   S:   <sender>           <simmelian>       <simmelianties>
#        <smalldiff>        <sociality>
#   T:   <threepath>        <transitive>      <triangles>=<triangle>
#        <triadcensus>      <tripercent>      <ttriple>=<ttriad>
#        <transitiveties>   <twopath
#==========================================================================

################################################################################
# The <InitErgmTerm.X> functions initialize each ergm term, X, by
#   1) checking the validity of X and its arguments via <check.ErgmTerm> and
#   2) setting appropiate values for each of the components in the returned list
# X is initialized for inclusion into a model that is specified by formula F and
# built via <ergm_model>
# 
# --PARAMETERS--
#   nw        : the network given in formula F
#   arglist   : the arguments given with term X in formula F
#   initialfit: whether the parameters for this term have been initially fit
#               (T or F); if FALSE, the ergm belongs to the curved exponential
#               family; if TRUE, the term X does not the ergm to be curved,
#               though other terms may; default=FALSE
#
# --IGNORED PARAMETERS--
#   ... : ignored, but necessary to accomodate other arguments
#         passed by <ergm_model>
#
# --RETURNED--
#   a list of term-specific elements required by the C changestats
#   functions and other R rountines; the first two components of this
#   list are required*, the remaining components are optional:
#     *name      : the name of term X; this is used to locate the C function
#                  calculating the change statistics for X, which will be
#                  'name' prepended with "d_"; for example if X=absdiff,
#                  'name'="absdiff", and the C function is "d_absdiff"
#     *coef.names: the vector of names for the coefficients (parameters)
#                  as they will be reported in the output
#     inputs     : the vector of (double-precision numeric) inputs that the 
#                  changestat function called d_<name> will require
#                  (see WHAT THE C CHANGESTAT FUNCTION RECEIVES below);
#                  this MUST be a vector!; thus, if the inputs are  matrices,
#                  they must be "flattened" to vectors; if they are categorical
#                  character-valued variables, they must be converted to numbers;
#                  optionally, 'inputs' may have an attribute named
#                  "ParamsBeforeCov",which is the number that used to be the
#                  old Element 1 and is needed for backwards compatability
#                  (see the old <InitErgm> for details); default=NULL
#     soname     : the name of the package containing the C function called
#                  d_'name'; default="ergm"
#     dependence : whether the addition of term X to the model makes the model
#                  into a dyadic dependence model (T or F); if all terms have
#                  'dependence' set FALSE, the model is assumed to be a
#                  dyadic independence model; default=TRUE
#    emptynwstats: the vector of values (if nonzero) for the statistics evaluated
#                  on the empty network; if all are zero for this term, this
#                  argument may be omitted.  Example:  If the degree0 term is
#                  among the statistics, this argument is unnecessary because
#                  degree0 = number of nodes for the empty network
#    minpar      : the vector of minimal valid values for each of the model's parameters
#    maxpar      : the vector of maximal valid values for each of the model's parameters
#    params      : a list of parameter values for curved exponential family model
#                  terms only; each item in the list should be named with the
#                  corresponding parameter name; those that coincide with the
#                  coef.names (used when initialfit=TRUE) will have their 'params'
#                  set by MPLE and their initial values in 'params' are ignored;
#                  otherwise, parameters should be given an initial value in this list
#    map         : a function taking two arguments, theta and length('params'), which
#                  gives the map from the canonical parameters, theta, to the curved
#                  parameters, eta; 'map' is only necessary for curved exponential
#                  family model terms
#   gradient     : a function taking two arguments, theta and length('params'), which
#                  gives the gradient of the eta map above as a p by q matrix, where
#                  p=length(theta), q=length(params); 'gradient' is only necessary
#                  for curved exponential family model terms
#
# WHAT THE C CHANGESTAT FUNCTION RECEIVES:
#                The changestat function, written in C and called d_'name',
#                will have access to 'inputs'; this array will be called INPUT_PARAMS
#                in the C code and its entries may accessed as INPUT_PARAMS[0],
#                INPUT_PARAMS[1], and so on; the size of INPUT_PARAMS=N_INPUT_PARAMS,
#                a value which is automatically set for you and which is available
#                inside the C function; thus INPUT_PARAMS[N_INPUT_PARAMS-1] is the last
#                element in the vector; note in particular that it is NOT necessary 
#                to add the number of inputs to 'inputs' since this is done automatically
#
################################################################################

#=======================InitErgmTerm utility functions============================#

GWDECAY <- list(
  map = function(x,n,...) {
    i <- 1:n
    x[1]*(exp(x[2])*(1-(1-exp(-x[2]))^i))
  },
  gradient = function(x,n,...) {
    i <- 1:n
    e2 <- exp(x[2])
    a <- 1-exp(-x[2])
    rbind((1-a^i)*e2, x[1] * ( (1-a^i)*e2 - i*a^(i-1) ) )
  },
  minpar = c(-Inf, 0)
)

.process_layers_degree <- function(nw, a, name, coef.names, inputs, emptynwstats=NULL){
  if(!is.null(a$Ls)){
    Ls <- .set_layer_namemap(a$Ls)
    if(is(Ls,"formula")) Ls <- list(Ls)
    auxiliaries <- .mk_.layer.net_auxform(Ls)
    if(!is.null(a$dir)){
      dir <- rep(a$dir, length.out=length(Ls))
      dir <- pmatch(a$dir, c("in","all","out"))-2
      if(any(is.na(dir) | dir==0)) stop("Invalid direction specification.")
    } else dir <- NULL
    inputs <- c(length(Ls), dir, inputs)
    if(!is.null(emptynwstats)) emptynwstats <- emptynwstats / length(unique(.peek_vattrv(nw, ".LayerID")))
    name <- paste0(name,"_ML_sum")
    coef.names <- .lspec_coef.namewrap(Ls)(coef.names)
  }else ergm_Init_abort("Use the non-layer version.")

  list(name = name, coef.names = coef.names, inputs = inputs, emptynwstats = emptynwstats, auxiliaries=auxiliaries, minval=0, maxval=network.size(nw), dependence=TRUE)
}

################################################################################

#' @templateVar name b1degreeL
#' @title Degree for the first mode in a bipartite (aka two-mode) network
#' @description Degree for the first mode in a bipartite (aka two-mode) network
#' @details This term adds one network statistic to the model for
#'   each element in `d` ; the \eqn{i} th such statistic equals the number of
#'   nodes of degree `d[i]` in the first mode of a bipartite network, i.e.
#'   with exactly `d[i]` edges. The first mode of a bipartite network object
#'   is sometimes known as the "actor" mode.
#'
#'   This term can only be used with undirected bipartite networks.
#'
#' @usage
#' # binary: b1degreeL(d, by=NULL, levels=NULL)
#' @param d a vector of distinct integers.
#' @param by a character string giving the name of an attribute in the network's vertex
#'   attribute list. If this is specified
#'   then each node's degree is tabulated only with other nodes having the same
#'   value of the `by` attribute.
#' @param levels a list of layer specifications. If given, degree of a node
#'   `i` is considered to be the number of edges in all layers,
#'   combined.
#'
#' @template ergmTerm-general
#'
#' @concept bipartite
#' @concept undirected
#' @concept categorical nodal attribute
#' @concept frequently-used
InitErgmTerm.b1degreeL <- function(nw, arglist, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm (nw, arglist, directed=FALSE, bipartite=TRUE,
                       varnames = c("d", "by", "levels", "Ls"),
                       vartypes = c("numeric", "character", "character,numeric,logical", "formula,list"),
                       defaultvalues = list(NULL, NULL, NULL, NULL),
                       required = c(TRUE, FALSE, FALSE, FALSE))
  ### Process the arguments
  nb1 <- get.network.attribute(nw, "bipartite")
  if (!is.null(a$by)) {  # CASE 1:  a$by GIVEN
    nodecov <- get.node.attr(nw, a$by)[seq_len(nb1)]
    u <- NVL(a$levels, sort(unique(nodecov)))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u,nomatch=length(u)+1) # Recode to numeric
    # Combine degree and u into 2xk matrix, where k=length(a$d)*length(u)
    lu <- length(u)
    du <- rbind(rep(a$d,lu), rep(1:lu, rep(length(a$d), lu)))
    emptynwstats <- rep(0, ncol(du))
    if (any(du[1,]==0)) { # Alter emptynwstats
      tmp <- du[2,du[1,]==0]
      for(i in 1:length(tmp)) 
        tmp[i] <- sum(nodecov==tmp[i])
      emptynwstats[du[1,]==0] <- tmp
    }
    name <- "b1degree_by_attr"
    coef.names <- paste("b1deg", du[1,], ".", a$by, u[du[2,]], sep="")
    inputs <- c(as.vector(du), nodecov)
  } else { # CASE 2:  a$by NOT GIVEN
    name <- "b1degree"
    coef.names <- paste("b1deg", a$d, sep="")
    inputs <- a$d
    emptynwstats <- rep(0, length(a$d))
    if (any(a$d==0)) { # alter emptynwstats
      emptynwstats[a$d==0] <- nb1
    }
  }

  c(.process_layers_degree(nw, a, name, coef.names, inputs, emptynwstats),
    minval = 0, maxval=nb1, conflicts.constraints="odegreedist")
}

################################################################################

#' @templateVar name b2degree
#' @title Degree for the second mode in a bipartite (aka two-mode) network
#' @description Degree for the second mode in a bipartite (aka two-mode) network
#' @details This term adds one network statistic to the model for
#'   each element in `d` ; the \eqn{i} th such statistic equals the number of
#'   nodes of degree `d[i]` in the second mode of a bipartite network, i.e.
#'   with exactly `d[i]` edges. The second mode of a bipartite network
#'   object is sometimes known as the "event" mode. The optional term
#'   `by` is a character string giving the name of an attribute in the
#'   network's vertex attribute list.
#'   This term can only be used with undirected bipartite networks.
#'
#' @usage
#' # binary: b2degree(d, by=NULL)
#' @param d a vector of distinct integers
#' @param by If this is specified
#'   then each node's degree is tabulated only with other nodes having the same
#'   value of the given attribute.
#' @param levels a list of layer specifications. If given, degree of a node
#'   `i` is considered to be the number of edges in all layers,
#'   combined.
#'
#' @template ergmTerm-general
#'
#' @concept bipartite
#' @concept undirected
#' @concept categorical nodal attribute
#' @concept frequently-used
InitErgmTerm.b2degreeL <- function(nw, arglist, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm (nw, arglist, directed=FALSE, bipartite=TRUE,
                       varnames = c("d", "by", "levels", "Ls"),
                       vartypes = c("numeric", "character", "character,numeric,logical", "formula,list"),
                       defaultvalues = list(NULL, NULL, NULL, NULL),
                       required = c(TRUE, FALSE, FALSE, FALSE))
  ### Process the arguments
  nb1 <- get.network.attribute(nw, "bipartite")
  n <- network.size(nw)
  if (!is.null(a$by)) {  # CASE 1:  a$by GIVEN
    nodecov <- get.node.attr(nw, a$by)[-seq_len(nb1)]
    u <- NVL(a$levels, sort(unique(nodecov)))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u,nomatch=length(u)+1) # Recode to numeric
    # Combine degree and u into 2xk matrix, where k=length(a$d)*length(u)
    lu <- length(u)
    du <- rbind(rep(a$d,lu), rep(1:lu, rep(length(a$d), lu)))
    emptynwstats <- rep(0, ncol(du))
    if (any(du[1,]==0)) { # Alter emptynwstats
      tmp <- du[2,du[1,]==0]
      for(i in 1:length(tmp)) 
        tmp[i] <- sum(nodecov==tmp[i])
      emptynwstats[du[1,]==0] <- tmp
    }
    name <- "b2degree_by_attr"
    coef.names <- paste("b2deg", du[1,], ".", a$by, u[du[2,]], sep="")
    inputs <- c(as.vector(du), nodecov)
  } else { # CASE 2:  a$by NOT GIVEN
    name <- "b2degree"
    coef.names <- paste("b2deg", a$d, sep="")
    inputs <- a$d
    emptynwstats <- rep(0, length(a$d))
    if (any(a$d==0)) { # alter emptynwstats
      emptynwstats[a$d==0] <- n-nb1
    }
  }

  c(.process_layers_degree(nw, a, name, coef.names, inputs, emptynwstats),
    minval = 0, maxval=network.size(nw)-nb1, conflicts.constraints="b2degreedist")
}

################################################################################

#' @templateVar name degreeL
#' @title Degree
#' @description Degree
#' @details This term adds one
#'   network statistic to the model for each element in `d` ; the \eqn{i} th
#'   such statistic equals the number of nodes in the network of degree
#'   `d[i]` , i.e. with exactly `d[i]` edges.
#'
#'   This term can only be used with undirected networks; for directed networks
#'   see `idegree` and `odegree` .
#'
#' @usage
#' # binary: degreeL(d, by=NULL, homophily=FALSE, levels=NULL, Ls=NULL)
#' @param d a vector of distinct integers
#' @param by a character string giving the name of an attribute in the
#'   network's vertex attribute list.
#' @param homophily If this is specified and `homophily` is `TRUE` ,
#'   then degrees are calculated using the subnetwork consisting of only
#'   edges whose endpoints have the same value of the `by` attribute.
#'   If `by` is specified and
#'   `homophily` is `FALSE` (the default), then separate degree
#'   statistics are calculated for nodes having each separate
#'   value of the attribute.
#' @param Ls if specified, degree of a node
#'   `i` is considered to be the number of edges in all layers,
#'   combined.
#'
#' @template ergmTerm-general
#'
#' @concept undirected
#' @concept categorical nodal attribute
#' @concept frequently-used
#' @concept directed
InitErgmTerm.degreeL<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=FALSE,
                      varnames = c("d", "by", "homophily", "levels", "Ls"),
                      vartypes = c("numeric", "character", "logical", "character,numeric,logical", "formula,list"),
                      defaultvalues = list(NULL, NULL, FALSE, NULL, NULL),
                      required = c(TRUE, FALSE, FALSE, FALSE, FALSE))
  d<-a$d; byarg <- a$by; homophily <- a$homophily
  emptynwstats<-NULL
  if(!is.null(byarg)) {
    nodecov <- get.node.attr(nw, byarg, "degree")
    u <- NVL(a$levels, sort(unique(nodecov)))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u,nomatch=length(u)+1) # Recode to numeric
    if (length(u)==1)
         stop ("Attribute given to degree() has only one value", call.=FALSE)
  }
  if(!is.null(byarg) && !homophily) {
    # Combine degree and u into 2xk matrix, where k=length(d)*length(u)
    lu <- length(u)
    du <- rbind(rep(d,lu), rep(1:lu, rep(length(d), lu)))
    if (any(du[1,]==0)) {
      emptynwstats <- rep(0, ncol(du))
      tmp <- du[2,du[1,]==0]
      for(i in 1:length(tmp)) tmp[i] <- sum(nodecov==tmp[i])
        emptynwstats[du[1,]==0] <- tmp
    }
  } else {
    if (any(d==0)) {
      emptynwstats <- rep(0, length(d))
      emptynwstats[d==0] <- network.size(nw)
    }
  }
  if(is.null(byarg)) {
    if(length(d)==0){return(NULL)}
    coef.names <- paste("degree",d,sep="")
    name <- "degree"
    inputs <- c(d)
  } else if (homophily) {
    if(length(d)==0){return(NULL)}
    # See comment in d_degree_w_homophily function
    coef.names <- paste("deg", d, ".homophily.",byarg, sep="")
    name <- "degree_w_homophily"
    inputs <- c(d, nodecov)
  } else {
    if(ncol(du)==0) {return(NULL)}
    #  No covariates here, so "ParamsBeforeCov" unnecessary
    # See comment in d_degree_by_attr function
    coef.names <- paste("deg", du[1,], ".", byarg,u[du[2,]], sep="")
    name <- "degree_by_attr"
    inputs <- c(as.vector(du), nodecov)
  }

  c(.process_layers_degree(nw, a, name, coef.names, inputs, emptynwstats),
    minval = 0, maxval=network.size(nw), conflicts.constraints="degreedist")
}

#=======================InitErgmTerm functions:  G============================#

################################################################################

#' @templateVar name gwb1degreeL
#' @title Geometrically weighted degree distribution for the first mode in a bipartite (aka two-mode) network
#' @description Geometrically weighted degree distribution for the first mode in a bipartite (aka two-mode) network
#' @details This term adds one network statistic to the model equal to the weighted
#'   degree distribution with decay controlled by the `decay` parameter, which should be non-negative,
#'   for nodes in the
#'   first mode of a bipartite network. The first mode of a bipartite network
#'   object is sometimes known as the "actor" mode.
#'
#'   This term can only be used with undirected bipartite
#'   networks.
#'
#' @usage
#' # binary: gwb1degreeL(decay, fixed=FALSE, cutoff=30, levels=NULL, Ls=NULL)
#' @param decay non-negative model parameter that is the same as theta_s in
#'   equation (14) in Hunter (2007).
#' @param fixed specify if the value supplied for `decay` may be fixed (if `fixed=TRUE` ),
#'   or it may be used as merely the starting value for the estimation
#'   in a curved exponential family model (the default).
#' @param attrname if specified, then separate degree
#'   statistics are calculated for nodes having each separate
#'   value of the attribute.
#' @param cutoff only relevant if `fixed=FALSE` . In that case it only uses this
#'   number of terms in computing the statistics to reduce the computational
#'   burden. Its default value can also be controlled by the `gw.cutoff` term option control parameter. (See [`control.ergm`] .)
#' @param levels a list of layer specifications. If given, degree of a node
#'   `i` is considered to be the number of edges in all layers,
#'   combined.
#' @param Ls TODO
#'
#' @template ergmTerm-general
#'
#' @concept bipartite
#' @concept undirected
#' @concept curved
InitErgmTerm.gwb1degreeL<-function(nw, arglist, gw.cutoff=30, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite=TRUE,                     
  # default for 'fixed' should be made 'FALSE' when the function can handle it!                    
                      varnames = c("decay", "fixed", "attrname","cutoff", "levels", "Ls"),
                      vartypes = c("numeric", "logical", "character","numeric", "character,numeric,logical", "formula,list"),
                      defaultvalues = list(NULL, FALSE, NULL, gw.cutoff, NULL, NULL),
                      required = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE))
  decay<-a$decay; fixed<-a$fixed; attrname<-a$attrname
  cutoff<-a$cutoff
  nb1 <- get.network.attribute(nw,"bipartite")
# d <- 1:(network.size(nw) - nb1)
  maxesp <- min(cutoff, network.size(nw)-nb1)

  d <- 1:maxesp
  if (!is.null(attrname) && !fixed) {
      stop("The gwb1degree term cannot yet handle a nonfixed decay ",
           "term with an attribute. Use fixed=TRUE.", call.=FALSE)
  }
  if(!fixed){# This is a curved exponential family model
    if(!is.null(a$decay)) warning("In term 'gwb1degree': decay parameter 'decay' passed with 'fixed=FALSE'. 'decay' will be ignored. To specify an initial value for 'decay', use the 'init' control parameter.", call.=FALSE)
    ld<-length(d)
    if(ld==0){return(NULL)}
    c(.process_layers_degree(nw, a, name="b1degree", coef.names=paste("gwb1degree#",d,sep=""), inputs=c(d)),
      params=list(list(gwb1degree=NULL,gwb1degree.decay=decay)),
      GWDECAY, conflicts.constraints="b1degreedist")
  } else {
    if(is.null(a$decay)) stop("Term 'gwb1degree' with 'fixed=TRUE' requires a decay parameter 'decay'.", call.=FALSE)

    if(!is.null(attrname)) {
      nodecov <- get.node.attr(nw, attrname, "gwb1degree")
      u <- NVL(a$levels, sort(unique(nodecov)))
      if(any(is.na(nodecov))){u<-c(u,NA)}
      nodecov <- match(nodecov,u,nomatch=length(u)+1) # Recode to numeric
      if (length(u)==1)
        stop ("Attribute given to gwb1degree() has only one value", call.=FALSE)
      # Combine degree and u into 2xk matrix, where k=length(d)*length(u)
      lu <- length(u)
      du <- rbind(rep(d,lu), rep(1:lu, rep(length(d), lu)))
      if(nrow(du)==0) {return(NULL)}
      #  No covariates here, so "ParamsBeforeCov" unnecessary
      # See comment in d_gwb1degree_by_attr function
      name <- "gwb1degree_by_attr"
      coef.names <- paste("gwb1deg", decay, ".",attrname, u, sep="")
      inputs <- c(decay, nodecov)
    }else{
      name <- "gwb1degree"
      coef.names <- paste("gwb1deg.fixed.",decay,sep="")
      inputs <- c(decay)
    }
    c(.process_layers_degree(nw, a, name=name, coef.names=coef.names, inputs=inputs), conflicts.constraints="b1degreedist")
  }
}



################################################################################

#' @templateVar name gwb2degreeL
#' @title Geometrically weighted degree distribution for the second mode in a bipartite (aka two-mode) network
#' @description Geometrically weighted degree distribution for the second mode in a bipartite (aka two-mode) network
#' @details This term adds one network statistic to the model equal to the weighted
#'   degree distribution with decay controlled by the which should be non-negative,
#'   for nodes in the
#'   second mode of a bipartite network. The second mode of a bipartite network
#'   object is sometimes known as the "event" mode.
#'
#'   This term can only be used with undirected bipartite
#'   networks.
#'
#' @usage
#' # binary: gwb2degreeL(decay, fixed=FALSE, attrname=NULL, cutoff=30, levels=NULL, Ls=NULL)
#' @param decay non-negative model parameter that is the same as theta_s in
#'   equation (14) in Hunter (2007).
#' @param fixed specify if the value supplied for `decay` may be fixed (if `fixed=TRUE` ),
#'   or it may be used as merely the starting value for the estimation
#'   in a curved exponential family model (the default).
#' @param attrname if specified, then separate degree
#'   statistics are calculated for nodes having each separate
#'   value of the attribute.
#' @param cutoff only relevant if `fixed=FALSE` . In that case it only uses this
#'   number of terms in computing the statistics to reduce the computational
#'   burden. Its default value can also be controlled by the `gw.cutoff` term option control parameter. (See [`control.ergm`] .)
#' @param levels a list of layer specifications. If given, degree of a node
#'   `i` is considered to be the number of edges in all layers,
#'   combined.
#' @param Ls TODO
#'
#' @template ergmTerm-general
#'
#' @concept bipartite
#' @concept undirected
#' @concept curved
InitErgmTerm.gwb2degreeL<-function(nw, arglist, gw.cutoff=30,  ...) {
  a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite=TRUE,
  # default for 'fixed' should be made 'FALSE' when the function can handle it!                    
                      varnames = c("decay", "fixed", "attrname","cutoff", "levels", "Ls"),
                      vartypes = c("numeric", "logical", "character", "numeric", "character,numeric,logical", "formula,list"),
                      defaultvalues = list(NULL, FALSE, NULL, gw.cutoff, NULL, NULL),
                      required = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE))
  decay<-a$decay; fixed<-a$fixed; attrname<-a$attrname
  cutoff<-a$cutoff
  nb1 <- get.network.attribute(nw,"bipartite")
# d <- 1:nb1
  maxesp <- min(cutoff,nb1)
  d <- 1:maxesp
  if (!is.null(attrname) && !fixed) {
      stop("The gwb2degree term cannot yet handle a nonfixed decay ",
           "term with an attribute. Use fixed=TRUE.", call.=FALSE) }
  if(!fixed){# This is a curved exponential family model
    if(!is.null(a$decay)) warning("In term 'gwb2degree': decay parameter 'decay' passed with 'fixed=FALSE'. 'decay' will be ignored. To specify an initial value for 'decay', use the 'init' control parameter.", call.=FALSE)
    ld<-length(d)
    if(ld==0){return(NULL)}
    c(.process_layers_degree(nw, a, name="b2degree", coef.names=paste("gwb2degree#",d,sep=""), inputs=c(d)),
      params=list(list(gwb2degree=NULL,gwb2degree.decay=decay)),
      GWDECAY, conflicts.constraints="b2degreedist")
  } else { 
    if(is.null(a$decay)) stop("Term 'gwb2degree' with 'fixed=TRUE' requires a decay parameter 'decay'.", call.=FALSE)

    if(!is.null(attrname)) {
      nodecov <- get.node.attr(nw, attrname, "gwb2degree")
      u <- NVL(a$levels, sort(unique(nodecov)))
      if(any(is.na(nodecov))){u<-c(u,NA)}
      nodecov <- match(nodecov,u,nomatch=length(u)+1) # Recode to numeric
      if (length(u)==1)
        stop ("Attribute given to gwb2degree() has only one value", call.=FALSE)
      # Combine degree and u into 2xk matrix, where k=length(d)*length(u)
      lu <- length(u)
      du <- rbind(rep(d,lu), rep(1:lu, rep(length(d), lu)))
      if(nrow(du)==0) {return(NULL)}
     #  No covariates here, so "ParamsBeforeCov" unnecessary
     # See comment in d_gwb2degree_by_attr function
      name <- "gwb2degree_by_attr"
      coef.names <- paste("gwb2deg", decay, ".", attrname, u, sep="")
      inputs <- c(decay, nodecov)
    }else{
      name <- "gwb2degree"
      coef.names <- paste("gwb2deg.fixed.",decay,sep="")
      inputs <- c(decay)
    }
    c(.process_layers_degree(nw, a, name=name, coef.names=coef.names, inputs=inputs), conflicts.constraints="b2degreedist")
  }
}

################################################################################

#' @templateVar name gwdegreeL
#' @title Geometrically weighted degree distribution
#' @description Geometrically weighted degree distribution
#' @details This term adds one network statistic to the model equal to the weighted
#'   degree distribution with decay controlled by the `decay` parameter.
#'
#'   This term can only be used with undirected networks.
#'
#' @usage
#' # binary: gwdegreeL(decay, fixed=FALSE, attrname=NULL, cutoff=30, levels=NULL)
#' @param decay non-negative model parameter that is the same as theta_s in
#'   equation (14) in Hunter (2007).
#' @param fixed specify if the value supplied for `decay` may be fixed (if `fixed=TRUE` ),
#'   or it may be used as merely the starting value for the estimation
#'   in a curved exponential family model (the default).
#' @param attrname if specified, then separate degree
#'   statistics are calculated for nodes having each separate
#'   value of the attribute.
#' @param cutoff only relevant if `fixed=FALSE` . In that case it only uses this
#'   number of terms in computing the statistics to reduce the computational
#'   burden. Its default value can also be controlled by the `gw.cutoff` term option control parameter. (See [`control.ergm`] .)
#' @param levels a list of layer specifications. If given, degree of a node
#'   `i` is considered to be the number of edges in all layers,
#'   combined.
#'
#' @template ergmTerm-general
#'
#' @concept undirected
#' @concept curved
#' @concept frequently-used
InitErgmTerm.gwdegreeL<-function(nw, arglist, gw.cutoff=30, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=FALSE,
                      varnames = c("decay", "fixed", "attrname","cutoff", "levels", "Ls"),
                      vartypes = c("numeric", "logical", "character", "numeric", "character,numeric,logical", "formula,list"),
                      defaultvalues = list(NULL, FALSE, NULL, gw.cutoff, NULL, NULL),
                      required = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE))
  decay<-a$decay; attrname<-a$attrname; fixed<-a$fixed  
  cutoff<-a$cutoff
# d <- 1:(network.size(nw)-1)
   maxesp <- min(cutoff,network.size(nw)-1)
  d <- 1:maxesp
  if (!is.null(attrname) && !fixed) {
    stop("The gwdegree term cannot yet handle a nonfixed decay ",
            "term with an attribute. Use fixed=TRUE.", call.=FALSE)
  }
  if(!fixed){ # This is a curved exponential family model
    if(!is.null(a$decay)) warning("In term 'gwdegree': decay parameter 'decay' passed with 'fixed=FALSE'. 'decay' will be ignored. To specify an initial value for 'decay', use the 'init' control parameter.", call.=FALSE)
    ld<-length(d)
    if(ld==0){return(NULL)}
    c(.process_layers_degree(nw, a, name="degree", coef.names=paste("gwdegree#",d,sep=""), inputs=c(d)),
      params=list(list(gwdegree=NULL,gwdegree.decay=decay)),
      GWDECAY, conflicts.constraints="degreedist")
  } else {
    if(is.null(a$decay)) stop("Term 'gwdegree' with 'fixed=TRUE' requires a decay parameter 'decay'.", call.=FALSE)

    if(!is.null(attrname)) {
      nodecov <- get.node.attr(nw, attrname, "gwdegree")
      u <- NVL(a$levels, sort(unique(nodecov)))
      if(any(is.na(nodecov))){u<-c(u,NA)}
      nodecov <- match(nodecov,u,nomatch=length(u)+1) # Recode to numeric
      if (length(u)==1)
        stop ("Attribute given to gwdegree() has only one value", call.=FALSE)
      # Combine degree and u into 2xk matrix, where k=length(d)*length(u)
      lu <- length(u)
      du <- rbind(rep(d,lu), rep(1:lu, rep(length(d), lu)))
      if(nrow(du)==0) {return(NULL)}
      #  No covariates here, so "ParamsBeforeCov" unnecessary
      name <- "gwdegree_by_attr"
      coef.names <- paste("gwdeg", decay, ".", attrname, u, sep="")
      inputs <- c(decay, nodecov)
    }else{
      name <- "gwdegree"
      coef.names <- paste("gwdeg.fixed.",decay,sep="")
      inputs <- c(decay)
    }
    c(.process_layers_degree(nw, a, name=name, coef.names=coef.names, inputs=inputs), conflicts.constraints="degreedist")
  }
}

################################################################################

#' @templateVar name gwidegreeL
#' @title Geometrically weighted in-degree distribution
#' @description Geometrically weighted in-degree distribution
#' @details This term adds one network statistic to the model
#'   equal to the weighted in-degree distribution with decay parameter. This
#'   term can only be used with directed networks.
#'
#' @usage
#' # binary: gwidegreeL(decay, fixed=FALSE, attrname=NULL, cutoff=30, levels=NULL, Ls=NULL)
#' @param decay non-negative model parameter that is the same as theta_s in
#'   equation (14) in Hunter (2007).
#' @param fixed specify if the value supplied for `decay` may be fixed (if `fixed=TRUE` ),
#'   or it may be used as merely the starting value for the estimation
#'   in a curved exponential family model (the default).
#' @param attrname if specified, then separate degree
#'   statistics are calculated for nodes having each separate
#'   value of the attribute.
#' @param cutoff only relevant if `fixed=FALSE` . In that case it only uses this
#'   number of terms in computing the statistics to reduce the computational
#'   burden. Its default value can also be controlled by the `gw.cutoff` term option control parameter. (See [`control.ergm`] .)
#' @param levels a list of layer specifications. If given, degree of a node
#'   `i` is considered to be the number of edges in all layers,
#'   combined.
#'
#' @template ergmTerm-general
#'
#' @concept directed
#' @concept curved
InitErgmTerm.gwidegreeL<-function(nw, arglist, gw.cutoff=30,  ...) {
  a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                      varnames = c("decay", "fixed", "attrname","cutoff", "levels", "Ls"),
                      vartypes = c("numeric", "logical", "character","numeric", "character,numeric,logical", "formula,list"),
                      defaultvalues = list(NULL, FALSE, NULL, gw.cutoff, NULL, NULL),
                      required = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE))
  decay<-a$decay; attrname<-a$attrname; fixed<-a$fixed  
  cutoff<-a$cutoff
# d <- 1:(network.size(nw)-1)
  maxesp <- min(cutoff,network.size(nw)-1)
  d <- 1:maxesp
  if (!is.null(attrname) && !fixed ) {
    stop("The gwidegree term cannot yet handle a nonfixed decay ",
            "term with an attribute. Use fixed=TRUE.", call.=FALSE)
    
  }
  if(!fixed){ # This is a curved exponential family model
    if(!is.null(a$decay)) warning("In term 'gwidegree': decay parameter 'decay' passed with 'fixed=FALSE'. 'decay' will be ignored. To specify an initial value for 'decay', use the 'init' control parameter.", call.=FALSE)
    ld<-length(d)
    if(ld==0){return(NULL)}
    c(.process_layers_degree(nw, a, name="idegree", coef.names=paste("gwidegree#",d,sep=""), inputs=c(d)),
      params=list(list(gwidegree=NULL,gwidegree.decay=decay)),
      GWDECAY, conflicts.constraints="idegreedist")
  } else { 
    if(is.null(a$decay)) stop("Term 'gwidegree' with 'fixed=TRUE' requires a decay parameter 'decay'.", call.=FALSE)

    if(!is.null(attrname)) {
      nodecov <- get.node.attr(nw, attrname, "gwidegree")
      u <- NVL(a$levels, sort(unique(nodecov)))
      if(any(is.na(nodecov))){u<-c(u,NA)}
      nodecov <- match(nodecov,u,nomatch=length(u)+1) # Recode to numeric
      if (length(u)==1)
        stop ("Attribute given to gwidegree() has only one value", call.=FALSE)
      # Combine degree and u into 2xk matrix, where k=length(d)*length(u)
      lu <- length(u)
      du <- rbind(rep(d,lu), rep(1:lu, rep(length(d), lu)))
      if(nrow(du)==0) {return(NULL)}
      #  No covariates here, so "ParamsBeforeCov" unnecessary
      name <- "gwidegree_by_attr"
      coef.names <- paste("gwideg", decay, ".", attrname, u, sep="")
      inputs <- c(decay, nodecov)
    }else{
      name <- "gwidegree"
      coef.names <- paste("gwideg.fixed.",decay,sep="")
      inputs <- c(decay)
    }
    c(.process_layers_degree(nw, a, name=name, coef.names=coef.names, inputs=inputs), conflicts.constraints="idegreedist")
  }
}

################################################################################

#' @templateVar name gwodegreeL
#' @title Geometrically weighted out-degree distribution
#' @description Geometrically weighted out-degree distribution
#' @details This term adds one network statistic to the model
#'   equal to the weighted out-degree distribution with decay parameter . This
#'   term can only be used with directed networks.
#'
#' @usage
#' # binary: gwodegreeL(decay, fixed=FALSE, attrname=NULL, cutoff=30, levels=NULL, Ls=NULL)
#' @param decay non-negative model parameter that is the same as theta_s in
#'   equation (14) in Hunter (2007).
#' @param fixed specify if the value supplied for `decay` may be fixed (if `fixed=TRUE` ),
#'   or it may be used as merely the starting value for the estimation
#'   in a curved exponential family model (the default).
#' @param attrname if specified, then separate degree
#'   statistics are calculated for nodes having each separate
#'   value of the attribute.
#' @param cutoff only relevant if `fixed=FALSE` . In that case it only uses this
#'   number of terms in computing the statistics to reduce the computational
#'   burden. Its default value can also be controlled by the `gw.cutoff` term option control parameter. (See [`control.ergm`] .)
#' @param levels a list of layer specifications. If given, degree of a node
#'   `i` is considered to be the number of edges in all layers,
#'   combined.
#' @param Ls TODO
#'
#' @template ergmTerm-general
#'
#' @concept directed
#' @concept curved
InitErgmTerm.gwodegreeL<-function(nw, arglist, gw.cutoff=30, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                      varnames = c("decay", "fixed", "attrname","cutoff", "levels", "Ls"),
                      vartypes = c("numeric", "logical", "character","numeric", "character,numeric,logical", "formula,list"),
                      defaultvalues = list(NULL, FALSE, NULL, gw.cutoff, NULL, NULL),
                      required = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE))
  decay<-a$decay; attrname<-a$attrname; fixed<-a$fixed  
  cutoff<-a$cutoff
# d <- 1:(network.size(nw)-1)
   maxesp <- min(cutoff,network.size(nw)-1)
  d <- 1:maxesp
  if (!is.null(attrname) && !fixed ) {
    stop("The gwodegree term cannot yet handle a nonfixed decay ",
            "term with an attribute. Use fixed=TRUE.", call.=FALSE)
    
  }
  if(!fixed){ # This is a curved exponential family model
    if(!is.null(a$decay)) warning("In term 'gwodegree': decay parameter 'decay' passed with 'fixed=FALSE'. 'decay' will be ignored. To specify an initial value for 'decay', use the 'init' control parameter.", call.=FALSE)
    ld<-length(d)
    if(ld==0){return(NULL)}
    c(.process_layers_degree(nw, a, name="odegree", coef.names=paste("gwodegree#",d,sep=""), inputs=c(d)),
      params=list(list(gwodegree=NULL,gwodegree.decay=decay)),
      GWDECAY, conflicts.constraints="odegreedist")
  } else {
    if(is.null(a$decay)) stop("Term 'gwodegree' with 'fixed=TRUE' requires a decay parameter 'decay'.", call.=FALSE)

    if(!is.null(attrname)) {
      nodecov <- get.node.attr(nw, attrname, "gwodegree")
      u <- NVL(a$levels, sort(unique(nodecov)))
      if(any(is.na(nodecov))){u<-c(u,NA)}
      nodecov <- match(nodecov,u,nomatch=length(u)+1) # Recode to numeric
      if (length(u)==1)
        stop ("Attribute given to gwodegree() has only one value", call.=FALSE)
      # Combine degree and u into 2xk matrix, where k=length(d)*length(u)
      lu <- length(u)
      du <- rbind(rep(d,lu), rep(1:lu, rep(length(d), lu)))
      if(nrow(du)==0) {return(NULL)}
      #  No covariates here, so "ParamsBeforeCov" unnecessary
      name <- "gwodegree_by_attr"
      coef.names <- paste("gwodeg", decay, ".", attrname, u, sep="")
      inputs <- c(decay, nodecov)
    }else{
      name <- "gwodegree"
      coef.names <- paste("gwodeg.fixed.",decay,sep="")
      inputs <- c(decay)
    }
    c(.process_layers_degree(nw, a, name=name, coef.names=coef.names, inputs=inputs), conflicts.constraints="odegreedist")
  }
}

################################################################################

#' @templateVar name idegreeL
#' @title In-degree
#' @description In-degree
#' @details This term adds one network statistic to
#'   the model for each element in `d` ; the \eqn{i} th such statistic equals
#'   the number of nodes in the network of in-degree `d[i]` , i.e. the number
#'   of nodes with exactly `d[i]` in-edges.
#'   This term can only be used with directed networks; for undirected networks
#'   see `degree` .
#'
#' @usage
#' # binary: idegreeL(d, by=NULL, homophily=FALSE, levels=NULL)
#' @param d a vector of distinct integers.
#' @param by an optional character string giving the name of an attribute in the
#'   network's vertex attribute list.
#' @param homophily only applied if by is specified. If set (`homophile == TRUE`),
#'   then degrees are calculated using the subnetwork consisting of only
#'   edges whose endpoints have the same value of the `by` attribute.
#'   Otherwise (the default), then separate degree
#'   statistics are calculated for nodes having each separate
#'   value of the attribute.
#' @param levels a list of layer specifications. If given, degree of a node
#'   `i` is considered to be the number of edges in all layers, combined.
#'
#' @template ergmTerm-general
#'
#' @concept directed
#' @concept categorical nodal attribute
#' @concept frequently-used
InitErgmTerm.idegreeL<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                      varnames = c("d", "by", "homophily", "levels", "Ls"),
                      vartypes = c("numeric", "character", "logical", "character,numeric,logical", "formula,list"),
                      defaultvalues = list(NULL, NULL, FALSE, NULL, NULL),
                      required = c(TRUE, FALSE, FALSE, FALSE, FALSE))
  d<-a$d; byarg <- a$by; homophily <- a$homophily
  emptynwstats<-NULL
  if(!is.null(byarg)) {
    nodecov <- get.node.attr(nw, byarg, "idegree")
    u <- NVL(a$levels, sort(unique(nodecov)))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u,nomatch=length(u)+1) # Recode to numeric
    if (length(u)==1)
         stop ("Attribute given to idegree() has only one value", call.=FALSE)
  }
  if(!is.null(byarg) && !homophily) {
    # Combine degree and u into 2xk matrix, where k=length(d)*length(u)
    lu <- length(u)
    du <- rbind(rep(d,lu), rep(1:lu, rep(length(d), lu)))
    if (any(du[1,]==0)) {
      emptynwstats <- rep(0, ncol(du))
      tmp <- du[2,du[1,]==0]
      for(i in 1:length(tmp)) tmp[i] <- sum(nodecov==tmp[i])
        emptynwstats[du[1,]==0] <- tmp
    }
  } else {
    if (any(d==0)) {
      emptynwstats <- rep(0, length(d))
      emptynwstats[d==0] <- network.size(nw)
    }
  }
  if(is.null(byarg)) {
    if(length(d)==0){return(NULL)}
    name <- "idegree"
    coef.names <- paste("idegree",d,sep="")
    inputs <- c(d)
  } else if (homophily) {
    if(length(d)==0){return(NULL)}
    name <- "idegree_w_homophily"
    coef.names <- paste("ideg", d, ".homophily.",byarg, sep="")
    inputs <- c(d, nodecov)
  } else {
    if(ncol(du)==0) {return(NULL)}
    #  No covariates here, so "ParamsBeforeCov" unnecessary
    name <- "idegree_by_attr"
    # See comment in d_idegree_by_attr function
    coef.names <- paste("ideg", du[1,], ".", byarg,u[du[2,]], sep="")
    inputs <- c(as.vector(du), nodecov)
  }

  c(.process_layers_degree(nw, a, name, coef.names, inputs, emptynwstats),
    minval = 0, maxval=network.size(nw), conflicts.constraints="idegreedist")
}

################################################################################

#' @templateVar name mutualL
#' @title Mutuality
#' @description Mutuality
#' @details In binary ERGMs, equal to the number of
#'   pairs of actors \eqn{i} and \eqn{j} for which \eqn{(i{\rightarrow}j)}{(i,j)}
#'   and \eqn{(j{\rightarrow}i)}{(j,i)} both exist. For valued ERGMs, equal to \eqn{\sum_{i<j} m(y_{i,j},y_{j,i})} ,
#'   where \eqn{m} is determined by `form` argument:
#'   - `"min"` for \eqn{\min(y_{i,j},y_{j,i})}
#'   - `"nabsdiff"` for \eqn{-|y_{i,j},y_{j,i}|}
#'   - `"product"` for \eqn{y_{i,j}y_{j,i}}
#'   - `"geometric"` for \eqn{\sqrt{y_{i,j}}\sqrt{y_{j,i}}}
#'
#'   See Krivitsky (2012) for a discussion of these statistics. `form="threshold"` simply
#'   computes the binary `mutuality` after thresholding at `threshold` .
#'
#'   This term can only be used with directed networks.
#'
#' @usage
#' # binary: mutualL(same=NULL, diff=FALSE, by=NULL, keep=NULL, Ls=NULL)
#' @param same optional argument. If passed the name of a vertex attribute,
#'   only mutual pairs that match on the attribute are counted. Only one of `same`
#'   or `by` may be used. If both parameters are passed, `same` takes precedent. This
#'   parameter is affected by `diff`.
#' @param diff separate counts for each unique matching value can be obtained by using
#'   `diff=TRUE` with `same`.
#' @param by each node is counted separately for each mutual pair in which it
#'   occurs and the counts are tabulated by unique values of the attribute if
#'   passed the name of a vertex attribute. This means that the sum of the mutual statistics when `by` is used
#'   will equal twice the standard mutual statistic. Only one of `same`
#'   or `by` may be used. If both parameters are passed, `same` takes precedent. This
#'   parameter is not affected by `diff`.
#' @param keep a numerical vector to specify which statistics should be kept whenever the `mutual` term would
#'   ordinarily result in multiple statistics.
#' @param Ls a list of two layer logic formulas or a single layer logic formula.
#'   If given, the statistic will count the number of dyads where a tie in
#'   `Ls[[1]]` reciprocates a tie in `Ls[[2]]` and vice
#'   versa. (Note that dyad that has mutual ties in `Ls[[1]]` and in
#'   `Ls[[2]]` will add 2 to this statistic.) If a formula is given,
#'   it is replicated.
#'
#' @template ergmTerm-general
#'
#' @concept directed
#' @concept frequently-used
#' @concept layer-aware
InitErgmTerm.mutualL<-function (nw, arglist, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, directed=TRUE, bipartite=NULL,
                      varnames = c("same", "by", "diff", "keep", "Ls"),
                      vartypes = c("character", "character", "logical", "numeric", "formula,list"),
                      defaultvalues = list(NULL, NULL, FALSE, NULL, NULL),
                      required = c(FALSE, FALSE, FALSE, FALSE, FALSE))
  ### Process the arguments
  if (!is.null(a$same) || !is.null(a$by)) {
    if (!is.null(a$same)) {
     attrname <- a$same
     if (!is.null(a$by)) 
       warning("Ignoring 'by' argument to mutual because 'same' exists", call.=FALSE)
    }else{
     attrname <- a$by
    }
    nodecov <- get.node.attr(nw, attrname)
    u <- NVL(a$levels, sort(unique(nodecov)))
    if (!is.null(a$keep)) {
      u <- u[a$keep]
    }
    #   Recode to numeric
    nodecov <- match(nodecov,u,nomatch=length(u)+1)
    # All of the "nomatch" should be given unique IDs so they never match:
    dontmatch <- nodecov==(length(u)+1)
    nodecov[dontmatch] <- length(u) + (1:sum(dontmatch))
    ui <- seq(along=u)
  }

  ### Construct the list to return
  if (!is.null(a$same) || !is.null(a$by)) {
    if (is.null(a$same)) {
      coef.names <- paste("mutual.by", attrname, u, sep=".")
      inputs <- c(ui, nodecov)
    }else{
     if (a$diff) {
      coef.names <- paste("mutual.same", attrname, u, sep=".")
      inputs <- c(ui, nodecov)
     }else{ 
      coef.names <- paste("mutual", attrname, sep=".")
      inputs <- nodecov
     }
    }
    if (is.null(a$same) && !is.null(a$by)) {
     name <- "mutual_by_attr"
    }else{
     name <- "mutual"
    }
  }else{
     name <- "mutual"
     coef.names <- "mutual"
     inputs <- NULL
  }

  maxval <- network.dyadcount(nw,FALSE)/2

  Ls <- .set_layer_namemap(a$Ls, nw)
  if(is(Ls,"formula")) Ls <- list(Ls)
  L1 <- Ls[[1]]
  L2 <- Ls[[2]]
  if(!is.null(L1) || !is.null(L2)){
    NVL(L1) <- L2
    NVL(L2) <- L1
    auxiliaries <- .mk_.layer.net_auxform(L1)
    aux2 <- .mk_.layer.net_auxform(L2)
    auxiliaries[[2]] <- call("+", auxiliaries[[2]], aux2[[2]])
    name <- paste(name, "ML", sep="_")
    coef.names <- .lspec_coef.namewrap(if(L1==L2) list(L1) else list(L1,L2))(coef.names)
    maxval <- maxval*2
  }else auxiliaries <- NULL
  
  list(name=name,                      #name: required
       coef.names = coef.names,        #coef.names: required
       inputs=inputs,
       auxiliaries = auxiliaries,
       minval = 0,
       maxval = maxval) 
}

################################################################################

#' @templateVar name odegreeL
#' @title Out-degree
#' @description Out-degree
#' @details This term adds one network statistic to
#'   the model for each element in `d` ; the \eqn{i} th such statistic equals
#'   the number of nodes in the network of out-degree `d[i]` , i.e. the
#'   number of nodes with exactly `d[i]` out-edges.
#'
#'   This term can only be used with directed networks; for undirected networks
#'   see `degree` .
#'
#'   If a list of layer specifications is given, degree of a node
#'   `i` is considered to be the number of edges in all layers,
#'   combined.
#'
#' @usage
#' # binary: odegreeL(d, by=NULL, homophily=FALSE, levels=NULL)
#' @param d a vector of distinct integers
#' @param by a character string giving the name of an attribute in the
#'   network's vertex attribute list.
#' @param homophily If this is specified and `homophily` is `TRUE` ,
#'   then degrees are calculated using the subnetwork consisting of only
#'   edges whose endpoints have the same value of the `by` attribute.
#'   If `by` is specified and
#'   `homophily` is `FALSE` (the default), then separate degree
#'   statistics are calculated for nodes having each separate
#'   value of the attribute.
#' @param levels list of layer specifications
#'
#' @template ergmTerm-general
#'
#' @concept directed
#' @concept categorical nodal attribute
#' @concept frequently-used
InitErgmTerm.odegreeL<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                      varnames = c("d", "by", "homophily", "levels", "Ls"),
                      vartypes = c("numeric", "character", "logical", "character,numeric,logical", "formula,list"),
                      defaultvalues = list(NULL, NULL, FALSE, NULL, NULL),
                      required = c(TRUE, FALSE, FALSE, FALSE, FALSE))
  d<-a$d; byarg <- a$by; homophily <- a$homophily
  emptynwstats<-NULL
  if(!is.null(byarg)) {
    nodecov <- get.node.attr(nw, byarg, "odegree")
    u <- NVL(a$levels, sort(unique(nodecov)))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u,nomatch=length(u)+1) # Recode to numeric
    if (length(u)==1)
         stop ("Attribute given to odegree() has only one value", call.=FALSE)
  }
  if(!is.null(byarg) && !homophily) {
    # Combine degree and u into 2xk matrix, where k=length(d)*length(u)
    lu <- length(u)
    du <- rbind(rep(d,lu), rep(1:lu, rep(length(d), lu)))
    if (any(du[1,]==0)) {
      emptynwstats <- rep(0, ncol(du))
      tmp <- du[2,du[1,]==0]
      for(i in 1:length(tmp)) tmp[i] <- sum(nodecov==tmp[i])
        emptynwstats[du[1,]==0] <- tmp
    }
  } else {
    if (any(d==0)) {
      emptynwstats <- rep(0, length(d))
      emptynwstats[d==0] <- network.size(nw)
    }
  }
  if(is.null(byarg)) {
    if(length(d)==0){return(NULL)}
    name <- "odegree"
    coef.names <- paste("odegree",d,sep="")
    inputs <- c(d)
  } else if (homophily) {
    if(length(d)==0){return(NULL)}
    name <- "odegree_w_homophily"
    # See comment in d_odegree_w_homophily function
    coef.names <- paste("odeg", d, ".homophily.",byarg, sep="")
    inputs <- c(d, nodecov)
  } else {
    if(ncol(du)==0) {return(NULL)}
    #  No covariates here, so "ParamsBeforeCov" unnecessary
    name <- "odegree_by_attr"
    # See comment in d_odegree_by_attr function
    coef.names <- paste("odeg", du[1,], ".", byarg,u[du[2,]], sep="")
    inputs <- c(as.vector(du), nodecov)
  }

  c(.process_layers_degree(nw, a, name, coef.names, inputs, emptynwstats),
    minval = 0, maxval=network.size(nw), conflicts.constraints="odegreedist")
}
