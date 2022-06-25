#  File R/InitErgmTerm.multinet.R in package ergm.multi, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2022 Statnet Commons
################################################################################
#' A multinetwork network representation.
#'
#' A function for specifying the LHS of a multi-network (a.k.a. multilevel) ERGM.
#'
#' @param ... network specification, in one of two formats:
#' 
#'   1. An (optionally named) list of networks with same directedness and bipartedness (but possibly different sizes).
#'
#'   1. Several networks as (optionally named) arguments.
#'
#' @return A network object with multinetwork metadata.
#'
#' @seealso [Help on model specification][ergmTerm] for specific terms
#' @seealso `vignette("Goeyvaerts_reproduction")` for a demonstration
#' @examples
#'
#' data(samplk)
#'
#' # Method 1: list of networks
#' monks <- Networks(list(samplk1, samplk2))
#' ergm(monks ~ N(~edges))
#'
#' # Method 2: networks as arguments
#' monks <- Networks(samplk1, samplk2)
#' ergm(monks ~ N(~edges))
#' 
#' @export
Networks <- function(...){
  args <- list(...)
  if(all(sapply(args, is, "network"))){
    nwl <- args
  }else if(is.list(args[[1]]) && all(sapply(args[[1]], is, "network"))){
    nwl <- args[[1]]
  }else stop("Unrecognized format for multinetwork specification. See help for information.")

  if(!.same_constraints(nwl, "constraints")) stop("Networks have differing constraint structures. This is not supported at this time.")
  if(!.same_constraints(nwl, "obs.constraints")) stop("Networks have differing observation processes. This is not supported at this time.")
  
  nw <- combine_networks(nwl, blockID.vattr=".NetworkID", blockName.vattr=".NetworkName", ignore.nattr = c(eval(formals(combine_networks)$ignore.nattr), "constraints", "obs.constraints", "ergm"), subnet.cache=TRUE)

  nw %n% "ergm" <- .combine_ergmlhs(nwl)

  nw %ergmlhs% "constraints" <-
      if(NVL(nwl[[1]] %ergmlhs% "constraints",base_env(~.))==base_env(~.))
        base_env(~blockdiag(".NetworkID"))
      else
        append_rhs.formula(nwl[[1]]%ergmlhs%"constraints", list(call("blockdiag",".NetworkID")), TRUE)
  if(!is.null(nwl[[1]]%ergmlhs%"obs.constraints")) nw %ergmlhs% "obs.constraints" <- nwl[[1]]%ergmlhs%"obs.constraints"

  nw
}

InitErgmTerm..subnets <- function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("attrname"),
                      vartypes = c("character"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))

  list(name="_subnets", coef.names=c(), iinputs=c(unlist(.block_vertexmap(nw, a$attrname))), dependence=FALSE)
}

#' An `as_tibble` method for combined networks.
#'
#' A method to obtain a network attribute table from a
#' [`combined_networks`] object, falling back to the
#' [network::as_tibble.network()] if vertex or edge attributes are
#' required.
#'
#' @param x  a [`combined_networks`] (inheriting from [`network::network`]).
#' @param attrnames a list (or a selection index) for attributes to obtain; for combined networks, defaults to all.
#' @param unit whether to obtain edge, vertex, or network attributes.
#' @param ... additional arguments, currently passed to unlist()].
#'
#' @seealso [network::as_tibble.network()]
#' @export
as_tibble.combined_networks<-function(x,attrnames=(match.arg(unit)%in%c("vertices","networks")), ..., unit=c("edges", "vertices", "networks")){
  unit <- match.arg(unit)
  if(unit!="networks") return(NextMethod())

  al <-
    if(is(x, "combined_networks") && !is.null(x %n% ".subnetattr")) (x %n% ".subnetattr")[[".NetworkID"]]
    else{
      # FIXME: Probably more efficient to use attrnames earlier, in order to save calls to get.network.attribute().
      xl <- if(is.network(x)) .split_constr_network(x, ".NetworkID", ".NetworkName") else x
      nattrs <- Reduce(union, lapply(xl, list.network.attributes))
      lapply(nattrs, function(nattr) lapply(xl, get.network.attribute, nattr)) %>% set_names(nattrs)
    }

  if(is.logical(attrnames) || is.numeric(attrnames)) attrnames <- na.omit(names(al)[attrnames])
  else intersect(attrnames, names(al))

  al[attrnames] %>% lapply(simplify_simple, ...) %>% as_tibble()
}

get_lminfo <- function(nattrs, lm=~1, subset=TRUE, contrasts=NULL, offset=NULL, weights=1){
  nn <- nrow(nattrs)
  
  subset <-
    if(mode(subset) %in% c("expression", "call")) eval(if(is(subset, "formula")) subset[[2]] else subset, envir = nattrs, enclos = environment(lm))
    else subset
  subset <- unwhich(switch(mode(subset),
                           logical = which(rep(subset, length.out = nn)),
                           numeric = subset),
                    nn)

  weights <- if(mode(weights) %in% c("expression", "call")) eval(if(is(weights, "formula")) weights[[2]] else weights, envir = nattrs, enclos = environment(lm))
             else weights
  weights <- rep(weights, length.out=nn)

  offset <- if(mode(offset) %in% c("expression", "call")) eval(if(is(offset, "formula")) offset[[2]] else offset, envir = nattrs, enclos = environment(lm))
             else offset

  subset[weights==0] <- FALSE
  
  ## model.frame
  # This loop is necessary because if the predictor vector for a
  # network is all 0s as is its offest, it needs to be dropped.
  subset.prev <- NULL
  while(!identical(subset, subset.prev)){
    nm <- sum(subset)
    xf <- do.call(stats::lm, list(lm, data=nattrs,contrast.arg=contrasts, offset = offset, subset=subset, weights=weights, na.action=na.fail, method="model.frame"), envir=environment(lm))
    xm <- model.matrix(attr(xf, "terms"), xf, contrasts=contrasts)
    offset <- model.offset(xf) %>% NVL(numeric(nrow(xm))) %>% matrix(nrow=nrow(xm)) # offset needs to have reliable dimension.
    subset.prev <- subset
    subset[subset] <- subset[subset] & !apply(cbind(xm,offset)==0, 1, all)
  }

  list(xf=xf, xm=xm, subset=subset, offset=offset)  
}

#' @import purrr
#' @import tibble
#' @templateVar name N
#' @title Evaluation on multiple networks
#' @description Evaluates the terms in `formula` on each of the networks joined
#'   using [`Networks`] function, and returns either a weighted
#'   sum or an [`lm`] -style linear model for the ERGM
#'   coefficients. Its syntax follows that of [`lm`] closely,
#'   with sensible defaults. The default formula ( `~1` ) sums the
#'   specified network statistics. If `lm` refers to any
#'   network attributes for which some networks have missing values, the
#'   term will stop with an error. This can be avoided by pre-filtering
#'   with `subset`, which controls which networks are affected by
#'   the term.
#'
#' @note Care should be taken to avoid multicollinearity when using this operator. When `lm` is given a model with intercept and a categorical predictor (including a [`logical`] one), it will use the first level (or `FALSE` ) as the baseline, but if the model is without intercept, it will use all levels of the first categorical predictor. This is typically what is wanted in a linear regression, but for the `N` operator, this can be problematic if the "intercept" effect is added by a different term. A workaround is to convert the categorical predictor to dummy variables before putting it into the `lm` formula.
#'
#' @usage
#' # binary: N(formula, lm=~1, subset=TRUE, weights=1, contrasts=NULL, offset=0, label=NULL)
#' @template ergmTerm-formula
#' @param lm a one-sided [lm()]-style formula whose RHS specifies the network-level predictors for the terms in the [ergm()] formula `formula`.
#' @param subset,contrasts see [lm()].
#' @param offset A constant, a vector of length equal to the number of networks, or a matrix whose number of rows is the number of networks and whose number of columns is the number of free parameters of the ERGM. It can be specified in `lm` as well.
#' @param weights reserved for future use; attempting to change it will cause an
#'   error: at this time, there is no way to assign sampling weights to
#'   networks.
#' @param label An optional parameter which will add a label to model
#'   parameters to help identify the term (which may have similar
#'   predictors but, say, a different network subset) in the output.
#'
#'
#' @template ergmTerm-general
#'
#' @concept operator
#' @concept directed
#' @concept undirected
InitErgmTerm.N <- function(nw, arglist, N.compact_stats=TRUE,...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula","lm","subset","weights","contrasts","offset","label"),
                      vartypes = c("formula","formula","formula,logical,numeric,expression,call","formula,logical,numeric,expression,call","list","formula,logical,numeric,expression,call","character"),
                      defaultvalues = list(NULL,~1,TRUE,1,NULL,NULL,NULL),
                      required = c(TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE))

  auxiliaries <- base_env(~.subnets(".NetworkID"))

  nwl <- .split_constr_network(nw, ".NetworkID", ".NetworkName")
  nwnames <- names(nwl)
  nn <- length(nwl)
  nattrs <- as_tibble(nw, unit="networks")

  lmi <- get_lminfo(nattrs, lm=a$lm, subset=a$subset, contrasts=a$contrasts, offset=a$offset, weights=a$weights)

  xf <- lmi$xf
  xm <- lmi$xm
  subset <- lmi$subset
  weights <- lmi$weights
  offset <- lmi$offset
  nm <- sum(subset)
  nwl <- nwl[subset]
  rm(lmi)

  ms <- lapply(nwl, function(nw1){
    ergm_model(a$formula, nw1, ..., offset.decorate=FALSE)
  })

  nparams <- ms %>% map_int(nparam, canonical=FALSE)
  nstats <-  ms %>% map_int(nparam, canonical=TRUE)

  # Check for MANOVA style matrix.
  if(!all_identical(nparams)) ergm_Init_abort("N() operator only supports models with the same numbers of parameters for every network. This may change in the future.")
  if(!all_identical(ms %>% map(param_names, canonical=FALSE))) ergm_Init_warn("Subnetwork models have different parameter names but the same parameter vector lengths; this may indicate specification problems.")
  nparam <- nparams[1]

  # Extract offsets and weights.
  offset <- model.offset(xf) %>% NVL(numeric(nm)) %>% matrix(nrow=nm, ncol=nparam) # offset is actually an nm*q matrix.
  weights <- model.weights(xf)
  if(!all(weights==1)) ergm_Init_abort("Network-level weights different from 1 are not supported at this time.")

  # Model parameters.
  params <- rep(list(NULL), nparam*ncol(xm))
  parnames <- colnames(xm)
  parnames <- ifelse(parnames=="(Intercept)", "1", parnames)
  names(params) <- ergm_mk_std_op_namewrap('N',a$label)(rep(param_names(ms[[1]], canonical=FALSE), each=ncol(xm)),parnames)
  if(with(ms[[1]]$etamap,
          any(mintheta[!offsettheta]!=-Inf) || any(maxtheta[!offsettheta]!=+Inf))){
    warning("Submodel specified to N() operator with a linear model formula has parameter constraints. They will be ignored.")
  }
  
  nstats.all <- integer(nn)
  nstats.all[subset] <- nstats # So networks not in subset get 0 stats.

  # Get the extended state and empty network stats wrapping. Note that
  # naming, curved models, etc., are handled "in-house" (for now).
  wms <- mapply(wrap.ergm_model, ms, nwl, SIMPLIFY=FALSE, USE.NAMES=FALSE)
  
  ## An important special case is when all models are linear and have
  ## the same number of stats. lm-offset does not work with this
  ## approach at this time, either *unless* all terms are offset by
  ## the same amount (that may vary over networks). Also,
  ## N.compact_stats can be set to FALSE to force the general case to
  ## be used.
  # TODO: A more refined detection than is.curved(), since currently
  # all operators show up as curved.
  # TODO: Check that multivariate offset is still not allowed in R.
  if(N.compact_stats &&
     all_identical(nstats) &&
     !any(ms%>%map_lgl(is.curved)) &&
     all(apply(offset,1,all_identical))){
    xm.all <- matrix(0, nn, ncol(xm))
    xm.all[subset,] <- xm
    
    if(!all(offset==0)){
      offset.all <- numeric(nn)
      offset.all[subset] <- offset[,1]
    }else offset.all <- NULL

    iinputs <- ncol(xm)+!is.null(offset.all)
    inputs <- c(t(cbind(xm.all,offset.all)))

    if(is.null(offset.all)){
      # This can be represented as a fully linear term.
      coef.names <- names(params)
      params <- etamap <- etagradient <- NULL
    }else{
      # Offset requires a bit of extra work.
      etamap <- function(x, n, ...){
        x %>% matrix(ncol=nparam) %>% rbind(1) %>% c()
      }
      etagradm <- list(cbind(diag(1, ncol(xm)), 0)) %>%
        rep(nparam) %>%
        Matrix::.bdiag() %>%
        as.matrix()
      etagradient <- function(x, n, ...){
        etagradm
      }
      coef.names <- names(params) %>% matrix(ncol=nparam) %>% rbind(paste0("offset",seq_len(nparam))) %>% c()
    }
    
    gss <- map(wms, "emptynwstats")
    gs <-
      if(all(map_lgl(gss, is.null))){ # Linear combination of 0s is 0.
        NULL
      }else{ # All numeric or NULL
        gs0 <- map_lgl(gss, is.null)
        lst(X = (cbind(xm,offset.all[subset]) %>% split(., row(.)))[!gs0],
            Y = gss[!gs0]) %>%
          pmap(outer) %>%
          reduce(`+`) %>%
          c()
    }

    list(name="MultiNet", coef.names = coef.names, inputs=inputs, iinputs=iinputs, submodels=ms, dependence=any(map_lgl(wms, "dependence")), emptynwstats = gs, auxiliaries = auxiliaries, map = etamap, gradient = etagradient, params = params)
        
  }else{
    Xl <- xm %>%
      split(., row(.)) %>% # Rows of xl as a list.
      map(rbind) %>% # List of row vectors.
      map(list) %>% # List of singleton lists of row vectors.
      map(rep, nparam) %>% # List of lists with appropriate parameter numbers.
      map(Matrix::.bdiag) %>% # nm-list of block-diagonal matrices. 
      map(as.matrix) # FIXME: Maintain representation as sparse matrices?
    
    # Xl can now be used as covariates to vectorized MANOVA-style
    # parameters.
    
    ol <- offset %>% split(., row(.)) # Rows of offset as a list.
    
    iinputs <- c(0,cumsum(nstats.all))
    
    etamap <- function(x, n, ...){
      Xl %>%
        map(`%*%`, c(x)) %>% # Evaluate X%*%c(x), where X is the predictor matrix and x is "theta".
        map(c) %>% # Output of the previous step is a matrix; convert to vector.
        map2(ol, `+`) %>% # Add on offset.
        map2(ms %>% map("etamap"), # Obtain the etamap.
             ergm.eta) %>% # Evaluate eta.
        #      map2(weights, `*`) %>% # Weight.
        unlist()
    }
    etagradient <- function(x, n, ...){
      Xl %>%
        map(`%*%`, c(x)) %>% # Evaluate X%*%c(x), where X is the predictor matrix and x is "theta".
        map(c) %>% # Output of the previous step is a matrix; convert to vector.
        map2(ol, `+`) %>% # Add on offset.
        map2(ms %>% map("etamap"), # Obtain the etamap.
             ergm.etagrad) %>% # Evaluate eta gradient.
        map2(Xl, ., crossprod) %>%
        #      map2(weights, `*`) %>% # Weight.
        do.call(cbind,.)
    }

    coef.names <- ergm_mk_std_op_namewrap(paste0('N#',rep(seq_len(nm), nstats)))(unlist(lapply(ms, param_names, canonical=TRUE)))

    # Empty network statistics.
    gss <- map(wms, "emptynwstats")
    gs <-
      if(all(map_lgl(gss, is.null))){ # Linear combination of 0s is 0.
        NULL
      }else{ # All numeric or NULL
        unlist(mapply(function(gs, nst) if(is.null(gs)) numeric(nst) else gs, gss, nstats, SIMPLIFY=FALSE))
      }

    list(name="MultiNets", coef.names = coef.names, iinputs=iinputs, submodels=ms, dependence=any(map_lgl(wms, "dependence")), emptynwstats = gs, auxiliaries = auxiliaries, map = etamap, gradient = etagradient, params = params)
  }
}


#' @templateVar name ByNetDStats
#' @title Evaluate the formula on each individual network
#' @description Evaluates the terms in `formula` on each of the
#'   networks joined using [`Networks`] function, and returns a vector
#'   of its statistics for each network. This term is used internally
#'   for model diagnostics, and should not, generally, be used in models.
#'
#' @usage
#' # binary: N(formula, subset=TRUE)
#' @template ergmTerm-formula
#' @param subset see [lm()].
#'
#' @template ergmTerm-general
#'
#' @concept operator
#' @concept directed
#' @concept undirected
#' @noRd
InitErgmTerm.ByNetDStats <- function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula", "subset"),
                      vartypes = c("formula","formula,logical,numeric,expression,call"),
                      defaultvalues = list(NULL,TRUE),
                      required = c(TRUE,FALSE))

  auxiliaries <- base_env(~.subnets(".NetworkID"))
  nattrs <- as_tibble(nw, unit="networks")

  lmi <- get_lminfo(nattrs, subset=a$subset)

  subset <- lmi$subset
  nn <- sum(subset)
  rm(lmi)
    
  m <- ergm_model(a$formula, nw, ..., offset.decorate=FALSE)

  coef.names <- ergm_mk_std_op_namewrap(paste0('N#',rep(seq_len(nn)[subset], each=nparam(m, canonical=TRUE))))(param_names(m, canonical=TRUE))
  iinputs <- cumsum(c(-1,subset))*nparam(m, canonical=TRUE)
  wm <- wrap.ergm_model(m, nw, NULL) # Only need a few things from wrap.
  list(name="ByNetDStats", coef.names = coef.names, iinputs=iinputs, submodel=m, auxiliaries=auxiliaries, dependence=wm$dependence, ext.encode=wm$ext.encode)
}
