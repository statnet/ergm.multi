#  File R/InitErgmReference.R in package ergm.multi, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2021 Statnet Commons

#' @name Bernoulli-ergmReference
#' @title Bernoulli-reference ERGM
#' @description Bernoulli-reference ERGM
#' @details Specifies each
#'   dyad's baseline distribution to be Bernoulli with probability of
#'   the tie being \eqn{0.5} . This is the only reference measure used
#'   in binary mode.
#'
#' @usage
#' # Bernoulli
#'
#' @template ergmReference-general
#'
NULL

#' @name DiscUnif-ergmReference
#' @title Discrete-Uniform-reference ERGM
#' @description Discrete-Uniform-reference ERGM
#' @details Specifies each dyad's baseline distribution to be discrete uniform
#'   between `a` and `b` (both inclusive): \eqn{h(y)=1} , with
#'   the support being
#'   `a` , `a` +1, , `b` -1, `b` . At this time, both `a` and
#'   `b` must be finite.
#'
#' @usage
#' # DiscUnif(a,b)
#'
#' @template ergmReference-general
#'
NULL

#' @name Unif-ergmReference
#' @title Coninuous-Uniform-reference ERGM
#' @description Coninuous-Uniform-reference ERGM
#' @details Specifies each dyad's baseline distribution to be continuous uniform
#'   between `a` and `b` (both inclusive): \eqn{h(y)=1} , with
#'   the support being
#'   `[a, b]`. At this time, both `a` and
#'   `b` must be finite.
#'
#' @usage
#' # Unif(a,b)
#'
#' @template ergmReference-general
#'
NULL

#' @name StdNormal-ergmReference
#' @title Standard-Normal-reference ERGM
#' @description Standard-Normal-reference ERGM
#' @details Specifies each dyad's baseline distribution to be the normal distribution
#'   with mean 0 and variance 1.
#'
#' @usage
#' # StdNormal
#'
#' @template ergmReference-general
#'
NULL
