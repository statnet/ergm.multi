#  File R/ergm.multi-package.R in package ergm.multi, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2022 Statnet Commons
################################################################################
#' @title \pkg{ergm.multi}: Fit, Simulate and Diagnose Exponential-Family Models for Multiple or Multilayer Networks
#'
#' @description \pkg{ergm.multi} is a collection of extensions and
#'   utilities for package \CRANpkg{ergm} to facilitate modeling of
#'   multilayer and multi-network models. Some experimental support
#'   for multimode networks is also implemented.
#'
#' @section Multilayer network models: Also known as multiplex,
#'   multirelational, or multivariate networks, in a multilayer network a pair
#'   of actors can have multiple simultaneous relations of different
#'   types. For example, in the [Lazega] lawyer data set included with
#'   this package, each pair of lawyers in the firm can have an advice
#'   relationship, a coworking relationship, a friendship
#'   relationship, or any combination thereof. Application of ERGMs to
#'   multilayer networks has a long history
#'   \insertCite{PaWa99l,LaPa99m}{ergm.multi}, and a number of \R
#'   packages exist for analysing and estimating them.
#'
#'   \pkg{ergm.multi} implements the general approach of
#'   \insertCite{KrKo20e;textual}{ergm.multi} for specifying
#'   multilayer ERGMs, including Layer Logic and the various
#'   cross-layer specifications. Its features include:
#'
#'   \describe{
#'
#'   \item{seamless integration with [ergm()]:}{Multilayer
#'     specification is contained entirely in an [ergm()]-style formula and can be
#'     nested with any other [ergm()] terms, including dynamic and multi-network.}
#'
#'   \item{unlimited layers:}{The number of layers in the modeled
#'     network is limited only by computing power.}
#'
#'   \item{flexibility and simplicity:}{*Any* valid binary ERGM can be
#'     specified for any later or a logical combination of layers
#'     using simple term operators.}
#'
#'   \item{heterogeneous layers:}{A network can have directed and
#'     undirected layers, which can be modeled jointly.}
#'
#'   \item{multimode/multilevel support (experimental):}{With some
#'     care, it is possible to specify models for unipartite and
#'     bipartie layers over different subsets of actors, which can be
#'     used to specify multimode models.}
#'
#'   }
#'
#'   See [Layer()] and [`ergmTerm?L`][L-ergmTerm] for examples.
#'
#'
#' @section Multi-network models: Joint modeling of independent
#'   samples of networks on disjoint sets of actors have a long
#'   history as well \insertCite{@ZiVa06m, @SlKo16m, @StSc19m, and @VeSl21e, for example}{ergm.multi}.
#'   \pkg{ergm.multi} facilitates
#'   fixed-effect models for samples of networks (possibly
#'   heterogeneous in size and composition), using a multivariate
#'   linear model for each network's ERGM parameters, with
#'   network-level attributes serving as predictors, as formulated by
#'   \insertCite{SlKo16m;textual}{ergm.multi} and
#'   \insertCite{KrCo22t;textual}{ergm.multi}.
#'
#'   Its features include:
#'
#'   \describe{
#'
#'   \item{seamless integration with [ergm()]:}{Multi-network model
#'     specification is contained entirely in an [ergm()]-style formula and can be
#'     nested with any other [ergm()] terms, including dynamic and multilayer.}
#'
#'   \item{flexibility and simplicity:}{*Any* valid binary or valued
#'     ERGM can be specified for the networks, using simple term
#'     operators and the network-level specification with an
#'     [lm()]-style formula.}
#'
#'   }
#'
#'   See [Networks()], [`ergmTerm?N`][N-ergmTerm] for specification,
#'   [gofN()] for diagnostic facilities, and
#'   `vignette("Goeyvaerts_reproduction")` for a demonstration.
#'
#' @name ergm.multi-package
#' @docType package
#' @author Pavel N. Krivitsky \email{pavel@@statnet.org}
#'
#' @references \insertAllCited{}
#' @keywords package models
NULL
