#  File R/data.R in package ergm.multi, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2022 Statnet Commons
################################################################################
#' A sample of within-household contact networks in Flanders and Brussels
#'
#' This is a list of 318 [`network`] objects derived from contact
#' diary data collected by by \insertCite{GoSa18h;textual}{ergm.multi}. The
#' study recruited households in Flanders and Brussels-Capital region
#' with at least one child 12 or under. The networks are symmetrized.
#'
#' @section Nonstandard Network Attributes: \describe{
#'
#' \item{`included`}{(logical) whether the network was included in
#' Goeyvaerts's analysis. (Two were excluded.)}
#'
#' \item{`weekday`}{(logical) whether the contact diary on which the
#' network is based was collected on a weekday, as opposed to
#' weekend.}
#'
#' }
#'
#' @section Nonstandard Vertex Attributes: \describe{
#'
#' \item{`age`}{(numeric) the household member's age.}
#'
#' \item{`gender`}{(character) the household member's gender (`"F"`/`"M"`).}
#'
#' \item{`role`}{(character) the household member's inferred role (`"Father"`/`"Mother"`/`"Child"`/`"Grandmother"`).}
#'
#' }
#'
#'
#' @usage data(Goeyvaerts)
#' @docType data
#' @section Licenses and Citation: When publishing results obtained
#'   using this data set, the original authors
#'   \insertCite{GoSa18h}{ergm.multi} should be cited, along with this
#'   \R package.
#' @references \insertAllCited{}
#'
#' @source The data were collected and by \insertCite{GoSa18h;textual}{ergm.multi} and
#'   curated by Pietro Coletti.
#'
#' @seealso `vignette("Goeyvaerts_reproduction")` for a vignette reproducing the Goeyvaerts analysis and performing diagnostics
#' @keywords datasets
"Goeyvaerts"


#' A network of advice, collaboration, and friendship in a law firm
#'
#' This dataset contains a [`network`] of relations of various types
#' among 71 lawyers (partners and associates) in a New England
#' (Northeastern US) corporate law firm referred to as \dQuote{SG&R}
#' collected 1988--1991 by \insertCite{La01c;textual}{ergm.multi}.
#'
#' All relations are directed.
#'
#' @section Nonstandard Vertex Attributes: \describe{
#'
#' \item{`age`}{(numeric) the lawyer's age.}
#'
#' \item{`gender`}{(character) the lawyer's gender (`"man"`/`"woman"`).}
#'
#' \item{`office`}{(character) in which of the firm's three offices the lawyer is based (`"Boston"`/`"Hartford"`/`"Providence"`).}
#'
#' \item{`practice`}{(character) which area of law the lawyer practices (`"corporate"`/`"litigation"`).}
#'
#' \item{`school`}{(character) from which law school the lawyer graduated (`"Harvard/Yale"`/`"UConnecticut"`/`"other"`).}
#'
#' \item{`seniority`}{(numeric) the lawyer's seniority rank in the firm (1 = high).}
#'
#' \item{`status`}{(character) the lawyer's status in the firm (`"associate"`/`"partner"`).}
#'
#' \item{`yrs_frm`}{(numeric) the number of years the lawyer has been with the firm.}
#'
#' }
#'
#' @section Nonstandard Edge Attributes: Each directed edge
#'   \eqn{i\rightarrow j}{i -> j} has the following attributes:
#'   \describe{
#'
#' \item{`advice`}{(logical) whether \eqn{i} has reported receiving advice from \eqn{j}. (Note that as defined, advice flows from head of the directed edge to the tail.)}
#'
#' \item{`coworker`}{(logical) whether \eqn{i} has reported receiving \eqn{j}'s assistance in preparing documents. (Note that as defined, assistance flows from head of the directed edge to the tail.)}
#'
#' \item{`friendship`}{(logical) whether \eqn{i} considers {j} a friend outside of work.}
#'
#' }
#'
#' @usage data(Lazega)
#' @docType data
#' @section Licenses and Citation: When publishing results obtained
#'   using this data set, the original author
#'   \insertCite{La01c}{ergm.multi} should be cited, along with this
#'   \R package.
#' @references \insertAllCited{}
#'
#' @source This version of the dataset was retrieved from the [RSiena
#'   web site](https://www.stats.ox.ac.uk/~snijders/siena/Lazega_lawyers_data.htm)
#'   and was compiled by Christopher Steven Marcum and Pavel
#'   N. Krivitsky for \insertCite{KrKo20e;textual}{ergm.multi}.
#'
#' @examples
#' \donttest{
#' data(Lazega)
#' # Construct a multilayer network for ergm(). (See `?Layer` for syntax.)
#' LLazega <- Layer(Lazega, c("advice", "coworker", "friendship"))
#' # Specify a layer logic model.
#' efit <- ergm(LLazega ~ L(~edges + mutual, ~advice) +
#'                        L(~edges + mutual, ~coworker) +
#'                        L(~edges + mutual, ~friendship) +
#'                        L(~edges + mutual, ~advice&coworker) +
#'                        L(~edges + mutual, ~advice&friendship) +
#'                        L(~edges + mutual, ~coworker&friendship))
#' summary(efit)
#' }
#' @keywords datasets
"Lazega"
