#  File R/InitErgmConstraint.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################

#' Terms used in Exponential Family Random Graph Models
#'
#' @name ergmTerm
#' @aliases ergm-terms ergm.terms terms-ergm terms.ergm InitErgmTerm InitErgmWtTerm
#' @docType package
#' @description The function [`ergm`] is used to fit exponential random graph
#' models, in which the probability of a given network, \eqn{y}, on a set of
#' nodes is \deqn{h(y) \exp\{\eta(\theta) \cdot }{h(y) exp{eta(theta).g(y)} /
#' c(theta),}\deqn{ g(y)\}/c(\theta)}{h(y) exp{eta(theta).g(y)} / c(theta),}
#' where \eqn{h(y)} is the reference measure (for valued network models),
#' \eqn{g(y)} is a vector of network statistics for \eqn{y}, \eqn{\eta(\theta)}
#' is a natural parameter vector of the same length (with
#' \eqn{\eta(\theta)\equiv\theta}{eta(theta)=theta} for most terms),
#' \eqn{\cdot}{.} is the dot product, and \eqn{c(\theta)} is the normalizing
#' constant for the distribution.
#'
#' The network statistics \eqn{g(y)} are entered as terms in the function call
#' to [`ergm`].  This page describes the possible terms (and hence
#' network statistics) included in [`ergm`][ergm-package] package.
#'
#' A cross-referenced HTML version of the term documentation is available via
#' `vignette('ergm-term-crossRef')` and terms can also be searched via
#' [`search.ergmTerms`].
#'
#' @section Multilayer networks:
#' In order to fit a model for multilayer
#' networks, first use [`Layer`] construct an LHS network that
#' [`ergm`] will understand as multilayered.
#'
#' Used on the formula directly, most, but not all, terms included in this
#' package will sum their statistics over the layers.
#'
#' Some terms are *layer-aware*. These terms have explicit (usually
#' optional) layer specification arguments. For other terms, an operator
#' `L(formula, Ls=~.)` can be used to specify to which layer(s) they
#' apply. Layer specification documentation follows.
#'
#' ## Layer Logic
#'
#' Each formula's right-hand side describes an observed layer *or* some
#' "logical" layer, whose ties are a function of corresponding ties in
#' observed layers.
#'
#' The observed layers can be referenced either by name or by number (i.e.,
#' order in which they were passed to [`Layer`]). When referencing by
#' number, enclose the number in quotation marks (e.g., "1") or
#' backticks (e.g., \dQuote{`1`}).
#'
#' \link[base:Arithmetic]{Arithmetical}, \link[base:Comparison]{relational},
#' and \link[base:Logic]{logical} operators can be used to combine them. All
#' listed operators are implemented, as well as functions [`abs`],
#' [`round`], and [`sign`]. Standard
#' \link[base:Syntax]{operator precedence} applies, so use of parentheses is
#' recommended to ensure the logical expression is what it looks like.
#'
#' For example, if LHS is \code{Layer(A=nwA, B=nwB)}, both \code{~`2`} and
#' \code{~B} refer to \code{nwB}, while \code{A&!B} refers to a
#' \dQuote{logical} layer that has ties that are in \code{nwA} but not in
#' \code{nwB}.
#'
#' Transpose function \code{\link{t}} applied to a directed layer will reverse
#' the direction of all relations (transposing the sociomatrix). Unlike the
#' others, it can only be used on an observed layer directly. For example,
#' \code{~t(`1`)&t(`2`)} is valid but \code{~t(`1`&`2`)} is not.
#'
#' At this time, logical expressions that produce complete graphs from empty
#' graph inputs (e.g., \code{A==B} or \code{!A}) are not supported.
#'
#' ## Summing layers
#'
#' Some terms combine multiple layers. I that case al ist of formulas may be
#' passed. For example, \code{Layer(nw1,nw2)~L(~edges, c(~`1`,~(`2`&!`1`)))}
#' produces the number of edges in layer 1 plus the number of edges in layer 2
#' but not in layer 1.
#'
#' For these formulas, one can specify the layer's weight on its left-handside.
#' For example, \code{Layer(nw1,nw2)~L(~edges, c(3~`1`,-1~(`2`&!`1`)))} will
#' produce three times the number of edges in layer 1, minus the number of
#' edges in layer 2 but not in layer 1.
#'
#' @section Nodal attribute levels:
#' Terms taking a categorical nodal covariate also take `levels`
#' argument. This can be used to control the set and the ordering of
#' attribute levels.
#'
#' @section Terms included in the [`ergm.multi`][ergm-package] package:
#' As noted above, a cross-referenced HTML version of the term documentation is
#' also available via `vignette('ergm-term-crossRef')` and terms
#' can also be searched via [`search.ergmTerms`].
#'
#' ## Term index
#' \if{latex}{\Sexpr[results=rd,stage=render]{ergm:::.formatIndexLatex(ergm:::.buildTermsDataframe("ergmTerm"))}}
#' \if{text}{\Sexpr[results=rd,stage=render]{ergm:::.formatIndexText(ergm:::.buildTermsDataframe("ergmTerm"))}}
#' \if{html}{\Sexpr[results=rd,stage=render]{ergm:::.formatIndexHtml(ergm:::.buildTermsDataframe("ergmTerm"))}}
#'
#' ## Frequently-used terms
#' \if{latex}{\Sexpr[results=rd,stage=render]{ergm:::.formatMatrixLatex(ergm:::.termMatrix("ergmTerm", categories=ergm:::FREQUENTLY_USED_TERM_CATEGORIES, only.include='frequently-used'))}}
#' \if{text}{\Sexpr[results=rd,stage=render]{ergm:::.formatMatrixText(ergm:::.termMatrix("ergmTerm", categories=ergm:::FREQUENTLY_USED_TERM_CATEGORIES, only.include='frequently-used'))}}
#' \if{html}{\Sexpr[results=rd,stage=render]{ergm:::.formatMatrixHtml(ergm:::.termMatrix("ergmTerm", categories=ergm:::FREQUENTLY_USED_TERM_CATEGORIES, only.include='frequently-used'))}}
#'
#' ## Operator terms
#' \if{latex}{\Sexpr[results=rd,stage=render]{ergm:::.formatMatrixLatex(ergm:::.termMatrix("ergmTerm", categories=ergm:::OPERATOR_CATEGORIES, only.include='operator'))}}
#' \if{text}{\Sexpr[results=rd,stage=render]{ergm:::.formatMatrixText(ergm:::.termMatrix("ergmTerm", categories=ergm:::OPERATOR_CATEGORIES, only.include='operator'))}}
#' \if{html}{\Sexpr[results=rd,stage=render]{ergm:::.formatMatrixHtml(ergm:::.termMatrix("ergmTerm", categories=ergm:::OPERATOR_CATEGORIES, only.include='operator'))}}
#'
#' ## All terms
#' \if{latex}{\Sexpr[results=rd,stage=render]{ergm:::.formatMatrixLatex(ergm:::.termMatrix("ergmTerm"))}}
#' \if{text}{\Sexpr[results=rd,stage=render]{ergm:::.formatMatrixText(ergm:::.termMatrix("ergmTerm"))}}
#' \if{html}{\Sexpr[results=rd,stage=render]{ergm:::.formatMatrixHtml(ergm:::.termMatrix("ergmTerm"))}}
#'
#' ## Terms by concept
#' \if{latex}{\Sexpr[results=rd,stage=render]{ergm:::.formatTocLatex(ergm:::.termToc("ergmTerm"))}}
#' \if{text}{\Sexpr[results=rd,stage=render]{ergm:::.formatTocText(ergm:::.termToc("ergmTerm"))}}
#' \if{html}{\Sexpr[results=rd,stage=render]{ergm:::.formatTocHtml(ergm:::.termToc("ergmTerm"))}}
#'
#' @seealso [`ergm`][ergm-package] package, [`search.ergmTerms`], [`ergm`], [`network`], [`%v%`], [`%n%`]
#'
#' @references
#' - Butts, CT.  (2008).  "A Relational Event Framework for Social
#' Action." *Sociological Methodology,* 38(1).
#'
#' - Davis, J.A. and Leinhardt, S.  (1972).  The Structure of Positive
#' Interpersonal Relations in Small Groups.  In J. Berger (Ed.),
#' *Sociological Theories in Progress, Volume 2*, 218--251.  Boston:
#' Houghton Mifflin.
#'
#' - Holland, P. W. and S. Leinhardt (1981). An exponential family of
#' probability distributions for directed graphs.  *Journal of the
#' American Statistical Association*, 76: 33--50.
#'
#' - Hunter, D. R. and M. S. Handcock (2006). Inference in curved
#' exponential family models for networks. *Journal of Computational and
#' Graphical Statistics*, 15: 565--583.
#'
#' - Hunter, D. R. (2007). Curved exponential family models for social
#' networks. *Social Networks*, 29: 216--230.
#'
#' - Krackhardt, D. and Handcock, M. S. (2007).  Heider versus Simmel:
#' Emergent Features in Dynamic Structures. *Lecture Notes in Computer
#' Science*, 4503, 14--27.
#'
#' - Krivitsky P. N. (2012). Exponential-Family Random Graph Models for
#' Valued Networks. *Electronic Journal of Statistics*, 2012, 6,
#' 1100-1128. \doi{10.1214/12-EJS696}
#'
#' - Robins, G; Pattison, P; and Wang, P.  (2009).  "Closure,
#' Connectivity, and Degree Distributions: Exponential Random Graph (p*) Models
#' for Directed Social Networks." *Social Networks,* 31:105-117.
#'
#' - Snijders T. A. B., G. G. van de Bunt, and C. E. G. Steglich.
#' Introduction to Stochastic Actor-Based Models for Network Dynamics.
#' *Social Networks*, 2010, 32(1), 44-60. \doi{10.1016/j.socnet.2009.02.004}
#'
#' - Morris M, Handcock MS, and Hunter DR. Specification of
#' Exponential-Family Random Graph Models: Terms and Computational Aspects.
#' *Journal of Statistical Software*, 2008, 24(4), 1-24.
#' \url{https://www.jstatsoft.org/v24/i04}
#'
#' - Snijders, T. A. B., P. E. Pattison, G. L. Robins, and M. S. Handcock
#' (2006). New specifications for exponential random graph models,
#' *Sociological Methodology*, 36(1): 99-153.
#' Sample Space Constraints for Exponential-Family Random Graph Models
#'
#' @keywords models
#'
#' @examples
#' \dontrun{
#' ergm(flomarriage ~ kstar(1:2) + absdiff("wealth") + triangle)
#'
#' ergm(molecule ~ edges + kstar(2:3) + triangle
#'                       + nodematch("atomic type",diff=TRUE)
#'                       + triangle + absdiff("atomic type"))
#' }
NULL

#' Reference Measures for Exponential-Family Random Graph Models
#'
#' @name ergmReference
#' @aliases ergm-references references-ergm ergm.references references.ergm
#' @docType package
#' @description This page describes the possible reference measures (baseline distributions)
#' for found in the [`ergm`][ergm-package] package, particularly the
#' default (Bernoulli) reference measure for binary ERGMs.
#'
#' The reference measure is specified on the RHS of a one-sided formula passed
#' as the `reference` argument to [`ergm`].  See the
#' [`ergm`] documentation for a complete description of how
#' reference measures are specified.
#'
#' @section Possible reference measures to represent baseline distributions in the [`ergm.multi`][ergm-package] package:
#' \if{latex}{\Sexpr[results=rd,stage=render]{ergm:::.formatIndexLatex(ergm:::.buildTermsDataframe("ergmReference"))}}
#' \if{text}{\Sexpr[results=rd,stage=render]{ergm:::.formatIndexText(ergm:::.buildTermsDataframe("ergmReference"))}}
#' \if{html}{\Sexpr[results=rd,stage=render]{ergm:::.formatIndexHtml(ergm:::.buildTermsDataframe("ergmReference"))}}
#'
#' ## All references
#' \if{latex}{\Sexpr[results=rd,stage=render]{ergm:::.formatMatrixLatex(ergm:::.termMatrix("ergmReference"))}}
#' \if{text}{\Sexpr[results=rd,stage=render]{ergm:::.formatMatrixText(ergm:::.termMatrix("ergmReference"))}}
#' \if{html}{\Sexpr[results=rd,stage=render]{ergm:::.formatMatrixHtml(ergm:::.termMatrix("ergmReference"))}}
#'
#' ## References by concept
#' \if{latex}{\Sexpr[results=rd,stage=render]{ergm:::.formatTocLatex(ergm:::.termToc("ergmReference"))}}
#' \if{text}{\Sexpr[results=rd,stage=render]{ergm:::.formatTocText(ergm:::.termToc("ergmReference"))}}
#' \if{html}{\Sexpr[results=rd,stage=render]{ergm:::.formatTocHtml(ergm:::.termToc("ergmReference"))}}
#'
#' @seealso [`ergm`][ergm-package], [`network`], `sna`, [`summary.ergm`], [`print.ergm`], `\%v\%`, `\%n\%`
#'
#' @references
#' - Hunter DR, Handcock MS, Butts CT, Goodreau SM, Morris M (2008b). \pkg{ergm}:
#' A Package to Fit, Simulate and Diagnose Exponential-Family Models for
#' Networks. *Journal of Statistical Software*, 24(3).
#' \url{https://www.jstatsoft.org/v24/i03/}.
#'
#' - Krivitsky PN (2012). Exponential-Family Random Graph Models for Valued
#' Networks. *Electronic Journal of Statistics*, 2012, 6, 1100-1128.
#' \doi{10.1214/12-EJS696}
#'
#' @keywords models
NULL

#' Constraints used in Exponential Family Random Graph Models
#'
#' @name ergmConstraint
#' @aliases ergm-constraints constraints-ergm ergm.constraints constraints.ergm ergm-hints hints
#' @docType package
#' @description [`ergm`] is used to fit exponential-family random graph models
#' (ERGMs), in which the probability of a given network, \eqn{y}, on a set of
#' nodes is \eqn{h(y) \exp\{\eta(\theta) \cdot g(y)\}/c(\theta)}, where
#' \eqn{h(y)} is the reference measure (usually \eqn{h(y)=1}), \eqn{g(y)} is a
#' vector of network statistics for \eqn{y}, \eqn{\eta(\theta)} is a natural
#' parameter vector of the same length (with \eqn{\eta(\theta)=\theta} for most
#' terms), and \eqn{c(\theta)} is the normalizing constant for the
#' distribution.
#'
#' This page describes the constraints (the networks \eqn{y} for which
#' \eqn{h(y)>0}) that are included with the [`ergm`][ergm-package]
#' package. Other packages may add new constraints.
#'
#' @section Constraints implemented in the [`ergm.multi`][ergm-package] package:
#' \if{latex}{\Sexpr[results=rd,stage=render]{ergm:::.formatIndexLatex(ergm:::.buildTermsDataframe("ergmConstraint"))}}
#' \if{text}{\Sexpr[results=rd,stage=render]{ergm:::.formatIndexText(ergm:::.buildTermsDataframe("ergmConstraint"))}}
#' \if{html}{\Sexpr[results=rd,stage=render]{ergm:::.formatIndexHtml(ergm:::.buildTermsDataframe("ergmConstraint"))}}
#'
#' ## All constraints
#' \if{latex}{\Sexpr[results=rd,stage=render]{ergm:::.formatMatrixLatex(ergm:::.termMatrix("ergmConstraint"))}}
#' \if{text}{\Sexpr[results=rd,stage=render]{ergm:::.formatMatrixText(ergm:::.termMatrix("ergmConstraint"))}}
#' \if{html}{\Sexpr[results=rd,stage=render]{ergm:::.formatMatrixHtml(ergm:::.termMatrix("ergmConstraint"))}}
#'
#' ## Constraints by concept
#' \if{latex}{\Sexpr[results=rd,stage=render]{ergm:::.formatTocLatex(ergm:::.termToc("ergmConstraint"))}}
#' \if{text}{\Sexpr[results=rd,stage=render]{ergm:::.formatTocText(ergm:::.termToc("ergmConstraint"))}}
#' \if{html}{\Sexpr[results=rd,stage=render]{ergm:::.formatTocHtml(ergm:::.termToc("ergmConstraint"))}}
#'
#' @references
#' - Goodreau SM, Handcock MS, Hunter DR, Butts CT, Morris M (2008a).  A
#' \pkg{statnet} Tutorial. *Journal of Statistical Software*, 24(8).
#' \url{https://www.jstatsoft.org/v24/i08/}.
#'
#' - Hunter, D. R. and Handcock, M. S. (2006) *Inference in curved
#' exponential family models for networks*, Journal of Computational and
#' Graphical Statistics.
#'
#' - Hunter DR, Handcock MS, Butts CT, Goodreau SM, Morris M (2008b).
#' \pkg{ergm}: A Package to Fit, Simulate and Diagnose Exponential-Family
#' Models for Networks.  *Journal of Statistical Software*, 24(3).
#' \url{https://www.jstatsoft.org/v24/i03/}.
#'
#' - Karwa V, Krivitsky PN, and Slavkovi\'c AB (2016). Sharing Social Network
#' Data: Differentially Private Estimation of Exponential-Family Random Graph
#' Models. *Journal of the Royal Statistical Society, Series C*, 66(3):
#' 481-500. \doi{10.1111/rssc.12185}
#'
#' - Krivitsky PN (2012). Exponential-Family Random Graph Models for Valued
#' Networks. *Electronic Journal of Statistics*, 6, 1100-1128.
#' \doi{10.1214/12-EJS696}
#'
#' - Morris M, Handcock MS, Hunter DR (2008).  Specification of
#' Exponential-Family Random Graph Models: Terms and Computational Aspects.
#' *Journal of Statistical Software*, 24(4). \url{https://www.jstatsoft.org/v24/i04/}.
#' @keywords models
NULL
