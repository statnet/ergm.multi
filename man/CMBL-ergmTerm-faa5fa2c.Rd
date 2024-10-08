% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/InitErgmTerm.multilayer.R
\name{CMBL-ergmTerm}
\alias{CMBL-ergmTerm}
\alias{InitErgmTerm.CMBL}
\title{Conway--Maxwell-Binomial dependence among layers}
\usage{
# binary: CMBL(Ls=~.)
}
\arguments{
\item{Ls}{a list (constructed by \code{\link[=list]{list()}} or \code{\link[=c]{c()}} of at least two Layer Logic specifications (c.f. Layer Logic section in the \code{\link[=Layer]{Layer()}} documentation).}
}
\description{
Models marginal dependence layers within each dyad by imposing
a Conway--Maxwell-Binomial (CMB) distribution on the number of
layers in each dyad that have a tie.

The term adds one statistic to the model, equalling the sum over
all the dyads in the network of \eqn{\log\{E!(R-E)!/R!\}} , where
\eqn{E} is the number of layers in \code{Ls} with an edge in that
dyad and \eqn{R} being the total number of layers in \code{Ls} .
}
\details{
A positive coefficient induces positive dependence and a negative
one induces negative dependence.
}
\seealso{
\code{\link[ergm:ergmTerm]{ergmTerm}} for index of model terms currently visible to the package.

\Sexpr[results=rd,stage=render]{ergm:::.formatTermKeywords("ergmTerm", "CMBL", "subsection")}
}
\concept{directed}
\concept{layer-aware}
\concept{undirected}
