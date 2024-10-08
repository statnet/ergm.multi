% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/InitErgmTerm.R
\name{b1degreeL-ergmTerm}
\alias{b1degreeL-ergmTerm}
\alias{InitErgmTerm.b1degreeL}
\title{Degree for the first mode in a bipartite (aka two-mode) network}
\usage{
# binary: b1degreeL(d, by=NULL, levels=NULL, Ls=NULL)
}
\arguments{
\item{d}{a vector of distinct integers.}

\item{by}{a character string giving the name of an attribute in the network's vertex
attribute list. If this is specified
then each node's degree is tabulated only with other nodes having the same
value of the \code{by} attribute.}

\item{levels}{if \code{by} is specified, which levels to consider.}

\item{Ls}{either a Layer Logic specification formula (c.f. Layer Logic section in the \code{\link[=Layer]{Layer()}} documentation) or a list thereof (constructed by \code{\link[=list]{list()}} or \code{\link[=c]{c()}}). If given, degree of a node}
}
\description{
This term adds one network statistic to the model for
each element in \code{d} ; the \eqn{i} th such statistic equals the number of
nodes of degree \code{d[i]} in the first mode of a bipartite network, i.e.
with exactly \code{d[i]} edges. The first mode of a bipartite network object
is sometimes known as the "actor" mode.

This term can only be used with undirected bipartite networks.
}
\seealso{
\code{\link[ergm:ergmTerm]{ergmTerm}} for index of model terms currently visible to the package.

\Sexpr[results=rd,stage=render]{ergm:::.formatTermKeywords("ergmTerm", "b1degreeL", "subsection")}
}
\concept{bipartite}
\concept{categorical nodal attribute}
\concept{frequently-used}
\concept{undirected}
