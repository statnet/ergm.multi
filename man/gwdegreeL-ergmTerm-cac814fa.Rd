% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/InitErgmTerm.R
\name{gwdegreeL-ergmTerm}
\alias{gwdegreeL-ergmTerm}
\alias{InitErgmTerm.gwdegreeL}
\title{Geometrically weighted degree distribution}
\usage{
# binary: gwdegreeL(decay, fixed=FALSE, attrname=NULL, cutoff=30, levels=NULL)
}
\arguments{
\item{decay}{non-negative model parameter that is the same as theta_s in
equation (14) in Hunter (2007).}

\item{fixed}{specify if the value supplied for \code{decay} may be fixed (if \code{fixed=TRUE} ),
or it may be used as merely the starting value for the estimation
in a curved exponential family model (the default).}

\item{attrname}{if specified, then separate degree
statistics are calculated for nodes having each separate
value of the attribute.}

\item{cutoff}{only relevant if \code{fixed=FALSE} . In that case it only uses this
number of terms in computing the statistics to reduce the computational
burden. Its default value can also be controlled by the \code{gw.cutoff} term option control parameter. (See \code{\link[ergm:control.ergm]{control.ergm}} .)}

\item{levels}{a list of layer specifications. If given, degree of a node
\code{i} is considered to be the number of edges in all layers,
combined.}
}
\description{
This term adds one network statistic to the model equal to the weighted
degree distribution with decay controlled by the \code{decay} parameter.

This term can only be used with undirected networks.
}
\seealso{
\code{\link[ergm:ergmTerm]{ergmTerm}} for index of model terms currently visible to the package.

\Sexpr[results=rd,stage=render]{ergm:::.formatTermKeywords("ergmTerm", "gwdegreeL", "subsection")}
}
\concept{curved}
\concept{frequently-used}
\concept{undirected}
