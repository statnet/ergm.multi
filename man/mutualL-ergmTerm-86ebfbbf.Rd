% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/InitErgmTerm.multilayer.R
\name{mutualL-ergmTerm}
\alias{mutualL-ergmTerm}
\alias{InitErgmTerm.mutualL}
\title{Mutuality}
\usage{
# binary: mutualL(same=NULL, diff=FALSE, by=NULL, keep=NULL, Ls=NULL)
}
\arguments{
\item{same}{optional argument. If passed the name of a vertex attribute,
only mutual pairs that match on the attribute are counted. Only one of \code{same}
or \code{by} may be used. If both parameters are passed, \code{same} takes precedent. This
parameter is affected by \code{diff}.}

\item{diff}{separate counts for each unique matching value can be obtained by using
\code{diff=TRUE} with \code{same}.}

\item{by}{each node is counted separately for each mutual pair in which it
occurs and the counts are tabulated by unique values of the attribute if
passed the name of a vertex attribute. This means that the sum of the mutual statistics when \code{by} is used
will equal twice the standard mutual statistic. Only one of \code{same}
or \code{by} may be used. If both parameters are passed, \code{same} takes precedent. This
parameter is not affected by \code{diff}.}

\item{keep}{a numerical vector to specify which statistics should be kept whenever the \code{mutual} term would
ordinarily result in multiple statistics.}

\item{Ls}{a list (constructed by \code{\link[=list]{list()}} or \code{\link[=c]{c()}} of one or two Layer Logic specifications (c.f. Layer Logic section in the \code{\link[=Layer]{Layer()}} documentation). If given, the statistic will count the number of dyads where a tie in \code{Ls[[1]]} reciprocates a tie in \code{Ls[[2]]} and vice versa. (Note that dyad that has mutual ties in \code{Ls[[1]]} and in \code{Ls[[2]]} will add 2 to this statistic.) If a formula is given, it is replicated.}
}
\description{
In binary ERGMs, equal to the number of
pairs of actors \eqn{i} and \eqn{j} for which \eqn{(i{\rightarrow}j)}{(i,j)}
and \eqn{(j{\rightarrow}i)}{(j,i)} both exist.
}
\details{
This term can only be used with directed networks.
}
\seealso{
\code{\link[ergm:ergmTerm]{ergmTerm}} for index of model terms currently visible to the package.

\Sexpr[results=rd,stage=render]{ergm:::.formatTermKeywords("ergmTerm", "mutualL", "subsection")}
}
\concept{directed}
\concept{frequently-used}
\concept{layer-aware}
