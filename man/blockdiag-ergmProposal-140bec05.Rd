% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/InitErgmProposal.blockdiag.R
\name{blockdiag-ergmProposal}
\alias{blockdiag-ergmProposal}
\alias{InitErgmProposal.blockdiag}
\title{A Metropolis--Hastings proposal for diagonal block constraints}
\description{
Typically used for \eqn{constraints= ~blockdiag}. Select a diagonal
block according to the weight, then randomly select a dayd within the
block for the toggle proposal.
}
\details{
\if{html}{\Sexpr[results=rd,stage=render]{ergm:::.formatProposalsHtml(ergm:::.buildProposalsList(proposal="blockdiag"))}}
\if{text}{\Sexpr[results=rd,stage=render]{ergm:::.formatProposalsText(ergm:::.buildProposalsList(proposal="blockdiag"))}}
\if{latex}{\Sexpr[results=rd,stage=render]{ergm:::.formatProposalsLatex(ergm:::.buildProposalsList(proposal="blockdiag"))}}
}
\seealso{
\code{\link[ergm:ergmProposal]{ergmProposal}} for index of proposals currently visible to the package.

\Sexpr[results=rd,stage=render]{ergm:::.formatTermKeywords("ergmProposal", "blockdiag", "subsection")}
}
\keyword{internal}
