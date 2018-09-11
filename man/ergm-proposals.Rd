%  File man/ergm-proposals.Rd in package ergm, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2003-2017 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{ergm-proposals}
\alias{ergm-proposals}
\alias{InitErgmProposal}
\alias{InitErgmProposal.blockdiag}
\alias{InitErgmProposal.blockdiagNonObserved}
\alias{InitErgmProposal.blockdiagTNT}

\alias{InitErgmProposal.blockdiagNonObservedTNT}

\keyword{internal}
\docType{methods}

\title{Metropolis-Hastings Proposal Methods for ERGM MCMC}
\description{
  \code{\link{ergm}}  uses a Metropolis-Hastings (MH) algorithm to control the behavior of the Markov Chain Monte Carlo (MCMC) for sampling networks.  The MCMC chain is intended to step around the sample space of possible networks, selecting a network at regular intervals to evaluate the statistics in the model.  For each MCMC step, \eqn{n} (\eqn{n=1} in the simple case) toggles are proposed to change the dyad(s) to the opposite value. The probability of accepting the proposed change is determined by the MH acceptance ratio.  The role of the different MH methods implemented in \code{\link{ergm}} is to vary how the sets of dyads are selected for toggle proposals.  This is used in some cases to improve the performance (speed and mixing) of the algorithm, and in other cases to constrain the sample space.
}

\section{MH proposal methods implemented in the \code{\link[=ergm-package]{ergm}} package}{

\bold{MH proposals for non-constrained ergm models}

  \describe{
  
        \item{InitErgmProposal.randomtoggle}{Propose a randomly selected dyad to toggle.
        }
        
        \item{InitErgmProposal.TNT}{Default MH algorithm. Stratifies the population of dyads
        	by edge status: those having ties and those having no ties (hence T/NT).  
        	This is useful for improving performance in sparse networks,
        	because it gives at least 50\% chance of proposing a toggle of an existing edge.
        }
             
}

\bold{MH proposals for constrained ergm models}


    \describe{

         \item{InitErgmProposal.blockdiag}{MHp for \eqn{constraints= ~blockdiag}. Select a diagonal block according to the weight, then randomly select a dayd within the block for the toggle proposal.
        }
        \item{InitErgmProposal.blockdiagNonObserved}{MHp for \eqn{constraints= ~blockdiag+observed}. Similar to InitErgmProposal.blockdiag, but
        applied only to missing dyads.
        }
        
         \item{InitErgmProposal.blockdiagNonObservedTNT}{Similar to InitErgmProposal.blockdiagNonObserved, except that it selects ties and non-ties for proposed toggles (in the block by construction) with equal probability.  Like the unconstrained TNT proposal, this is useful for improving performance in sparse networks.
        }
        
        
        \item{InitErgmProposal.blockdiagTNT}{MHp for \eqn{constraints= ~blockdiag}. Similar to InitErgmProposal.blockdiag, except that it selects ties
        and non-ties for proposed toggles (in the block by construction) with equal probability.  Like the unconstrained TNT proposal, this is useful for improving performance in sparse networks.
        }

       }}

\seealso{\code{\link[=ergm-package]{ergm}} package, \code{\link{ergm}}, \code{\link{ergm-constraints}},  \code{\link{ergm_proposal}}}

  
\references{

Goodreau SM, Handcock MS, Hunter DR, Butts CT, Morris M (2008a).
 A \pkg{statnet} Tutorial.
 \emph{Journal of Statistical Software}, 24(8).
 \url{http://www.jstatsoft.org/v24/i08/}.

Hunter, D. R. and Handcock, M. S. (2006)
    \emph{Inference in curved exponential family models for networks},
   Journal of Computational and Graphical Statistics.

Hunter DR, Handcock MS, Butts CT, Goodreau SM, Morris M (2008b).
 \pkg{ergm}: A Package to Fit, Simulate and Diagnose
  Exponential-Family Models for Networks.
 \emph{Journal of Statistical Software}, 24(3).
 \url{http://www.jstatsoft.org/v24/i03/}.

Krivitsky PN (2012). Exponential-Family Random Graph Models for Valued
Networks. \emph{Electronic Journal of Statistics}, 2012, 6,
1100-1128. \doi{10.1214/12-EJS696}
 
Morris M, Handcock MS, Hunter DR (2008).
 Specification of Exponential-Family Random Graph Models:
  Terms and Computational Aspects.
 \emph{Journal of Statistical Software}, 24(4).
 \url{http://www.jstatsoft.org/v24/i04/}.
  
}

\keyword{models}



