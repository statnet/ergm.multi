%  File inst/NEWS.Rd in package ergm.multi, part of the
%  Statnet suite of packages for network analysis, https://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  https://statnet.org/attribution .
%
%  Copyright 2003-2024 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{NEWS}
\title{NEWS file for the \code{\link[=ergm-package]{ergm}} package}
\encoding{UTF-8}

% When a release is forked, copy the sections to a new release and
% delete unused sections.

%% \section{Changes in version X.Y.Z}{

%%   \subsection{NEW FEATURES}{
%%     \itemize{
%%       \item{
%%       }
%%     }
%%   }


%%   \subsection{BUG FIXES}{
%%     \itemize{
%%       \item{
%%       }
%%     }
%%   }


%%   \subsection{OTHER USER-VISIBLE CHANGES}{
%%     \itemize{
%%       \item{
%%       }
%%     }
%%   }

  %% \subsection{C-LEVEL FACILITIES}{
  %%   \itemize{
  %%   }
  %% }

  % \subsection{UTILITIES}{
  %   \itemize{
  %     \item Likewise, "UTILITIES" is probably unnecessary.
  %   }
  % }

  % \subsection{INSTALLATION}{
  %   \itemize{
  %     \item ... as is "INSTALLATION".
  %   }
  % }

  % \subsection{PACKAGE INSTALLATION}{
  %   \itemize{
  %     \item ditto.
  %   }
  % }

%% }


\section{Changes in version 0.2.1}{

  \subsection{NEW FEATURES}{
    \itemize{
      \item{
        \code{\link{marg_cond_sim}()} now has an argument for how to handle when the constrained sample has a higher variance than the unconstrained. (Previously, it would always throw an error.)
      }
    }
  }


  \subsection{BUG FIXES}{
    \itemize{
      \item{
        Fixes for documentation issues flagged by CRAN and compatibility with the latest \CRANpkg{roxygen2}.
      }
    }
  }


  \subsection{OTHER USER-VISIBLE CHANGES}{
    \itemize{
      \item{
        Some citations updated.
      }
      \item{
        When the optimized block-diagonal proposals encounter an initial network with off-block-diagonal edges, they no longer throw an error but rather use the \CRANpkg{ergm} 4.6.0+ proposal declination API.
      }
    }
  }

}



\section{Changes in version 0.2.0}{

  \subsection{NEW FEATURES}{
    \itemize{
      \item \code{\link{gofN}} objects now how \code{\link{augment}} and \code{\link{autoplot}} methods for Tidyverse and \CRANpkg{ggplot2} diagnostics.

      \item \code{\link{plot.gofN}()} and \code{\link{autoplot.gofN}()} are now more flexible in how they label extreme points via the \code{id.label=} argument.
    }
  }

  \subsection{BUG FIXES}{
    \itemize{
      \item Fixed how parallel processing settings are handled in \code{\link{gofN}()} and \code{\link{marg_cond_sim}()}.
    }
  }

}



\section{Changes in version 0.1.2}{
  \subsection{BUG FIXES}{
    \itemize{
      \item \code{\link{marg_cond_sim}()} now expects \code{\link{control.gofN}()} as its \code{control=} list.

      \item Some unit test tolerances have been loosed, making false positives less likely.

      \item Miscellaneous documentation fixes.
    }
  }

  \subsection{OTHER USER-VISIBLE CHANGES}{
    \itemize{
      \item Tweaks to the \code{\link{combine_networks}()} API to make it more robust, space-efficient, and flexible; in particular to make it possible for \CRANpkg{tergm} to depend on it rather than maintaining its own copy.
    }
  }
}


\section{Changes in version 0.1.1}{

  \subsection{NEW FEATURES}{
    \itemize{
      \item{
        In the \code{\link{Layer}(nw, attrlist)} specification, \code{attrlist=} can now be a named vector, in which case the layers are renamed accordingly.
      }
    }
  }


  \subsection{BUG FIXES}{
    \itemize{
      \item{
        Tolerance for some tests of the \code{\link{N}()} operator has been increased. This fixes a sporadic issue on M-series MacOS systems due to using a different linear algebra library.
      }
      \item{
        A bug in \code{\link{gofN}()} for valued ERGMs has been fixed.
      }
    }
  }


  \subsection{OTHER USER-VISIBLE CHANGES}{
    \itemize{
      \item{
        A number of optimisations in multilayer ERGM code.
      }
      \item{
        [Layer] Logic backend now uses integer arithmetic, for a significant speed improvement, at the cost of behaviour slightly inconsistent with \R's, such as \code{/} behaving as \code{\%/\%} (integer division).
      }
    }
  }
}


\section{Changes in version 0.1.0}{

  \subsection{NEW FEATURES}{
    \itemize{
      \item{
        Initial public release.
      }
    }
  }

}
