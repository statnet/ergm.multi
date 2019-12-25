\name{ergm-terms}
\alias{ergm-terms}
\alias{ergm.terms}
\alias{terms-ergm}
\alias{terms.ergm}
\alias{InitErgmTerm}

\alias{L}
\alias{N}
\alias{b1degreeL}
\alias{b1degrangeL}
\alias{b2degreeL}
\alias{b2degrangeL}
\alias{degreeL}
\alias{despL}
\alias{dgwdspL}
\alias{dgwespL}
\alias{dgwnspL}
\alias{dnspL}
\alias{gwb1degreeL}
\alias{gwb2degreeL}
\alias{gwdegreeL}
\alias{gwidegreeL}
\alias{gwdegreeL}
\alias{gwodegreeL}
\alias{idegreeL}
\alias{degreeL}
\alias{CMBL}
\alias{mutualL}
\alias{twostarL}
\docType{methods}
\title{Terms used in Exponential Family Random Graph Models}
\description{
  The function \code{\link{ergm}} is used to fit exponential 
  random graph models, in which
  the probability of a given network, \eqn{y}, on a set of nodes is 
  \deqn{h(y) \exp\{\eta(\theta) \cdot
    g(y)\}/c(\theta)}{h(y) exp\{eta(theta).g(y)\} / c(theta),} where
  \eqn{h(y)} is the reference measure (for valued network models),
  \eqn{g(y)} is a vector of network statistics for \eqn{y},
  \eqn{\eta(\theta)} is a natural parameter vector of the same 
  length (with \eqn{\eta(\theta)\equiv\theta}{eta(theta)=theta} for most
  terms), \eqn{\cdot}{.} is the dot product, and \eqn{c(\theta)} is the
  normalizing constant for the distribution.

  The network statistics \eqn{g(y)} are entered as terms in the 
  function call to \code{\link{ergm}}.

  This page describes the possible terms (and hence network statistics)
  included in \code{\link[=ergm-package]{ergm}} package. Other packages
  may add their own terms, and package
  \code{\link[ergm.userterms:ergm.userterms-package]{ergm.userterms}} provides tools
  for implementing them.

  The current recommendation for any package implementing additional
  terms is to create a help file with a name or alias \code{ergm-terms},
  so that \code{help("ergm-terms")} will list ERGM terms available from
  all loaded packages.
}

\section{Multilayer networks}{
  In order to fit a model for multilayer networks, first use
  \code{\link{Layer}} construct an LHS network that \code{\link{ergm}}
  will understand as multilayered.

  Used on the formula directly, most, but not all, terms included in
  this package will sum their statistics over the layers.

  Some terms are \emph{layer-aware}. These terms have explicit (usually
  optional) layer specification arguments. For other terms, an operator
  \code{L(formula, Ls=~.)} can be used to specify to which layer(s) they
  apply. Layer specification documentation follows.

  \subsection{Layer Logic}{

    Each formula's right-hand side describes an observed layer \emph{or} some
    \dQuote{logical} layer, whose ties are a function of corresponding
    ties in observed layers.

    The observed layers can be referenced either by name or by number
    (i.e., order in which they were passed to \code{\link{Layer}}). When
    referencing by number, enclose the number in quotation marks (e.g.,
    \dQuote{"1"}) or backticks (e.g., \dQuote{`1`}).

    \link[base:Arithmetic]{Arithmetical},
    \link[base:Comparison]{relational}, and \link[base:Logic]{logical}
    operators can be used to combine them. All listed operators are
    implemented, as well as functions \code{\link{abs}},
    \code{\link{round}}, and \code{\link{sign}}. Standard \R
    \link[base:Syntax]{operator precedence} applies, so use of
    parentheses is recommended to ensure the logical expression is what
    it looks like.

    For example, if LHS is \code{Layer(A=nwA, B=nwB)}, both \code{~`2`}
    and \code{~B} refer to \code{nwB}, while \code{A&!B} refers to a
    \dQuote{logical} layer that has ties that are in \code{nwA} but not
    in \code{nwB}.

    Transpose function \code{\link{t}} applied to a directed layer will
    reverse the direction of all relations (transposing the
    sociomatrix). Unlike the others, it can only be used on
    an observed layer directly. For example, \code{~t(`1`)&t(`2`)} is
    valid but \code{~t(`1`&`2`)} is not.

    At this time, logical expressions that produce complete graphs from
    empty graph inputs (e.g., \code{A==B} or \code{!A}) are not
    supported.

  }

  \subsection{Summing layers}{

    Some terms combine multiple layers. I that case al ist of formulas
    may be passed. For example, \code{Layer(nw1,nw2)~L(~edges,
    c(~`1`,~(`2`&!`1`)))} produces the number of edges in layer 1 plus the
    number of edges in layer 2 but not in layer 1.

    For these formulas, one can specify the layer's weight on its
    left-handside. For example, \code{Layer(nw1,nw2)~L(~edges,
    c(3~`1`,-1~(`2`&!`1`)))} will produce three times the number of
    edges in layer 1, minus the number of edges in layer 2 but not in
    layer 1.
    
  }
}

\section{Nodal attribute levels}{
  Terms taking a categorical nodal covariate also take \code{levels}
  argument. This can be used to control the set and the ordering of
  attribute levels.
}


\section{Terms to represent network statistics included in the
  \code{\link[=ergm-package]{ergm}} package}{
  A cross-referenced html version of the term documentation is is available via \code{vignette('ergm-term-crossRef')} and terms can also be searched via \code{\link{search.ergmTerms}}.
  \describe{
    %%%%%%%%%%%%%%%%%%%%%%555
    % This section below is specially formatted so that it can be parsed to build additional 
    % documentation about ergm terms. (parsing functions are in build_term_index.R)
    % The list of term tags must appear on the same line as the term name to be parsed correctly
    %%beginTerms%%

  \item{\code{L(formula, Ls=~.)} (binary) (operator) (layer-aware)}{\emph{Evaluation on layers}:
    Evaluates the terms in \code{formula} on an actual or logical layers
    specified in formula \code{Ls} and sums the results. \code{Ls} is a
    one-sided formula with terms separated by \code{+} operators, and
    corresponding statistics for each layer specification is summed
    up. A dot (\code{.}) stands for summing all layers.
  }

  \item{\code{N(formula, lm=~1, subset=TRUE, weights=1, contrasts=NULL, offset=0, label=NULL)} (binary) (operator), \code{N(formula, lm=~1, subset=TRUE, weights=1, contrasts=NULL, offset=0, label=NULL)} (valued) (operator)}{\emph{Evaluation on multiple networks}:

    Evaluates the terms in \code{formula} on each of the networks joined
    using \code{\link{Networks}} function, and returns either a weighted
    sum or an \code{\link{lm}}-style linear model for the ERGM
    coefficients. Its syntax follows that of \code{\link{lm}} closely,
    with sensible defaults. The default formula (\code{~1}) sums the
    specified network statistics. If \code{lm} refers to any
    network attributes for which some networks have missing values, the
    term will stop with an error. This can be avoided by pre-filtering
    with \code{subset}, which controls which networks are affected by
    the term.

    Care should be taken to avoid multicollinearity when using this operator. When \code{lm} is given a model with intercept and a categorical predictor (including a \code{\link{logical}} one), it will use the first level (or \code{FALSE}) as the baseline, but if the model is without intercept, it will use all levels of the first categorical predictor. This is typically what is wanted in a linear regression, but for the \code{N} operator, this can be problematic if the \dQuote{intercept} effect is added by a different term. A workaround is to convert the categorical predictor to dummy variables before putting it into the \code{lm} formula.

    \code{offset} (which can be specified in \code{lm} as well) can be a constant, a vector of length equal to the number of networks, or a matrix whose number of rows is the number of networks and whose number of columns is the number of free parameters of the ERGM.

    Note that although \code{weights} is an argument in the term, it is
    reserved for future use; attempting to change it will cause an
    error: at this time, there is no way to assign sampling weights to
    networks.

    An optional parameter \code{label} will add a label to model
    parameters to help identify the term (which may have similar
    predictors but, say, a different network subset) in the output.
  }


    \item{\code{b1degreeL(d, by=NULL, levels=NULL)}  (binary)  (bipartite)  (undirected) (categorical nodal attribute) (frequently-used)
    }{\emph{Degree for the first mode in a
    bipartite (aka two-mode) network:} The \code{d} argument is a vector of
    distinct integers. This term adds one network statistic to the model for
    each element in \code{d}; the \eqn{i}th such statistic equals the number of
    nodes of degree \code{d[i]} in the first mode of a bipartite network, i.e.
    with exactly \code{d[i]} edges. The first mode of a bipartite network object
    is sometimes known as the "actor" mode. The optional argument \code{by} is
    a character string giving the name of an attribute in the network's vertex
    attribute list. If this is specified 
    then each node's degree is tabulated only with other nodes having the same
    value of the \code{by} attribute.
    This term can only be used with undirected bipartite networks.
    
    If a list of layer specifications is given, degree of a node
    \code{i} is considered to be the number of edges in all layers,
    combined.
}

    \item{\code{b2degrangeL(from, to=+Inf, by=NULL, homophily=FALSE, levels=NULL)}  (binary)  (bipartite)  (undirected)
    }{\emph{Degree range for the second mode in a bipartite
    (a.k.a. two-mode) network:} 
    The \code{from} and \code{to} arguments are vectors of distinct
    integers (or \code{+Inf}, for \code{to} (its default)). If one of the vectors has
    length 1, it is recycled to the length of the other. Otherwise, they
    must have the same length. This term adds one
    network statistic to the model for each element of \code{from} (or \code{to}); the \eqn{i}th
    such statistic equals the number of nodes of the second mode
    ("events") in the network of degree
    greater than or equal to 
    \code{from[i]} but strictly less than \code{to[i]}, i.e. with
    edge count
    in semiopen interval \code{[from,to)}. The optional argument
    \code{by} is a character string giving the name of an attribute in the
    network's vertex attribute list. 
    If this is specified and \code{homophily} is \code{TRUE},
    then degrees are calculated using the subnetwork consisting of only
    edges whose endpoints have the same value of the \code{by} attribute.
    If \code{by} is specified and
    \code{homophily} is \code{FALSE} (the default), then separate degree range
    statistics are calculated for nodes having each separate
    value of the attribute.
    
    This term can only be used with bipartite networks; for directed networks
    see \code{idegrange} and \code{odegrange}. For undirected networks,
    see \code{degrange}, and see \code{b1degrange}
    for degrees of the first mode ("actors").}
  
    \item{\code{b2degree(d, by=NULL)}  (binary)  (bipartite)  (undirected)  (categorical nodal attribute) (frequently-used)
    }{\emph{Degree for the second mode in a
    bipartite (aka two-mode) network:} The \code{d} argument is a vector of
    distinct integers. This term adds one network statistic to the model for
    each element in \code{d}; the \eqn{i}th such statistic equals the number of
    nodes of degree \code{d[i]} in the second mode of a bipartite network, i.e.
    with exactly \code{d[i]} edges. The second mode of a bipartite network
    object is sometimes known as the "event" mode. The optional term
    \code{by} is a character string giving the name of an attribute in the
    network's vertex attribute list. If this is specified 
    then each node's degree is tabulated only with other nodes having the same
    value of the \code{by} attribute.
    This term can only be used with undirected bipartite networks.

    If a list of layer specifications is given, degree of a node
    \code{i} is considered to be the number of edges in all layers,
    combined.
}
    
    \item{\code{ddspL(d, type="OTP", Ls.path=NULL, L.in_order=FALSE)} (binary) (directed) (undirected) (layer-aware)}{  \emph{Directed dyadwise shared partners:} This term adds one network statistic to the model for each element in \code{d} where the \eqn{i}th such statistic equals the number of \emph{dyads} in the network with exactly \code{d[i]} shared partners. This term can only be used with directed networks.  Multiple shared partner definitions are possible; the \code{type} argument may be used to select the type of shared partner to be counted (see below for type codes).  By default, outgoing two-paths are employed.
    
    While there is only one shared partner configuration in the undirected case, nine distinct configurations are possible for directed graphs.  Currently, edgewise shared partner terms may be defined with respect to five of these configurations; they are defined here as follows (using terminology from Butts (2008) and the \code{relevent} package):

  \describe{
    \item{Outgoing Two-path (OTP)}{vertex \eqn{k} is an OTP shared partner of ordered pair \eqn{(i,j)} iff \eqn{i \to k \to j}{i->k->j}. Also known as "transitive shared partner".}
    \item{Incoming Two-path (ITP)}{vertex \eqn{k} is an ITP shared partner of ordered pair \eqn{(i,j)} iff \eqn{j \to k \to i}{j->k->i}. Also known as "cyclical shared partner"}
    \item{Outgoing Shared Partner (OSP)}{vertex \eqn{k} is an OSP shared partner of ordered pair \eqn{(i,j)} iff \eqn{i \to k, j \to k}{i->k, j->k}.}
    \item{Incoming Shared Partner (ISP)}{vertex \eqn{k} is an ISP shared partner of ordered pair \eqn{(i,j)} iff \eqn{k \to i, k \to j}{k->i, k->j}.}
  }

  Note that Robins et al. (2009) define closely related statistics to several of the above, using slightly different terminology.

  For a multilayer network, a vector of one or two formulas
   \code{Ls.path} provides the layer specifications for the ties of the
   2-path or the shared partnership. (If only one formula is given the
   layers are assumed to be the same.) If \code{L.in_order==TRUE}, the
   first tie of the two-path must be the first element of \code{Ls.path}
   and the second must be the second; otherwise, any ordering counts, provided there is exactly one of each.

   This and related terms take an additional term option (via the
   \code{control$term.option} list), \code{cache.sp} controlling whether the
   implementation will cache the number of shared partners for each dyad
   in the network. This has a memory cost but may significuantly speed
   up calculation, particularly for denser networks.
 }

  \item{\code{degreeL(d, by=NULL, homophily=FALSE, levels=NULL, Ls=NULL)}  (binary)  (undirected)  (categorical nodal attribute) (frequently-used)}{\emph{Degree:} 
    The \code{d} argument is a vector of distinct integers. This term adds one
    network statistic to the model for each element in \code{d}; the \eqn{i}th
    such statistic equals the number of nodes in the network of degree
    \code{d[i]}, i.e. with exactly \code{d[i]} edges. The optional argument
    \code{by} is a character string giving the name of an attribute in the
    network's vertex attribute list. 
    If this is specified and \code{homophily} is \code{TRUE},
    then degrees are calculated using the subnetwork consisting of only
    edges whose endpoints have the same value of the \code{by} attribute.
    If \code{by} is specified and
    \code{homophily} is \code{FALSE} (the default), then separate degree 
    statistics are calculated for nodes having each separate
    value of the attribute.
    This term can only be used with undirected networks; for directed networks
    see \code{idegree} and \code{odegree}.

    If a list of layer specifications is given, degree of a node
    \code{i} is considered to be the number of edges in all layers,
    combined.
  }

  \item{\code{despL(d, type="OTP", L.base=NULL, Ls.path=NULL, L.in_order=FALSE)} (binary) (directed) (undirected) (layer-aware)}{  \emph{Directed edgewise shared partners:} This term adds one network statistic to the model for each element in \code{d} where the \eqn{i}th such statistic equals the number of \emph{edges} in the network with exactly \code{d[i]} shared partners. This term can only be used with directed networks.  Multiple shared partner definitions are possible; the \code{type} argument may be used to select the type of shared partner to be counted (see \code{ddsp} for type codes).  By default, outgoing two-paths are employed.

    For a multilayer network, a vector of one or two formulas
   \code{Ls.path} provides the layer specifications for the ties of the
   2-path or the shared partnership. (If only one formula is given the
   layers are assumed to be the same.) If \code{L.in_order==TRUE}, the
   first tie of the two-path must be the first element of \code{Ls.path}
   and the second must be the second; otherwise, any ordering counts,
   provided there is exactly one of each. (For types \code{"OSP"} and \code{"ISP"}, the \dQuote{first} tie is considered to be the one incident on the tail of the base tie.)

   \code{L.base} gives the layer specification for the base.

  }
    
  \item{\code{dgwdspL(decay, fixed=FALSE, cutoff=30, type="OTP", Ls.path=NULL, L.in_order=FALSE)} (binary) (directed) (undirected) (layer-aware)}{  \emph{Geometrically weighted dyadwise shared partner distribution:} This term adds one network statistic to the model equal to the geometrically weighted dyadwise shared partner distribution with decay parameter \code{decay} parameter, which should be non-negative. (this parameter was called \code{alpha} prior to \code{ergm 3.7}). The value supplied for this parameter may be fixed (if \code{fixed=TRUE}), or it may be used instead as the starting value for the estimation of \code{decay} in a curved exponential family model (when \code{fixed=FALSE}, the default) (see Hunter and Handcock, 2006).  Note that the GWDSP statistic is equal to the sum of GWNSP plus GWESP. For a directed network, multiple shared partner definitions are possible; the \code{type} argument may be used to select the type of shared partner to employ (see \code{ddsp} for definitions).  By default, outgoing two-paths are employed.  The optional argument \code{cutoff} sets the number of underlying DSP terms to use in computing the statistics to reduce the computational burden.  Its default value can also be controlled by the \code{gw.cutoff} term option control parameter. (See \code{\link{control.ergm}}.)

    See term \code{ddsp} for layer specification information.
   }
   
   \item{\code{dgwespL(decay, fixed=FALSE, cutoff=30, type="OTP", L.base=NULL, Ls.path=NULL, L.in_order=FALSE)} (binary) (directed) (undirected) (layer-aware)}{ \emph{Geometrically weighted edgewise shared partner distribution:} This term adds a statistic equal to the geometrically weighted \emph{edgewise} (not dyadwise) shared partner distribution with decay parameter \code{decay} parameter, which should be non-negative. (this parameter was called \code{alpha} prior to \code{ergm 3.7}). The value supplied for this parameter may be fixed (if \code{fixed=TRUE}), or it may be used instead as the starting value for the estimation of \code{decay} in a curved exponential family model (when \code{fixed=FALSE}, the default) (see Hunter and Handcock, 2006). For a directed network, multiple shared partner definitions are possible; the \code{type} argument may be used to select the type of shared partner to employ (see \code{ddsp} for definitions).  By default, outgoing two-paths are employed.  The optional argument \code{cutoff} sets the number of underlying ESP terms to use in computing the statistics to reduce the computational burden.  Its default value can also be controlled by the \code{gw.cutoff} term option control parameter. (See \code{\link{control.ergm}}.)

     See term \code{desp} for layer specification information.
   }
  
   \item{\code{dgwnspL(decay, fixed=FALSE, cutoff=30, type="OTP", L.base=NULL, Ls.path=NULL, L.in_order=FALSE)} (binary) (directed) (undirected) (layer-aware)}{ \emph{Geometrically weighted non-edgewise shared partner distribution:} This term is just like gwesp and gwdsp except it adds a statistic equal to the geometrically weighted nonedgewise (that is, over dyads that do not have an edge) shared partner distribution with decay parameter \code{decay} parameter, which should be non-negative. (this parameter was called \code{alpha} prior to \code{ergm 3.7}). The value supplied for this parameter may be fixed (if \code{fixed=TRUE}), or it may be used instead as the starting value for the estimation of \code{decay} in a curved exponential family model (when \code{fixed=FALSE}, the default) (see Hunter and Handcock, 2006). For a directed network, multiple shared partner definitions are possible; the \code{type} argument may be used to select the type of shared partner to employ (see \code{ddsp} for definitions).  By default, outgoing two-paths are employed.  The optional argument \code{cutoff} sets the number of underlying NSP terms to use in computing the statistics to reduce the computational burden.  Its default value can also be controlled by the \code{gw.cutoff} term option control parameter. (See \code{\link{control.ergm}}.)

     See term \code{desp} for layer specification information.
   }
    
   \item{\code{dnspL(d, type="OTP", L.base=NULL, Ls.path=NULL, L.in_order=FALSE)} (binary) (directed) (undirected) (layer-aware)}{  \emph{Directed non-edgewise shared partners:} This term adds one network statistic to the model for each element in \code{d} where the \eqn{i}th such statistic equals the number of \emph{non-edges} in the network with exactly \code{d[i]} shared partners. This term can only be used with directed networks.  Multiple shared partner definitions are possible; the \code{type} argument may be used to select the type of shared partner to be counted (see \code{ddsp} for type codes).  By default, outgoing two-paths are employed.

     See term \code{desp} for layer specification information.     
   }

    \item{\code{gwb1degreeL(decay, fixed=FALSE, cutoff=30, levels=NULL, Ls=NULL)}  (binary)  (bipartite)  (undirected) (curved)}{\emph{Geometrically weighted
    degree distribution for the first mode in a bipartite (aka two-mode)
    network:} 
    This term adds one network statistic to the model equal to the weighted
    degree distribution with decay controlled by the \code{decay} parameter, which should be non-negative,
    for nodes in the
    first mode of a bipartite network. The first mode of a bipartite network
    object is sometimes known as the "actor" mode. 
    The \code{decay} parameter is the same as theta_s in
    equation (14) in Hunter (2007).  The value supplied for
    this parameter may be fixed (if \code{fixed=TRUE}),
    or it may be used as merely the starting value for the estimation
    in a curved exponential family model (the default).
    The optional argument \code{cutoff}
    is only relevant if \code{fixed=FALSE}. In that case it only uses this 
    number of terms in computing the statistics to reduce the computational
    burden.  Its default value can also be controlled by the \code{gw.cutoff} term option control parameter. (See \code{\link{control.ergm}}.) If \code{attrname} is specified then separate degree
    statistics are calculated for nodes having each separate
    value of the attribute.
    This term can only be used with undirected bipartite
    networks.

    If a list of layer specifications is given, degree of a node
    \code{i} is considered to be the number of edges in all layers,
    combined.
  }
    
    \item{\code{gwb2degreeL(decay, fixed=FALSE, attrname=NULL, cutoff=30, levels=NULL, Ls=NULL)}  (binary)  (bipartite)  (undirected) (curved)}{\emph{Geometrically weighted
    degree distribution for the second mode in a bipartite (aka two-mode)
    network:} 
    This term adds one network statistic to the model equal to the weighted
    degree distribution with decay controlled by the which should be non-negative,
    for nodes in the
    second mode of a bipartite network. The second mode of a bipartite network
    object is sometimes known as the "event" mode. 
    The \code{decay} parameter is the same as theta_s in
    equation (14) in Hunter (2007).  The value supplied for
    this parameter may be fixed (if \code{fixed=TRUE}),
    or it may be used as merely the starting value for the estimation
    in a curved exponential family model (the default).
    The optional argument \code{cutoff}
    is only relevant if \code{fixed=FALSE}. In that case it only uses this 
    number of terms in computing the statistics to reduce the computational
    burden.  Its default value can also be controlled by the \code{gw.cutoff} term option control parameter. (See \code{\link{control.ergm}}.) If \code{attrname} is specified then separate degree
    statistics are calculated for nodes having each separate
    value of the attribute.
    This term can only be used with undirected bipartite
    networks.

    If a list of layer specifications is given, degree of a node
    \code{i} is considered to be the number of edges in all layers,
    combined.
  }
    
    \item{\code{gwdegreeL(decay, fixed=FALSE, attrname=NULL, cutoff=30, levels=NULL)}  (binary)  (undirected)  (curved) (frequently-used)
    }{\emph{Geometrically weighted
    degree distribution:} 
    This term adds one network statistic to the model equal to the weighted
    degree distribution with decay controlled by the \code{decay} parameter.
    The \code{decay} parameter is the same as theta_s in
    equation (14) in Hunter (2007).  The value supplied for
    this parameter may be fixed (if \code{fixed=TRUE}),
    or it may be used instead as the starting value for the estimation of \code{decay}
    in a curved exponential family model (when \code{fixed=FALSE}, the default) (see Hunter and Handcock, 2006).
    The optional argument \code{cutoff}
    is only relevant if \code{fixed=FALSE}. In that case it only uses this 
    number of terms in computing the statistics to reduce the computational
    burden.  Its default value can also be controlled by the \code{gw.cutoff} term option control parameter. (See \code{\link{control.ergm}}.)    This term can only be used with undirected networks.

    If a list of layer specifications is given, degree of a node
    \code{i} is considered to be the number of edges in all layers,
    combined.

    If \code{attrname} is specified then separate degree
    statistics are calculated for nodes having each separate
    value of the attribute.  This term can only be used with undirected networks.}

    \item{\code{gwidegreeL(decay, fixed=FALSE, attrname=NULL, cutoff=30, levels=NULL, Ls=NULL)}  (binary)  (directed) (curved)}{\emph{Geometrically weighted in-degree distribution:} This term adds one network statistic to the model
    equal to the weighted in-degree distribution with decay parameter
    \code{decay} parameter, which should be non-negative. (this parameter was called \code{alpha} prior to \code{ergm 3.7}). The value supplied for
    this parameter may be fixed (if \code{fixed=TRUE}),
    or it may be used instead as the starting value for the estimation of \code{decay}
    in a curved exponential family model (when \code{fixed=FALSE}, the default) (see Hunter and Handcock, 2006). This
    term can only be used with directed networks.
    The optional argument \code{cutoff}
    is only relevant if \code{fixed=FALSE}. In that case it only uses this 
    number of terms in computing the statistics to reduce the computational
    burden.  Its default value can also be controlled by the \code{gw.cutoff} term option control parameter. (See \code{\link{control.ergm}}.)

    If a list of layer specifications is given, degree of a node
    \code{i} is considered to be the number of edges in all layers,
    combined.

    If \code{attrname} is specified then separate degree
    statistics are calculated for nodes having each separate
    value of the attribute.
  }


  \item{\code{gwdegreeLL(decay, fixed=FALSE, attrname=NULL, cutoff=30, levels=NULL, Ls, dir)}  (binary)  (directed) (curved)}{\emph{Geometrically weighted
    layer-combination-degree distribution:} This term adds one network statistic to the model
    equal to the weighted layer-cumulative degree distribution with decay parameter
    \code{decay} parameter, which should be non-negative. (this parameter was called \code{alpha} prior to \code{ergm 3.7}). The value supplied for
    this parameter may be fixed (if \code{fixed=TRUE}),
    or it may be used instead as the starting value for the estimation of \code{decay}
    in a curved exponential family model (when \code{fixed=FALSE}, the default) (see Hunter and Handcock, 2006). This
    term can only be used with directed networks.
    The optional argument \code{cutoff}
    is only relevant if \code{fixed=FALSE}. In that case it only uses this 
    number of terms in computing the statistics to reduce the computational
    burden.  Its default value can also be controlled by the \code{gw.cutoff} term option control parameter. (See \code{\link{control.ergm}}.)

    In this term, the layer-cumulative degree frequency equals
    the number of nodes in the network whose total degree across all
    layers listed in \code{Ls} combined (with direction being
    corresponding element of vector \code{dir}, \code{"in"} or \code{"out"}), is \code{d[i]}.

    If \code{attrname} is specified then separate degree
    statistics are calculated for nodes having each separate
    value of the attribute.}

    \item{\code{gwodegreeL(decay, fixed=FALSE, attrname=NULL, cutoff=30, levels=NULL, Ls=NULL)}  (binary)  (directed) (curved)}{\emph{Geometrically weighted out-degree distribution:} This term adds one network statistic to the model
    equal to the weighted out-degree distribution with decay parameter
    \code{decay} parameter, which should be non-negative. (this parameter was called \code{alpha} prior to \code{ergm 3.7}). The value supplied for
    this parameter may be fixed (if \code{fixed=TRUE}),
    or it may be used instead as the starting value for the estimation of \code{decay}
    in a curved exponential family model (when \code{fixed=FALSE}, the default) (see Hunter and Handcock, 2006). This
    term can only be used with directed networks.
    The optional argument \code{cutoff}
    is only relevant if \code{fixed=FALSE}. In that case it only uses this 
    number of terms in computing the statistics to reduce the computational
    burden.  Its default value can also be controlled by the \code{gw.cutoff} term option control parameter. (See \code{\link{control.ergm}}.)

    If a list of layer specifications is given, degree of a node
    \code{i} is considered to be the number of edges in all layers,
    combined.

    If \code{attrname} is specified then separate degree
    statistics are calculated for nodes having each separate
    value of the attribute.}

  \item{\code{idegreeL(d, by=NULL, homophily=FALSE, levels=NULL)}  (binary)  (directed)  (categorical nodal attribute) (frequently-used)}{\emph{In-degree:} The \code{d} argument
    is a vector of distinct integers. This term adds one network statistic to
    the model for each element in \code{d}; the \eqn{i}th such statistic equals
    the number of nodes in the network of in-degree \code{d[i]}, i.e. the number
    of nodes with exactly \code{d[i]} in-edges. The optional term
    \code{by} is a character string giving the name of an attribute in the
    network's vertex attribute list. 
    If this is specified and \code{homophily} is \code{TRUE},
    then degrees are calculated using the subnetwork consisting of only
    edges whose endpoints have the same value of the \code{by} attribute.
    If \code{by} is specified and
    \code{homophily} is \code{FALSE} (the default), then separate degree 
    statistics are calculated for nodes having each separate
    value of the attribute.
    This term can only be used with directed networks; for undirected networks
    see \code{degree}.

    If a list of layer specifications is given, degree of a node
    \code{i} is considered to be the number of edges in all layers,
    combined.
  }

  \item{\code{degreeL(d, by=NULL, Ls, dir)}  (binary)  (directed)  (categorical nodal attribute) (frequently-used)}{\emph{Layer-combinaton-degree:} The \code{d} argument
    is a vector of distinct integers. This term adds one network statistic to
    the model for each element in \code{d}; the \eqn{i}th such statistic equals
    the number of nodes in the network whose total degree across all
    layers listed in \code{Ls} combined (with direction being
    corresponding element of vector \code{dir}, \code{"in"} or \code{"out"}), is \code{d[i]}.
    If \code{by} is specified, then separate degree 
    statistics are calculated for nodes having each separate
    value of the attribute.
    This term can only be used with directed networks; for undirected networks
    see \code{degree}.
  }

    \item{\code{CMBL(Ls=~.)} (binary) (directed) (undirected) (layer-aware)}{\emph{Conway--Maxwell-Binomial dependence among layers}:

      Models marginal dependence layers within each dyad by imposing
      a Conway--Maxwell-Binomial (CMB) distribution on the number of
      layers in each dyad that have a tie.

      \code{Ls} is a list of one-sided formulas indicating the layers
      whose dependence is being modeled. A dot (\code{.}) stands for all
      layers in the network.

      The term adds one statistic to the model, equalling the sum over
      all the dyads in the network of \eqn{\log\{E!(R-E)!/R!\}}, where
      \eqn{E} is the number of layers in \code{Ls} with an edge in that
      dyad and \eqn{R} being the total number of layers in \code{Ls}.

      A positive coefficient induces positive dependence and a negative
      one induces negative dependence.
  }

    \item{\code{mutualL(same=NULL, diff=FALSE, by=NULL, keep=NULL, Ls=NULL)} (binary)  (directed)  (frequently-used) (layer-aware), \code{mutual(form="min",threshold=0)}  (valued)  (directed)
     }{\emph{Mutuality:} 
    In binary ERGMs, equal to the number of
    pairs of actors \eqn{i} and \eqn{j} for which \eqn{(i{\rightarrow}j)}{(i,j)}
    and \eqn{(j{\rightarrow}i)}{(j,i)} both exist. For valued ERGMs, equal to \eqn{\sum_{i<j} m(y_{i,j},y_{j,i})},
      where \eqn{m} is determined by \code{form} argument: \code{"min"}
      for \eqn{\min(y_{i,j},y_{j,i})}, \code{"nabsdiff"} for
      \eqn{-|y_{i,j},y_{j,i}|}, \code{"product"} for
      \eqn{y_{i,j}y_{j,i}}, and \code{"geometric"} for
      \eqn{\sqrt{y_{i,j}}\sqrt{y_{j,i}}}. See Krivitsky (2012) for a
      discussion of these statistics. \code{form="threshold"} simply
      computes the binary \code{mutuality} after
      thresholding at \code{threshold}.

    This term can only be used
    with directed networks. The binary version also has the following capabilities: if the optional \code{same} argument is passed
    the name of a vertex attribute,
    only mutual pairs that match on the attribute are counted;
    separate counts for each unique matching value can be obtained by using
    \code{diff=TRUE} with \code{same}; and if \code{by} is passed the name of a vertex attribute,
    then each node is counted separately for each mutual pair in which it 
    occurs and the counts are tabulated by unique values of the attribute.
    This means that the sum of the mutual statistics when \code{by} is used
    will equal twice the standard mutual statistic.  Only one of \code{same}
    or \code{by} may be used, and only the former is affected by \code{diff};
    if both \code{same} and \code{by} are passed, \code{by} is ignored.
    Finally, if \code{keep} is passed a numerical vector, this vector of integers
    tells which statistics should be kept whenever the \code{mutual} term would
    ordinarily result in multiple statistics.

    If layer specification is given in argument \code{Ls} (as either a
    list of two layer logic formulas or a single layer logic formula),
    the statistic will count the number of dyads where a tie in
    \code{Ls[[1]]} reciprocates a tie in \code{Ls[[2]]} and vice
    versa. (Note that dyad that has mutual ties in \code{Ls[[1]]} and in
    \code{Ls[[2]]} will add 2 to this statistic.) If a formula is given,
    it is replicated.
  }

    \item{\code{odegreeL(d, by=NULL, homophily=FALSE, levels=NULL)}  (binary)  (directed)  (categorical nodal attribute) (frequently-used)}{\emph{Out-degree:} The \code{d} argument
    is a vector of distinct integers. This term adds one network statistic to
    the model for each element in \code{d}; the \eqn{i}th such statistic equals
    the number of nodes in the network of out-degree \code{d[i]}, i.e. the
    number of nodes with exactly \code{d[i]} out-edges. The optional argument
    \code{by} is a character string giving the name of an attribute in the
    network's vertex attribute list. 
    If this is specified and \code{homophily} is \code{TRUE},
    then degrees are calculated using the subnetwork consisting of only
    edges whose endpoints have the same value of the \code{by} attribute.
    If \code{by} is specified and
    \code{homophily} is \code{FALSE} (the default), then separate degree 
    statistics are calculated for nodes having each separate
    value of the attribute.
    This term can only be used with directed networks; for undirected networks
    see \code{degree}.

    If a list of layer specifications is given, degree of a node
    \code{i} is considered to be the number of edges in all layers,
    combined.
}

  \item{\code{twostarL(Ls, type, distinct=TRUE}  (binary)  (directed) (layer-aware)}{\emph{Multilayer two-star:} 
    This term adds one statistic to the model, equal to the number of
    cross-layer two-stars or two-paths in
    the network. Argument \code{Ls} specifies the two layers of
    interest, and argument \code{type} determines what is counted:
    \describe{
      \item{\code{"out"}}{Number of configurations
	\eqn{(i{\rightarrow}j), (i{\rightarrow}k)}{(i,j), (i,k)}, where
	\eqn{(i{\rightarrow}j)}{(i,j)} is in logical layer \code{Ls[[1]]}
	and \eqn{(i{\rightarrow}k)}{(i,k)} is in logical layer \code{Ls[[2]]}.}
      \item{\code{"in"}}{Number of configurations
	\eqn{(j{\rightarrow}i), (k{\rightarrow}i)}{(j,i), (k,i)}, where
	\eqn{(j{\rightarrow}i)}{(j,i)} is in logical layer \code{Ls[[1]]}
	and \eqn{(k{\rightarrow}i)}{(k,i)} is in logical layer \code{Ls[[2]]}.}
      \item{\code{"path"}}{Number of configurations
	\eqn{(j{\rightarrow}i), (i{\rightarrow}k)}{(j,i), (i,k)}, where
	\eqn{(j{\rightarrow}i)}{(j,i)} is in logical layer \code{Ls[[1]]}
	and \eqn{(i{\rightarrow}k)}{(i,k)} is in logical layer
	\code{Ls[[2]]}.}
    }

    If \code{distinct==TRUE}, \eqn{j} and \eqn{k} above are required to
    be distinct. That is, the constituent edges may not be coincident or
    reciprocal.
  }
 }
}
\references{
\itemize{
  \item Bomiriya, R. P, Bansal, S., and Hunter, D. R. (2014).  Modeling
  Homophily in ERGMs for Bipartite Networks.  Submitted.
  
  \item Butts, CT.  (2008).  \dQuote{A Relational Event Framework for Social Action.} \emph{Sociological Methodology,} 38(1).
  
  \item Davis, J.A. and Leinhardt, S.  (1972).  The Structure of Positive Interpersonal
  Relations in Small Groups.  In J. Berger (Ed.), \emph{Sociological Theories in Progress,
  Volume 2}, 218--251.  Boston: Houghton Mifflin.
  
  \item Holland, P. W. and S. Leinhardt (1981). An exponential family
  of probability distributions for directed graphs.
  \emph{Journal of the American Statistical Association}, 76:
  33--50.
  
  \item Hunter, D. R. and M. S. Handcock (2006). Inference in curved
  exponential family models for networks. \emph{Journal of Computational
  and Graphical Statistics}, 15: 565--583.
  
  \item Hunter, D. R. (2007). Curved exponential family models for social
  networks. \emph{Social Networks}, 29: 216--230.
  
  \item Krackhardt, D. and Handcock, M. S. (2007).  Heider versus Simmel: Emergent Features
  in Dynamic Structures. \emph{Lecture Notes in Computer Science}, 4503, 14--27.

  \item Krivitsky P. N. (2012). Exponential-Family Random Graph Models for Valued
    Networks. \emph{Electronic Journal of Statistics}, 2012, 6,
    1100-1128. \doi{10.1214/12-EJS696}
    
  \item Robins, G; Pattison, P; and Wang, P.  (2009).  \dQuote{Closure, Connectivity, and Degree Distributions: Exponential Random Graph (p*) Models for Directed Social Networks.}  \emph{Social Networks,} 31:105-117.

  \item Snijders T. A. B., G. G. van de Bunt, and
    C. E. G. Steglich. Introduction to Stochastic Actor-Based Models for
    Network Dynamics. \emph{Social Networks}, 2010, 32(1),
    44-60. \doi{10.1016/j.socnet.2009.02.004}

  \item Morris M, Handcock MS, and Hunter DR. Specification of Exponential-Family Random Graph Models: Terms and Computational Aspects. \emph{Journal of Statistical Software}, 2008, 24(4), 1-24. \url{http://www.jstatsoft.org/v24/i04}
    
    
  \item Snijders, T. A. B., P. E. Pattison, G. L. Robins, and
  M. S. Handcock (2006). New specifications for exponential random graph
  models, \emph{Sociological Methodology}, 36(1): 99-153. 
  }
}
\seealso{\code{\link[=ergm-package]{ergm}} package, \code{\link{search.ergmTerms}}, \code{\link{ergm}}, \code{\link{network}}, \code{\link{\%v\%}}, \code{\link{\%n\%}}}
\examples{
\dontrun{
ergm(flomarriage ~ kstar(1:2) + absdiff("wealth") + triangle)

ergm(molecule ~ edges + kstar(2:3) + triangle
                      + nodematch("atomic type",diff=TRUE)
                      + triangle + absdiff("atomic type"))
}
% TODO: Write a valued example.
}
\keyword{models}
