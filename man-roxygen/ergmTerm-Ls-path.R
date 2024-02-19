#  File man-roxygen/ergmTerm-Ls-path.R in package ergm.multi, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2024 Statnet Commons
################################################################################
#' @param Ls.path,L.in_order a vector of one or two formulas `Ls.path`
#'   provides the Layer Logic (c.f. Layer Logic section in the
#'   [Layer()] documentation) specifications for the ties of the
#'   2-path or the shared partnership. (If only one formula is given
#'   the layers are assumed to be the same.) If `L.in_order==TRUE` ,
#'   the first tie of the two-path must be the first element of
#'   `Ls.path` and the second must be the second; otherwise, any
#'   ordering counts, provided there is exactly one of each. (For
#'   types `"OSP"` and `"ISP"` , the first tie is considered to be the
#'   one incident on the tail of the base tie.)
