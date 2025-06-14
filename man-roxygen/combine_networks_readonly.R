#  File man-roxygen/combine_networks_readonly.R in package ergm.multi, part of
#  the Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################
#' @return Due to certain optimizations, the resulting [`network`] object's network and vertex attributes should be treated as read-only: do not modify them. If you need to change existing attributes or add new ones, do so on the input networks and call `<%= combiner %>(...)` again.
