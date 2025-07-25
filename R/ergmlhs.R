#  File R/ergmlhs.R in package ergm.multi, part of the Statnet suite of
#  packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################

#' Combine the [`%ergmlhs%`] settings of a list of networks
#'
#' This is a helper function to go through the [`network`] objects in the list
#' and accumulate their [`%ergmlhs%`] settings, printing a message if
#' the settings clash. Later networks' settings overwrite the earlier.
#'
#' @param nwl a [`list`] of [`network`] objects whose settings are to
#'   be combined.
#' @param ignore.settings a [`character`] vector of setting names to
#'   be ignored.
#'
#' @return A [`list`] of settings, suitable for being assigned as the
#'   `ergm` network attribute.
#'
#' @keywords internal
#' @export
combine_ergmlhs <- function(nwl, ignore.settings=c()){
  ergml <- nwl %>% map(`%n%`, "ergm") %>% map(NVL, list())
  l <- ergml[[1]]

  for(i in seq_along(ergml)[-1]){
    l2 <- ergml[[i]]
    all_names <- setdiff(union(names(l),names(l2)), ignore.settings)
    for(name in all_names){
      if(!identical(l[[name]], l2[[name]])){
        message("Setting ", sQuote(name), " differs between network ", i, " and some prior network and will overwrite it. Is probably benign (e.g., environment of a formula that does not reference its environment), or the networks may have inconsistent ERGM settings.")
        l[[name]] <- l2[[name]]
      }
    }
  }
  l
}
