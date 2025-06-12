## These will eventually be moved to 'statnet.common'.

dbl_along <- function(x) numeric(length(x))
int_along <- function(x) integer(length(x))
chr_along <- function(x) character(length(x))
lgl_along <- function(x) logical(length(x))

## This should be in 'network'.

b1.size <- function(x) if(is.bipartite(x)) x %n% "bipartite" else FALSE
b2.size <- function(x) if(is.bipartite(x)) network.size(x) - b1.size(x) else FALSE
