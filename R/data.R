#' A sample of within-household contact networks in Flanders and Brussels
#'
#' This is a list of 318 [`network`] objects derived from contact
#' diary data collected by by Goeyvaerts et al. (2018). The study
#' recruited households in Flanders and Brussels-Capital region with
#' at least one child 12 or under. The networks are symmetrized.
#'
#' @section Nonstandard Network Attributes: \describe{
#'
#' \item{`included`}{(logical) whether the network was included in
#' Goeyvaerts's analysis. (Two were excluded.)}
#'
#' \item{`weekday`}{(logical) whether the contact diary on which the
#' network is based was collected on a weekday, as opposed to
#' weekend.}
#'
#' }
#'
#' @section Nonstandard Vertex Attributes: \describe{
#'
#' \item{`age`}{(numeric) the household member's age.}
#'
#' \item{`gender`}{(character) the household member's gender (`"F"`/`"M"`).}
#'
#' \item{`role`}{(character) the household member's inferred role (`"Father"`/`"Mother"`/`"Child"`/`"Grandmother"`).}
#'
#' }
#'
#'
#' @usage
#' data(Goeyvaerts)
#' @docType data
#' @section Licenses and Citation: When publishing results obtained using this
#' data set, the original authors (Goeyvaerts et al., 2018)
#' should be cited, along with this \R package.
#' @references
#'
#' Nele Goeyvaerts, Eva Santermans, Gail Potter, Andrea Torneri, Kim
#' V. Kerckhove, Lander Willem, Marc Aerts, Philippe Beutels & Niel
#' Hens (2018) Household Members Do Not Contact Each Other at Random:
#' Implications for Infectious Disease Modelling. *Proceedings of the
#' Royal Society B: Biological Sciences*, 285(1893):
#' 20182201. \doi{10.1098/rspb.2018.2201}
#'
#' @source The data were collected and by Goeyvaerts et al. and
#'   curated by Pietro Coletti.
#'
#' @seealso `vignette("Goeyvaerts_reproduction")` for a vignette reproducing the Goeyvaerts analysis and performing diagnostics
#' @keywords datasets
"Goeyvaerts"
