# `ergm.multi`: Fit, Simulate and Diagnose Exponential-Family Models for Multiple or Multilayer Networks

[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/ergm.multi?color=2ED968)](https://cranlogs.r-pkg.org/)
[![cran version](https://www.r-pkg.org/badges/version/ergm.multi)](https://cran.r-project.org/package=ergm.multi)
[![R build status](https://github.com/statnet/ergm.multi/workflows/R-CMD-check/badge.svg)](https://github.com/statnet/ergm.multi/actions)
[![R-universe](https://statnet.r-universe.dev/ergm.multi/badges/version)](https://statnet.r-universe.dev/ergm.multi)

A set of extensions for the 'ergm' package to fit multilayer/multiplex/multirelational networks and samples of multiple networks. 'ergm.multi' is a part of the Statnet suite of packages for network analysis. See Krivitsky, Koehly, and Marcum (2020) <doi:10.1007/s11336-020-09720-7> and Krivitsky, Coletti, and Hens (2023) <doi:10.1080/01621459.2023.2242627>.

## Public and Private repositories

To facilitate open development of the package while giving the core developers an opportunity to publish on their developments before opening them up for general use, this project comprises two repositories:
* A public repository `statnet/ergm.multi`
* A private repository `statnet/ergm.multi-private`

The intention is that all developments in `statnet/ergm.multi-private` will eventually make their way into `statnet/ergm.multi` and onto CRAN.

Developers and Contributing Users to the Statnet Project should read https://statnet.org/private/ for information about the relationship between the public and the private repository and the workflows involved.

## Latest Windows and MacOS binaries

[R-Universe](https://r-universe.dev) builds a set of binaries after every commit to the main branch of the repository. We strongly encourage testing against them before filing a bug report, as they may contain fixes that have not yet been sent to CRAN. To obtain the binaries from r-universe, navigate to the package page at https://statnet.r-universe.dev/ergm.multi .
You may also want to install the corresponding latest binaries for packages on which `ergm.multi` depends, in particular [`statnet.common`](https://github.com/statnet/statnet.common), [`network`](https://github.com/statnet/network), [`networkLite`](https://github.com/EpiModel/networkLite), and [`ergm`](https://github.com/statnet/ergm).
