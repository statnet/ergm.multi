# `ergm.multi`: Fit, Simulate and Diagnose Exponential-Family Models for Multiple or Multilayer Networks

[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/ergm.multi?color=2ED968)](https://cranlogs.r-pkg.org/)
[![cran version](https://www.r-pkg.org/badges/version/ergm.multi)](https://cran.r-project.org/package=ergm.multi)
[![R build status](https://github.com/statnet/ergm.multi/workflows/R-CMD-check/badge.svg)](https://github.com/statnet/ergm.multi/actions)

A set of extensions for the 'ergm' package to fit multilayer/multiplex/multirelational networks and samples of multiple networks. 'ergm' is a part of the Statnet suite of packages for network analysis.

## Public and Private repositories

To facilitate open development of the package while giving the core developers an opportunity to publish on their developments before opening them up for general use, this project comprises two repositories:
* A public repository `statnet/ergm.multi`
* A private repository `statnet/ergm.multi-private`

The intention is that all developments in `statnet/ergm.multi-private` will eventually make their way into `statnet/ergm.multi` and onto CRAN.

Developers and Contributing Users to the Statnet Project should read https://statnet.github.io/private/ for information about the relationship between the public and the private repository and the workflows involved.

## Latest Windows and MacOS binaries

A set of binaries is built after every commit to the repository. We strongly encourage testing against them before filing a bug report, as they may contain fixes that have not yet been sent to CRAN. They can be downloaded through the following links:

* [MacOS binary (a `.tgz` file in a `.zip` file)](https://nightly.link/statnet/ergm.multi/workflows/R-CMD-check.yaml/master/macOS-rrelease-binaries.zip)
* [Windows binary (a `.zip` file in a `.zip` file)](https://nightly.link/statnet/ergm.multi/workflows/R-CMD-check.yaml/master/Windows-rrelease-binaries.zip)

You will need to extract the MacOS `.tgz` or the Windows `.zip` file from the outer `.zip` file before installing. These binaries are usually built under the latest version of R and their operating system and may not work under other versions.

You may also want to install the corresponding latest binaries for packages on which `ergm.multi` depends, in particular [`statnet.common`](https://github.com/statnet/statnet.common), [`network`](https://github.com/statnet/network), and [`ergm`](https://github.com/statnet/ergm).
