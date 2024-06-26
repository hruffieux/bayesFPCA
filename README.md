<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- First time: run usethis::use_readme_rmd() to create a pre-commit hook that 
prevents from committing if the README.Rmd has changed, but has not been 
re-knitted to generate an updated README.md -->

## bayesFPCA - Bayesian Functional Principal Component Analysis Suite <img src="man/figures/bayesFPCA_logo.png" align="right" height="150"/>

<!-- Run for the R CMD checks, run usethis::use_github_actions() to set up the pipeline, possibly modify the .yaml file and then: -->
<!-- [![R build status](https://github.com/hruffieux/bayesFPCA/workflows/R-CMD-check/badge.svg)](https://github.com/hruffieux/bayesFPCA/actions) -->
<!-- [![](https://travis-ci.org/hruffieux/bayesFPCA.svg?branch=master)](https://travis-ci.org/hruffieux/bayesFPCA) -->

[![License: GPL
v3](https://img.shields.io/badge/license-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![](https://img.shields.io/badge/devel%20version-0.1.0-blue.svg)](https://github.com/hruffieux/bayesFPCA)
[![](https://img.shields.io/github/languages/code-size/hruffieux/bayesFPCA.svg)](https://github.com/hruffieux/bayesFPCA)

## Overview

**bayesFPCA** is an R package providing tools for univariate and
multivariate functional principal component analysis (FPCA) in the
Bayesian setting. It provides efficient variational inference
implementations (mean-field and variational message passing).

## Installation

Then, to install the package in R, run the following command:

``` r
if(!require(remotes)) install.packages("remotes")
remotes::install_github("hruffieux/bayesFPCA")
```

## Main functions

The two main functions to perform univariate and multivariate FPCA
inference are `run_mfvb_fpca()` for the mean-field variational Bayes
implementation and `run_vmp_fpca()` for the variational message passing
implementation.

The package also provides a number of functions to generate and display
functional data as well as the resulting FPCA estimates.

## License and authors

This software uses the GPL v3 license, see [LICENSE](LICENSE). Authors
and copyright are provided in [DESCRIPTION](DESCRIPTION).

## Issues

To report an issue, please use the [bayesFPCA issue
tracker](https://github.com/hruffieux/bayesFPCA/issues) at github.com.
