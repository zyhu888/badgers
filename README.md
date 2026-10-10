**badgers**: Bayesian Graph-Informed Regression with Measurement Error
in Predictors
================

[![License:
GPL-3](https://img.shields.io/badge/License-GPL--3-green.svg)](LICENSE.md)

## badgers

This R package fits BADGERS, a Bayesian regression model for
high-dimensional predictors that are observed with measurement error
through replicate measurements and that are linked by a known or
partially known graph. The graph hyperparameters `(mu, nu)` are selected
by minimizing the validation mean squared prediction error.

## Installation

`badgers` contains compiled C++ code, so a C++ toolchain is required
(Xcode command line tools on macOS, Rtools on Windows). You can install
`badgers` from GitHub with:

``` r
# install.packages("devtools")
devtools::install_github("zyhu888/badgers")
```

## Main Functions

- `badgers()` is the primary user-facing function. It fits a grid of
  `(mu, nu)` values and returns the run with the smallest validation
  MSPE.

- `badgers_mcmc()` runs the MCMC sampler for a fixed `(mu, nu)`.

- `dg()` simulates replicate measurements of high-dimensional
  predictors, a sparse response and a true or noisy working graph from
  user-specified parameters.

- `dg_batch()` generates and saves many data sets with `dg()`.

## Parallel Runs

`badgers()` evaluates the `(mu, nu)` grid serially. For large analyses,
run fixed-`(mu, nu)` jobs in parallel with `badgers_mcmc()` and combine
the MSPE results afterward. See the parallel-runs vignette and the
installed example:

``` r
system.file("examples", "parallel-runs.R", package = "badgers")
```

## Citation

References will be added after the manuscript is available.
