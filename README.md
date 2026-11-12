**badgers**: **B**ayesian **AD**aptive **G**raph-guided
**E**rrors-in-variables **R**egression **S**olution
================

[![License:
GPL-3](https://img.shields.io/badge/License-GPL--3-green.svg)](LICENSE.md)

## badgers

This R package fits BADGERS (**B**ayesian **AD**aptive **G**raph-guided
**E**rrors-in-variables **R**egression **S**olution), a graph-guided
Bayesian regression model for high-dimensional data with replicated
covariates and measurement error. Latent covariates are inferred from
the replicate measurements, which may differ in number across
predictors, and a working graph links the predictor-specific shrinkage
parameters through a constrained Wishart prior. The log-shrinkage mean
and scale `(mu, nu)` are selected by minimizing the validation mean
squared prediction error.

## Installation

`badgers` contains compiled C++ code, so a C++ toolchain is required
(Xcode command line tools on macOS, Rtools on Windows). You can install
`badgers` from a local source checkout with:

``` sh
R CMD INSTALL .
```

Once the GitHub repository is publicly accessible, install it with:

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
  user-specified parameters. Its defaults are the simulation settings of
  the paper.

- `dg_batch()` generates and saves many data sets with `dg()`.

## Data Format

Each of the training, validation and test sets is a response vector of
length `n` and a three-dimensional array of replicates with dimension
`num_rep` by `n` by `p`. Predictors may have different numbers of
replicates: put the replicates of predictor `j` in the first rows of
`W[, , j]` and pad the remaining rows with `NA`, using the same pattern
for every observation. The working graph is a `p` by `p` 0/1 adjacency
matrix.

## Reproducible Example

``` r
library(badgers)
sim <- dg(seed = 123, p = 8, g = 2, gsm = 2,
          qmin = 2, qmax = 3, ntrain = 12, ntune = 12, ntest = 12)

set.seed(456)
fit <- badgers(
  ytrain = sim$y, Wtrain = sim$W,
  ytune = sim$ytune, Wtune = sim$Wtune,
  ytest = sim$ytest, Wtest = sim$Wtest,
  adj_graph = sim$graph, mu_grid = c(1, 3), nu_grid = 1,
  n_iter = 100, burnin = 50, verbose = FALSE
)
print(fit)
#> BADGERS fit
#> Selected mu: 3
#> Selected nu: 1
#> Tuning MSPE: 4.363395
#> Test MSPE: 5.974529
#>
#> Nonzero coefficients by credible interval:
#>   95%: 0
#>   90%: 0
#>   80%: 1
#>
#> Tuning MSPE over the (mu, nu) grid:
#>    nu
#> mu         1
#>   1 7.912984
#>   3 4.363395
fit$best_result$beta_mean
#> [1]  0.53581721  0.53580288  0.02255819 -0.03233261  0.01734638  0.01843885
#> [7] -0.01275880 -0.01409319
which(fit$best_result$CI95 != 0)
#> integer(0)
```

This short chain demonstrates the interface. For inference, run longer,
independent chains and assess mixing, effective sample sizes, and
sensitivity. The defaults are 10000 iterations and 5000 burn-in; these
counts do not ensure convergence. Proposal adaptation continues after
burn-in.

`CI95`, `CI90`, and `CI80` contain posterior means set to zero when the
corresponding equal-tailed credible interval includes zero. They contain
coefficient estimates, not interval endpoints. Full coefficient draws,
including burn-in, are in `fit$best_result$BETA`.

The model has no intercept or automatic centering/scaling and assumes
zero-mean latent covariates. Estimate any preprocessing using training
data and apply it consistently to validation and test data. Missing
responses and non-finite observed values are rejected. Adjacency
matrices must be symmetric; the diagonal is ignored. Select parameters
using validation MSPE and reserve test outcomes for final evaluation.
Failed grid points produce a warning and are recorded in
`fit$grid_errors`.

See `vignette("badgers-introduction", package = "badgers")` for the
workflow.

## Parallel Runs

`badgers()` evaluates the `(mu, nu)` grid serially. For large analyses,
run fixed-`(mu, nu)` jobs in parallel with `badgers_mcmc()` and combine
the MSPE results afterward. The example uses the `parallel` package and
runs on Windows, macOS and Linux. See the parallel-runs vignette and the
installed example:

``` r
# path of the installed script
system.file("examples", "parallel-runs.R", package = "badgers")

# open it in the editor
file.edit(system.file("examples", "parallel-runs.R", package = "badgers"))

# or run it directly
source(system.file("examples", "parallel-runs.R", package = "badgers"))
```

The same script is also available in the repository as
`inst/examples/parallel-runs.R`, and it is shown in `?badgers`.

## Citation

If you use `badgers` in your research, please cite:

``` bibtex
@unpublished{hu2026graph,
  title  = {A Graph-Guided Bayesian Regression Model for High-Dimensional Data with Replicated Covariates and Measurement Error},
  author = {Hu, Zeyu and Qu, Kewen and Long, Qi and Yan, Jun and Li, Wenrui},
  note   = {Manuscript},
  year   = {2026}
}
```

## References

Hu, Z., Qu, K., Long, Q., Yan, J., & Li, W. (2026). A Graph-Guided
Bayesian Regression Model for High-Dimensional Data with Replicated
Covariates and Measurement Error. Manuscript.
