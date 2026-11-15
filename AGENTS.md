# Repository Guidelines

## Project Structure & Module Organization

`badgers` is an R package for Bayesian graph-guided regression with measurement error. `R/badgers.R` contains input validation, fitting wrappers, and simulation utilities; `R/print.R` defines fit printing. `src/badgers.cpp` implements the Rcpp/Armadillo sampler. `man/` contains function help, `vignettes/` contains tutorials, and `inst/examples/parallel-runs.R` demonstrates portable PSOCK parallel execution. Package metadata and exports live in `DESCRIPTION` and `NAMESPACE`. No dedicated test directory currently exists.

## Build, Test, and Development Commands

Use R >= 4.0.0 and a C++ toolchain, with the dependencies listed in `DESCRIPTION` installed. Run these commands from the repository root:

- `R CMD INSTALL .`: compile and install the package locally.
- `R CMD build .`: build the source archive, including vignettes.
- `R CMD check badgers_0.1.0.tar.gz`: check the built archive; adjust the filename when the version changes.
- `Rscript -e 'Rcpp::compileAttributes()'`: regenerate bindings after changing exported C++ signatures.
- `Rscript -e 'rmarkdown::render("README.Rmd")'`: regenerate `README.md` after editing its source.

## Coding Style & Naming Conventions

Follow the RStudio project settings: two-space indentation, spaces instead of tabs, UTF-8, Unix line endings, and no trailing whitespace. Match surrounding R and C++ style; use descriptive snake_case names for new functions and arguments. Preserve established public names and result fields such as `MSPE_tune` and `CI95`. No formatter or linter is configured.

Do not hand-edit `R/RcppExports.R` or `src/RcppExports.cpp`. Keep help files and examples synchronized with API changes. Preserve `-ffp-contract=off` in both Makevars files and the sampler's numerical accumulation behavior: these support reproducibility against R arithmetic.

## Testing Guidelines

No testing framework or coverage threshold is configured. Run package checks and exercise affected examples. For sampler changes, use reproducible seeds and small simulations; check input validation, output dimensions, credible intervals, and tuning MSPE. Cover unequal replicate counts and consistent trailing `NA` padding. Keep smoke runs short; research runs use longer chains.

## Commit & Pull Request Guidelines

Recent commits use concise imperative subjects, such as “Add GPL-3 license” and “Trim parallel example output.” Follow that convention and keep changes focused. In pull requests, describe the behavior changed, relevant issues, validation commands and results, and any numerical or API effects. Include documentation updates when user-facing behavior changes.
