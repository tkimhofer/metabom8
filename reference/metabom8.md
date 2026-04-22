# metabom8: A High-Performance R Package for Metabolomics Modeling and Analysis

metabom8 (pronounced *metabo-mate*) provides pipelines for 1D NMR data
import, preprocessing, multivariate modeling (PCA, OPLS), metabolite
identification, and visualisation. Core functions are accelerated in C++
via Rcpp, RcppArmadillo, and RcppEigen for improved computational
performance.

## Features

- Import, preprocessing, and analysis of 1D NMR spectra.

- **Principal Components Analysis (PCA)** with *back-scaled* loadings
  (projected back to the spectral domain and visualized as spectra).

- **Orthogonal Partial Least Squares (OPLS)** fitted iteratively via the
  NIPALS algorithm, with an *automatic stopping criterion* to determine
  the optimal number of components.

- Automatic model selection using cross-validated performance (\\R^2\\,
  \\Q^2\\, and cross-validated AUC for classification), with safeguards
  against overfitting.

- Robust statistical validation for small to large sample sizes:
  stratified Monte Carlo cross-validation and k-fold CV.

- Model diagnostics and validation (DModX, permutation testing).

- Metabolite identification via STOCSY and STORM.

- Native C++ acceleration via RcppArmadillo and RcppEigen.

## Vignettes

- Getting Started: `vignette("Getting Started")`

## See also

Useful links:

- <https://tkimhofer.github.io/metabom8/>

- Report bugs at <https://github.com/tkimhofer/metabom8/issues>

## Author

**Maintainer**: Torben Kimhofer <tkimhofer@gmail.com>
([ORCID](https://orcid.org/0000-0001-7158-9930))
