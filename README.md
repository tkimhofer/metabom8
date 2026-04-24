<!-- badges: start -->
  [![R-CMD-check](https://github.com/tkimhofer/metabom8/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/tkimhofer/metabom8/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

# metabom8 

**A High-Performance R Package for Metabolomics Modeling and Analysis**


`metabom8` (*metabo-mate*) provides:

- 1D NMR data import & preprocessing  
- Multivariate modeling (PCA, PLS, O-PLS)  
- Metabolite identification  

Core functionality is accelerated using C++ (`Rcpp`, `Armadillo`, `Eigen`) for improved performance.

---

## 🛠️ Features

- NMR data import and preprocessing, including: 
  - calibration, baseline correction, alignment, normalisation
  - provenance logging
- Multivariate statistics:  
  - Principal Components Analysis (PCA) 
  - Partial Least Squares (PLS) 
  - Orthogonal Partial Least Squares (OPLS)
- Automatic selection of the optimal number of components for supervised models
- Robust cross-validation framework, for large and small sample sizes (e.g., Monte Carlo CV)
- OPLS model diagnostics & validation (DModX, permutation testing)
- Tools supporting metabolite identification (STOCSY, STORM)
- Native C++ acceleration via `RcppArmadillo` and `RcppEigen`

---

## 📦 Installation

```r
# Install from Bioconductor (available in Bioc 3.23, release scheduled for 29/04/026)
install.packages("BiocManager")
BiocManager::install("metabom8")

# Install development version (until Bioconductor release)
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("tkimhofer/metabom8")
```

---

## 🚀 Quick Start

```r
library(metabom8)

# Load example data
data(hiit_raw, package = "metabom8")

# plot spectra interactively with plotly
plot_spec(hiit_raw)

# piped preprocessing 
hiit <- hiit_raw |> 
  calibrate(type = "tsp") |>
  excise() |>
  correct_baseline(method='asls') |>
  align_spectra() |>
  pqn()

## Provenance logging
print_provenance(hiit)


data(covid, package = "metabom8")
X <- covid$X
Y <- covid$an$type

# Modelling Context
uv = uv_scaling(center = TRUE)
mc_cv <- balanced_mc(k=15, split=2/3, type="DA")

# PCA
pca_model <- pca(covid$X, scaling = uv, ncomp=2)

# PLS
pls_model <- pls(X, Y, scaling = uv, validation_strategy=mc_cv)
show(pls_model)

# OPLS
opls_model <- opls(X, Y, scaling = uv, validation_strategy=mc_cv)

Tx <- scores(opls_model)
Px <- loadings(opls_model)
vip <- vip(opls_model)
```

---

## ⚡ Performance

Benchmarking OPLS modelling against another widely used R implementations under comparable conditions:

```r
bench::mark(
  metabom8 = metabom8::opls(X, Y, scaling=uv, validation_strategy=kfold),
  ropls    = ropls::opls(X, Y, predI = 1, orthoI = 1),
  check = FALSE
)
```

| Method     | Median Time | Iter/sec | Memory Use |
|------------|-------------|----------|------------|
| **metabom8** | 1.61 sec     | 0.617    | 247.33 MB  |
| **ropls**    | 9.12 sec     | 0.110    | 2.09 GB    |

> `metabom8` is based on C++ accelerated linear algebra and provides a >5× speed-up over `ropls`*.  

\**[`ropls`](https://bioconductor.org/packages/ropls) (uv scaling, 7-fold CV; 1 predictive and 1 orthogonal component).*

---

## 📘 Documentation

Comprehensive documentation and vignettes are available at:  
🔗 https://tkimhofer.github.io/metabom8/

---

## 📝 License & Citation

MIT License © Torben Kimhofer  
If `metabom8` or any of its components contributes to your work, please ensure appropriate citation. See [`CITATION`](https://github.com/tkimhofer/metabom8/blob/master/inst/CITATION) for citation details.

---

## 🙋 Getting Help

- Found a bug? Open an [issue](https://github.com/tkimhofer/metabom8/issues)
- Want to contribute? Fork the repo and submit a pull request

---

*Built with 💙 by [@tkimhofer](https://github.com/tkimhofer)*

