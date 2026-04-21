<!-- badges: start -->
  [![R-CMD-check](https://github.com/tkimhofer/metabom8/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/tkimhofer/metabom8/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

# metabom8 

**An R library for NMR-based metabolic profiling**  

`metabom8` (pronounced *metabo-mate*) provides pipelines for 1D NMR data import, preprocessing, multivariate modeling (PCA, O-PLS), metabolite identification, and visualization — with core functions accelerated using C++ (`Rcpp`, `Armadillo`, `Eigen`) for improved computational performance.

---

## 🛠️ Features

- Preprocessing and analysis of 1D NMR spectra
- Prinicpal Component Analysis (PCA) and backscaled model loadings
- Orthogonal Partial Least Squares (OPLS) modeling via the NIPALS algorithm
- OPLS with automatic selection of the optimal number of components based on R2, Q2 and cross-valdated AUC
- Robust statistical validation for small to large sample sizes: stratified Monte Carlo CV and k-fold CV
- Custom model visualisations using `ggplot2` and `plotly`
- Model diagnostics & validation (DModX, Permutation Testing)
- Metabolite identification via STOCSY and STORM
- Native C++ acceleration via `RcppArmadillo` and `RcppEigen`

---

## 📦 Installation

```r
# official release (to be added with bioconductor release)
install.packages("BiocManager")
BiocManager::install('metabom8')
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

Benchmarked against the `ropls` package under identical conditions using the `bench` package:

```r
bench::mark(
  metabom8 = metabom8::opls(X, Y, ...),
  ropls    = ropls::opls(X, Y, predI = 1, orthoI = 1, ...),
  check = FALSE
)
```

| Method     | Median Time | Iter/sec | Memory Use |
|------------|-------------|----------|------------|
| **metabom8** | 1.61 sec     | 0.621    | 2.81 GB     |
| **ropls**    | 3.95 sec     | 0.253    | 1.61 GB     |

> `metabom8` provides a ~2.5× speed-up using C++ linear algebra (Eigen, Armadillo) when compared to widely used implementations [^1].

[^1]: Benchmark comparison uses the [`ropls`](https://bioconductor.org/packages/ropls) package (Bioconductor), a widely used OPLS implementation (parameters: uv scaling, 7-fold CV, 1 predictive + 1 orthogonal component).
---

## 📘 Documentation

Comprehensive documentation and vignettes are available at:  
🔗 https://tkimhofer.github.io/metabom8/

---

## 🔗 Related Packages

- [`nmrdata`](https://github.com/tkimhofer/nmrdata): Example dataset for NMR spectral analysis (used for performance comparison)
  
If you find `metabom8` useful, please consider giving it a ⭐ — it makes it easier for others to discover the project!

---

## 📝 License & Citation

MIT License © Torben Kimhofer  
If `metabom8` or any of its components contributes to your work, please ensure appropriate citation. See [`CITATION`](https://github.com/tkimhofer/metabom8/blob/master/inst/CITATION) for citation details.

---

## 🙋 Getting Help

- Found a bug? Open an [issue](https://github.com/tkimhofer/metabom8/issues)
- Want to contribute? Fork the repo and submit a pull request
- Code of Conduct: see [`CODE_OF_CONDUCT.md`](https://github.com/tkimhofer/metabom8/blob/master/CODE_OF_CONDUCT.md)

---

*Built with 💙 by [@tkimhofer](https://github.com/tkimhofer)*

