
# metabom8 

**An R library for NMR-based metabolic profiling**  

`metabom8` (pronounced *metabo-mate*) provides pipelines for 1D NMR data import, preprocessing, multivariate modeling (PCA, O-PLS), metabolite identification, and visualization — with core functions accelerated using C++ (`Rcpp`, `Armadillo`, `Eigen`) for improved computational performance.

---

## 🛠️ Features

- Preprocessing and analysis of 1D NMR and MS spectral data
- Prinicpal Component Analysis (PCA) for unsupervised analysis 
- Orthogonal Partial Least Squares (OPLS) modeling with automatic selection of the optimal number of components
- Robust statistical validation: k-fold and stratified Monte Carlo CV
- Metabolite identification via STOCSY and STORM
- Custom plotting functions using `ggplot2` and `plotly`
- Native C++ acceleration via `RcppArmadillo` and `RcppEigen`

---

## 📦 Installation

```r
install.packages("remotes")
remotes::install_github("tkimhofer/metabom8")
```

---

## 🚀 Quick Start

```r
library(metabom8)
library(nmrdata)

# Load example data
data(bariatric, package = "nmrdata")

idx <- bariatric$an$Class %in% c("Pre-op", "RYGB")
X <- bariatric$X.pqn[idx, ]
Y <- bariatric$an$Class[idx]

# Fit an OPLS model using Monte Carlo Cross-Valdation
model <- opls(
  X = X,
  Y = Y,
  center = TRUE,
  scale = "UV",
  cv = list(method = "MC", k = 7, split=2/3)
)

# Plot scores and loadings
plotscores(model)
plotload(model)
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

- [MVA vignette](https://tkimhofer.github.io/metabom8/articles/MVA.html)
- [Preprocessing pipeline](https://tkimhofer.github.io/metabom8/articles/PreProc.html)
- [Function reference](https://tkimhofer.github.io/metabom8/reference/)

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

