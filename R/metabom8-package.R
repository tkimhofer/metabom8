#' metabom8: A High-Performance R Package for Metabolomics Modeling and Analysis
#'
#' @description
#' metabom8 (pronounced \emph{metabo-mate}) provides pipelines for 1D NMR data import,
#' preprocessing, multivariate modeling (PCA, OPLS), metabolite identification,
#' and visualisation. Core functions are accelerated in C++ via \pkg{Rcpp},
#' \pkg{RcppArmadillo}, and \pkg{RcppEigen} for improved computational performance.
#'
#' @section Features:
#' \itemize{
#'   \item Import, preprocessing, and analysis of 1D NMR spectra.
#'   \item \strong{Principal Components Analysis (PCA)} with \emph{back-scaled} loadings
#'         (projected back to the spectral domain and visualized as spectra).
#'   \item \strong{Orthogonal Partial Least Squares (OPLS)} fitted iteratively via the NIPALS algorithm,
#'         with an \emph{automatic stopping criterion} to determine the optimal number of components.
#'   \item Automatic model selection using cross-validated performance (\eqn{R^2}, \eqn{Q^2},
#'         and cross-validated AUC for classification), with safeguards against overfitting.
#'   \item Robust statistical validation for small to large sample sizes:
#'         stratified Monte Carlo cross-validation and k-fold CV.
#'   \item Custom model visualizations using \pkg{ggplot2} and \pkg{plotly}.
#'   \item Model diagnostics and validation (DModX, permutation testing).
#'   \item Metabolite identification via STOCSY and STORM.
#'   \item Native C++ acceleration via \pkg{RcppArmadillo} and \pkg{RcppEigen}.
#' }
#'
#' @section Vignettes:
#' \itemize{
#'   \item Data import and preprocessing: \code{vignette("Preproc")}
#'   \item Multivariate analysis and metabolite identification: \code{vignette("MVA")}
#'   \item All vignettes: \code{browseVignettes("metabom8")}
#' }
#'
#' @name metabom8
"_PACKAGE"
