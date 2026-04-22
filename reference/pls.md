# Fit a Partial Least Squares (PLS) model

Fits a supervised Partial Least Squares (PLS) model using a NIPALS-based
algorithm with optional cross-validation and automatic component
selection.

## Usage

``` r
pls(X, Y, scaling, validation_strategy, maxPCo = 5)
```

## Arguments

- X:

  Numeric matrix of predictors (rows = samples, columns = variables).

- Y:

  Numeric matrix or factor vector of responses.

- scaling:

  A scaling strategy object (e.g., `uv_scaling(center = TRUE)`),
  specifying model-internal centering and/or scaling applied during
  fitting. This does not modify the original spectral matrix.

- validation_strategy:

  A cross-validation strategy object defining how resampling is
  performed (e.g., k-fold, Monte Carlo).

- maxPCo:

  Integer. Maximum number of predictive components to evaluate.

## Value

An object of class `m8_model` containing the fitted PLS model,
cross-validation results, and performance statistics.

## Details

Model components are estimated sequentially using a NIPALS-based
algorithm. For each component, cross-validated performance metrics
(e.g., Q², R², classification AUC) are computed according to the
supplied `validation_strategy`. Component extraction stops when the
`stopRule` indicates overfitting or when `maxPCo` is reached.

Scaling specified via `scaling` is applied internally during model
fitting and does not alter the input matrix `X`. Spectral preprocessing
steps (e.g., alignment, baseline correction) should be performed prior
to model fitting.

The returned model object stores:

- Fitted component models

- Cross-validation results

- Performance metrics (R², Q², AUC)

- Model control parameters

- Input data provenance metadata

- Session information for reproducibility

## See also

[`opls`](https://tkimhofer.github.io/metabom8/reference/opls.md),
[`uv_scaling`](https://tkimhofer.github.io/metabom8/reference/uv_scaling.md)

Other modelling:
[`opls()`](https://tkimhofer.github.io/metabom8/reference/opls.md),
[`pca()`](https://tkimhofer.github.io/metabom8/reference/pca.md)

## Examples

``` r
data(covid)

cv <- balanced_mc(k=5, split=2/3)
scaling <- uv_scaling(center=TRUE)
model <-pls(X=covid$X, Y=covid$an$type, scaling, cv)
#> Performing discriminant analysis.
#> A PLS-DA model with 1 component was fitted.

show(model)
#> 
#> m8_model <pls>
#> ----------------------------------------
#> Dimensions  : 10 samples x 27819 variables
#> Mode        : classification 
#> Preprocess  : center | UV 
#> Components  : 1 (2 tested)
#> Validation  : BalancedMonteCarlo (k = 5) 
#> Stop rule   : cv_improvement_negligible 
#> ----------------------------------------
#> Use summary() for performance metrics.
#> 
summary(model)
#> $perf
#>   comp AUCt AUCv       R2X selected        fit               stop_reason
#> 1    1    1    1 0.2502977     TRUE     fitted                      <NA>
#> 2    2    1    1 0.2179150    FALSE not fitted cv_improvement_negligible
#>   stop_metric stop_delta
#> 1          NA         NA
#> 2           1          0
#> 
#> $engine
#> [1] "pls"
#> 
#> $y_type
#> [1] "DA"
#> 

Tp <- scores(model)
Pp <- loadings(model)
```
