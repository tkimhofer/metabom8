# Stratified k-fold cross-validation index generator

Generates a list of training set indices for stratified k-fold
cross-validation (CV). Stratification is performed based on the first
column of `Y`. For regression tasks, `Y` is binned by quantiles to
emulate class balance. Ensures class proportions are preserved across
folds when possible.

## Usage

``` r
.kFoldStratified(k, stratified)
```

## Arguments

- k:

  Integer. Number of folds.

- stratified:

  List with three elements:

  - `type`: Character. Either `"R"` for regression or `"DA"` for
    discriminant analysis.

  - `Y`: Matrix. Outcome matrix with a single column.

  - `probs`: Numeric vector. Probabilities used for stratification of
    numeric `Y` (regression only).

## Value

A list of length `k`, each element containing the training set row
indices (integers) for one CV fold. Returns `NULL` and a warning if
stratification is not feasible due to class imbalance.
