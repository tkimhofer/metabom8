# Generate Cross-Validation (CV) Training Set Indices

Generates a list of row indices representing training sets for various
CV strategies including stratified and Monte Carlo methods. Handles both
regression (R) and discriminant analysis (DA), and multi-column Y.

## Usage

``` r
.cvSetsMethod(
  Y,
  type,
  method = "k-fold_stratified",
  k = 7,
  split = 2/3,
  probs = NULL
)
```

## Arguments

- Y:

  matrix or data.frame. Outcome matrix. If multi-column, only the first
  column will be used.

- type:

  character. Type of analysis: "R" for regression or "DA" for
  discriminant analysis. Use suffix "-mY" to indicate multi-column Y.

- method:

  character. One of: `"k-fold"`, `"k-fold_stratified"`, `"MC"`,
  `"MC_balanced"`.

- k:

  integer. Number of CV iterations/folds.

- split:

  numeric. For Monte Carlo methods, the fraction of samples to use in
  the training set (0 \< split \< 1).

- probs:

  Numeric vector. Probabilities used for stratification of numeric `Y`
  (regression only).

## Value

A list of integer vectors. Each vector contains row indices representing
a training set.

## Details

- In `"k-fold_stratified"` and `"MC_balanced"` modes, stratification is
  applied to preserve class proportions or outcome distribution.

- If `type` is regression (`"R"`), quantile-based binning is applied to
  stratify continuous Y values.

- If `type` is discriminant analysis (`"DA"`), class labels are used
  directly.

- If `Y` has multiple columns, only the first column is used to create
  resampling sets.
