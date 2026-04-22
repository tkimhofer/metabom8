# Check Dimensions of X and Y Matrices (PLS Context)

Checks that the input matrices `X` and `Y` have compatible dimensions
for Partial Least Squares (PLS) analysis. Specifically, the number of
observations (rows) must match, and the number of variables (columns) in
`X` must be greater than in `Y`.

## Usage

``` r
.checkDimXY(X, Y)
```

## Arguments

- X:

  Numeric matrix. Predictor data matrix (observations x variables).

- Y:

  Numeric matrix. Response data matrix (observations x variables).

## Value

Invisibly returns `NULL` if checks pass. Throws an error if dimensions
are incompatible.
