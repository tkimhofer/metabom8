# Ensure Input is a Row Matrix

Ensures that the input object `X` is in row-matrix form. If `X` is a
numeric vector (i.e., has no dimensions), it is converted to a matrix
with 1 row and `length(X)` columns. This is helpful for standardizing
inputs before applying matrix operations.

## Usage

``` r
.dimX(X)
```

## Arguments

- X:

  Numeric vector, matrix, or array. Typically NMR data where rows
  represent spectra.

## Value

A numeric matrix with one row if input was a vector; otherwise, returns
`X` unchanged.
