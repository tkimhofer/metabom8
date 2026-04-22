# Generate Monte Carlo Cross-Validation (MCCV) Training Indices

Generates a list of training set indices for Monte Carlo
Cross-Validation.

## Usage

``` r
.mc(k, Y, split)
```

## Arguments

- k:

  Integer. Number of MCCV iterations (i.e., training sets to generate).

- Y:

  A matrix (n x p), where rows are observations and columns are
  outcomes. Only the number of rows is used here.

- split:

  Numeric. Fraction of the data to include in each training set. Should
  be between 0 and 1.

## Value

A list of length `k`, each containing a vector of training set indices.

## Details

For each fold, a random sample (without replacement) of size
`floor(nrow(Y) * split)` is drawn.
