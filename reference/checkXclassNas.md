# Check X matrix (PLS context)

Validates that input `X` is a numeric matrix and does not contain
missing, NaN, or infinite values.

## Usage

``` r
.checkXclassNas(X)
```

## Arguments

- X:

  Input data (univariate or multivariate), formatted as a matrix.

## Value

`NULL` if all checks pass; otherwise throws an informative error.
