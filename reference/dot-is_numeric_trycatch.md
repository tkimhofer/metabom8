# Check if input is numeric-like using tryCatch

Tests whether a vector `x` can be safely coerced to numeric without
producing `NA`s (excluding original `NA`s) or warnings.

This function attempts to coerce the input to numeric and returns `TRUE`
if all non-missing elements convert without coercion failure; otherwise
`FALSE`. It catches warnings and errors during coercion and treats those
as non-numeric.

## Usage

``` r
.is_numeric_trycatch(x)
```

## Arguments

- x:

  A vector to test for numeric coercion.

## Value

Logical `TRUE` if `x` is numeric-like, `FALSE` otherwise.
