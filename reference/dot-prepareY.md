# Prepare Response Vector for OPLS/OPLS-DA

Converts a response vector `Y` into a numeric matrix suitable for
regression or classification models.

- Numeric input is returned as a column matrix.

- Factor or character input:

  - If binary (2 levels): returns a single numeric column (values 1 and
    2).

  - If multiclass (\>2 levels): returns a dummy matrix with 1 for
    presence and -1 for absence. Also returns a mapping between original
    labels and numeric values (only for categorical input).

## Usage

``` r
.prepareY(Y)
```

## Arguments

- Y:

  A numeric, factor, or character vector.

## Value

A list with:

- `[[1]]`:

  Numeric matrix of processed response.

- `[[2]]`:

  Data frame mapping class labels to encoding (empty for numeric input).
