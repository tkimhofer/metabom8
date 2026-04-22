# Check Y for missing values (PLS context)

Checks input Y for NA/NaN/Inf values, determines analysis type
(Regression or Discriminant Analysis), and returns cleaned Y matrix,
associated levels, and type.

## Usage

``` r
.checkYclassNas(Y)
```

## Arguments

- Y:

  A vector, matrix, or data.frame. Target variable.

## Value

A list: cleaned Y matrix, levels (if applicable), and type ('R' or 'DA')
