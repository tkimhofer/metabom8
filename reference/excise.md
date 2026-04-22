# Excise Chemical Shift Regions from 1D NMR Spectra

Removes specified chemical shift regions from 1D \\^1\\H NMR spectra. By
default, commonly excluded metabolomics regions are removed
(upfield/downfield noise, water, urea).

## Usage

``` r
excise(x, ppm = NULL, regions = NULL)
```

## Arguments

- x:

  Numeric matrix or vector. Spectra in rows and variables in columns.

- ppm:

  Numeric vector of chemical shift positions (ppm). If omitted, `ppm` is
  inferred from `colnames(X)`.

- regions:

  Named list of numeric vectors (length 2), specifying ppm regions to
  remove. Each element must define lower and upper bounds.

## Value

Numeric matrix with specified chemical shift regions removed.

## Details

Default regions removed (ppm):

- Upfield noise: `[min(ppm), 0.25]`

- Residual water: `[4.5, 5.2]`

- Urea region: `[5.5, 6.0]`

- Downfield noise: `[9.7, max(ppm)]`

Removed regions are recorded in `attr(X, "m8_prep")` if present. Updated
ppm values are stored in column names and in `attr(X, "m8_axis")$ppm`.

## Examples

``` r
set.seed(1)
ppm <- seq(0, 10, length.out = 1000)
X <- matrix(rnorm(100 * length(ppm)), nrow = 100)

Xe <- excise(X, ppm)

dim(Xe)
#> [1] 100 825
names(attributes(Xe))
#> [1] "dim"      "dimnames" "m8_axis"  "m8_prep" 
```
