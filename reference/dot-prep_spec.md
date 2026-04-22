# Prepare Spectral Data for Plotting

Internal helper that validates input, enforces ppm ordering, subsets the
selected shift region, and returns matrix-form data.

## Usage

``` r
.prep_spec(x, ppm, shift)
```

## Arguments

- x:

  Numeric vector or matrix of spectra.

- ppm:

  Numeric chemical shift vector.

- shift:

  Length-2 numeric vector specifying ppm window.

## Value

A list with elements:

- ppm:

  Subsetted chemical shift axis

- X:

  Matrix of subsetted spectra
