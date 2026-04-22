# Check Consistency of Spectral Data and ppm Vector

Validates that the NMR data matrix and the chemical shift reference
vector (`ppm`) are consistent. Specifically, it checks that

- `ppm` contains no NA or infinite values.

- The number of columns in `X` matches the length of `ppm`.

This function is typically used as a preliminary data validation step
before downstream spectral analysis.

## Usage

``` r
.check_X_ppm(X, ppm)
```

## Arguments

- X:

  Numeric matrix. NMR spectral matrix with samples in rows and variables
  in columns (e.g., intensities).

- ppm:

  Numeric vector. Chemical shift reference values corresponding to
  columns of `X`.

## Value

Logical. Returns `TRUE` if inputs are valid, `FALSE` otherwise.
