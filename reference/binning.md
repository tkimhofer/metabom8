# Spectral data binning

Equidistant binning of spectra by summarising intensities within ppm
bins.

## Usage

``` r
binning(X, ppm = NULL, width = NULL, npoints = NULL, fun = sum)
```

## Arguments

- X:

  Numeric matrix or data frame with spectra in rows, or a named list as
  returned by
  [`read1d`](https://tkimhofer.github.io/metabom8/reference/read1d.md)/`read1d_proc`
  containing `X` and `ppm`.

- ppm:

  Numeric vector of chemical shift positions (length must match
  `ncol(X)`). If `NULL`, `ppm` is inferred in the following order:

  1.  `X$ppm` if `X` is a list input,

  2.  `attr(X, "m8_axis")$ppm` (if present),

  3.  numeric `colnames(X)` (if present).

- width:

  Numeric. Bin size in ppm, or `NULL` if `npoints` is specified.

- npoints:

  Integer. Desired number of bins per spectrum, or `NULL` if `width` is
  specified. If both are provided, `npoints` is used.

- fun:

  Function. Summary function applied to intensities within each bin.
  Must return a single numeric value (e.g. `sum`, `mean`, `max`).

## Value

Numeric matrix with spectra in rows and binned ppm variables in columns.

## Details

If present, preprocessing provenance is appended to `attr(X, "m8_prep")`
using `.m8_stamp()`. The ppm axis is updated in `attr(X, "m8_axis")$ppm`
and column names are set to the bin centres.

When `width` is specified, spectra are interpolated onto a regular ppm
grid and then aggregated within bins (`interp = TRUE` in provenance).
When `npoints` is specified, aggregation is performed by index bins on
the original grid (`interp = FALSE` in provenance).

## See also

Other preprocessing:
[`align_segment()`](https://tkimhofer.github.io/metabom8/reference/align_segment.md),
[`align_spectra()`](https://tkimhofer.github.io/metabom8/reference/align_spectra.md),
[`calibrate()`](https://tkimhofer.github.io/metabom8/reference/calibrate.md),
[`correct_baseline()`](https://tkimhofer.github.io/metabom8/reference/correct_baseline.md),
[`correct_lw()`](https://tkimhofer.github.io/metabom8/reference/correct_lw.md),
[`pqn()`](https://tkimhofer.github.io/metabom8/reference/pqn.md),
[`print_preprocessing()`](https://tkimhofer.github.io/metabom8/reference/print_preprocessing.md)

## Examples

``` r
set.seed(1)
X <- matrix(rnorm(2 * 100), nrow = 2)
ppm <- round(seq(10, 0.5, length.out = 100), 3)
colnames(X) <- ppm

Xb <- binning(X, ppm, width = 0.5)
Xb_mean <- binning(X, ppm, width = 0.5, fun = mean)
dim(Xb)
#> [1]  2 19
```
