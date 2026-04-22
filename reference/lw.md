# Full Width at Half Maximum (FWHM) Estimation

Estimates the full width at half maximum (FWHM; line width) of a
singlet-like peak within a specified chemical-shift range for each
spectrum.

## Usage

``` r
lw(X, ppm = NULL, shift = c(-0.1, 0.1), sf)
```

## Arguments

- X:

  Numeric matrix (spectra in rows) *or* a named list as returned by
  [`read1d`](https://tkimhofer.github.io/metabom8/reference/read1d.md)/`read1d_proc`
  containing `X`, `ppm`, and `meta`.

- ppm:

  Numeric vector of chemical shift values (ppm) corresponding to columns
  of `X`. If `NULL`, `ppm` is inferred in the following order:

  1.  `attr(X, "m8_axis")$ppm` (if present),

  2.  numeric `colnames(X)` (if present).

- shift:

  Numeric vector of length 2. Chemical shift range containing the peak
  (e.g., `c(-0.1, 0.1)` for TSP).

- sf:

  Spectrometer frequency in MHz. Either a single numeric (recycled
  across spectra) or a numeric vector of length `nrow(X)` (one value per
  spectrum). If `X` is a list input and `sf` is missing, `sf` is taken
  from `X$meta$a_SFO1` when available. If `sf` is missing for a matrix
  input, the function attempts to use `attr(X, "m8_meta")$a_SFO1` when
  present.

## Value

Numeric vector of FWHM values in Hz (length `nrow(X)`).

## Details

For each spectrum, the function:

1.  extracts the region defined by `shift`,

2.  finds the peak apex within that region,

3.  computes the half-height level relative to the local baseline
    (minimum in the window),

4.  estimates the left and right half-height crossing points by linear
    interpolation,

5.  converts the width from ppm to Hz using `sf` (MHz), i.e.
    `Hz = ppm * sf`.

If no valid half-height crossings can be found (e.g., very low SNR or
truncated peak), `NA` is returned for that spectrum.

## Note

The ppm axis may be increasing or decreasing; FWHM is computed as an
absolute width and is therefore independent of axis direction.

## Examples

``` r
# Simulated NMR peaks with different linewidths
ppm <- seq(-0.2, 0.2, length.out = 1000)

# generate peaks with increasing width
sds <- seq(0.01, 0.03, length.out = 10)

X <- t(sapply(sds, function(s)
  dnorm(ppm, mean = 0, sd = s)
))

sf <- 600  # spectrometer frequency in MHz

fwhm_vals <- lw(X, ppm = ppm, shift = c(-0.1, 0.1), sf = sf)

plot(sds, fwhm_vals,
  xlab = "Gaussian sd",
  ylab = "Estimated FWHM (Hz)",
  pch = 16)
```
