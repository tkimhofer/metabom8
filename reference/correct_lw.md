# Linewidth correction by scaling spectra to a reference linewidth

Applies a multiplicative correction to each spectrum to compensate for
linewidth-induced peak height variation, using a reference peak region
(e.g. TSP). The correction assumes an empirical power-law relationship
between peak height and linewidth:

## Usage

``` r
correct_lw(
  x,
  lw_ref = 0.8,
  beta = 0.6,
  shift = c(-0.1, 0.1),
  ppm = NULL,
  sf = NULL,
  estimate_beta = TRUE,
  beta_min_n = 6,
  only_if_improves = TRUE,
  notes = NULL
)
```

## Arguments

- x:

  A numeric matrix (samples x variables) or a named list containing `X`,
  `ppm`, and optionally `meta`.

- lw_ref:

  Numeric reference linewidth (FWHM) to which spectra are normalized. If
  `NULL`, the median linewidth across samples is used.

- beta:

  Numeric exponent governing the linewidth-to-height relationship.
  Ignored if `estimate_beta = TRUE` and sufficient samples are
  available.

- shift:

  Numeric vector of length 2 specifying the ppm region used for
  linewidth estimation and peak height measurement (e.g. `c(-0.1, 0.1)`
  for TSP).

- ppm:

  Optional numeric ppm axis. If `NULL`, inferred from
  `attr(X, "m8_axis")$ppm` or numeric `colnames(X)`.

- sf:

  Optional spectrometer frequency (MHz). If `NULL`, obtained from
  `meta$a_SFO1` in list input or `attr(X, "m8_meta")$a_SFO1`.

- estimate_beta:

  Logical; if `TRUE`, estimate `beta` from a log-log regression of
  reference peak height versus FWHM.

- beta_min_n:

  Minimum number of valid samples required to estimate `beta`.

- only_if_improves:

  Logical; if `TRUE`, apply correction only if the CV of the reference
  peak height decreases after correction.

- notes:

  Optional character string appended to the preprocessing log.

## Value

An object of the same type as input:

- If input is a matrix, returns a corrected matrix with attributes
  preserved and updated.

- If input is a list, returns the same list with corrected `$X`.

## Details

\$\$I \propto \mathrm{FWHM}^{-\beta}\$\$

Spectra are scaled such that all samples are adjusted to a common
reference linewidth `lw_ref`:

\$\$X\_{\mathrm{adj}} = X \cdot
\left(\frac{\mathrm{FWHM}}{lw\_{\mathrm{ref}}}\right)^{\beta}\$\$

where FWHM is estimated per spectrum within the specified `shift` region
(typically containing an internal reference signal).

If `estimate_beta = TRUE`, the exponent `beta` is estimated from a
log-log regression between peak height (max within `shift`) and FWHM
across samples. Otherwise, the user-supplied `beta` is used.

The correction is only applied if it reduces the coefficient of
variation (CV) of the reference peak height across samples when
`only_if_improves = TRUE`.

Input may be either:

- a numeric matrix `X` (samples x variables), with ppm axis supplied via
  `attr(X, "m8_axis")$ppm` or numeric `colnames(X)`, or

- a named list with elements `X`, `ppm`, and optionally `meta`.

Attributes `"m8_prep"`, `"m8_axis"`, and `"m8_meta"` are preserved and
the applied correction is appended to the preprocessing log via
`.m8_stamp()`.

This correction removes systematic peak height variation due to spectral
broadening (e.g. shimming differences) under the assumption that the
internal reference signal has constant concentration across samples. The
exponent `beta` reflects the empirical sensitivity of peak height to
linewidth in the current acquisition and processing regime.

This adjustment improves comparability of signal intensities across
spectra by reducing linewidth-induced amplitude bias, thereby enhancing
the detectability of small effects in downstream multivariate analyses
(e.g. PLS-type models).

## See also

[`lw`](https://tkimhofer.github.io/metabom8/reference/lw.md)

Other preprocessing:
[`align_segment()`](https://tkimhofer.github.io/metabom8/reference/align_segment.md),
[`align_spectra()`](https://tkimhofer.github.io/metabom8/reference/align_spectra.md),
[`binning()`](https://tkimhofer.github.io/metabom8/reference/binning.md),
[`calibrate()`](https://tkimhofer.github.io/metabom8/reference/calibrate.md),
[`correct_baseline()`](https://tkimhofer.github.io/metabom8/reference/correct_baseline.md),
[`pqn()`](https://tkimhofer.github.io/metabom8/reference/pqn.md),
[`print_preprocessing()`](https://tkimhofer.github.io/metabom8/reference/print_preprocessing.md)

## Examples

``` r
data(covid_raw)
# matrix input with ppm in m8_axis
Xcor <- correct_lw(covid_raw, shift = c(-0.1, 0.1))
```
