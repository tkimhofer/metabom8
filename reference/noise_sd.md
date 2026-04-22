# Estimate Noise Standard Deviation in 1D NMR Spectra

Estimates the noise standard deviation (\\\sigma\\) for each spectrum
from a signal-free ppm region. The default estimator is robust
(MAD-based) and suitable for signal-to-noise calculations.

## Usage

``` r
noise_sd(
  X,
  ppm = NULL,
  where = c(14.6, 14.7),
  method = c("mad", "sd", "p95"),
  baseline_correct = FALSE,
  lambda = 10000,
  min_points = 50L,
  ns = NULL,
  normalise_scans = FALSE
)
```

## Arguments

- X:

  Numeric matrix (spectra in rows), numeric vector (single spectrum), or
  a named list as returned by
  [`read1d`](https://tkimhofer.github.io/metabom8/reference/read1d.md)/`read1d_proc`
  containing `X`, `ppm`, and `meta`.

- ppm:

  Numeric vector of chemical shift values (ppm) corresponding to columns
  of `X`. If `NULL`, `ppm` is inferred in the following order:

  1.  `X$ppm` if `X` is a list input,

  2.  `attr(X, "m8_axis")$ppm` (if present),

  3.  numeric `colnames(X)` (if present).

- where:

  Numeric vector of length 2. ppm range used for noise estimation.
  Should be free of metabolite signals (default `c(14.6, 14.7)`).

- method:

  Character. Noise estimator: `"mad"` (default), `"sd"`, or `"p95"`
  (legacy amplitude).

- baseline_correct:

  Logical. If `TRUE`, subtract a smooth baseline in the noise window
  using asymmetric least squares
  ([`asysm`](https://rdrr.io/pkg/ptw/man/asysm.html)). Default `FALSE`.

- lambda:

  Numeric. Smoothing parameter for
  [`asysm`](https://rdrr.io/pkg/ptw/man/asysm.html) when
  `baseline_correct = TRUE`.

- min_points:

  Integer. Minimum number of points required in the noise window. May be
  a single value (recycled across spectra) or a numeric vector of length
  `nrow(X)`. If `X` is a list input as returned by
  [`read1d`](https://tkimhofer.github.io/metabom8/reference/read1d.md)/`read1d_proc`
  and `ns` is `NULL`, the function attempts to extract the number of
  scans from `X$meta$a_NS` when available.

- ns:

  Number of scans / transients (required if normalise_scans=TRUE)

- normalise_scans:

  Logical. If `TRUE`, noise estimates are multiplied by \\\sqrt{NS}\\ to
  account for the theoretical scaling of noise with the number of scans
  (\\\sigma \propto 1/\sqrt{NS}\\). This is useful when comparing noise
  levels across spectra acquired with different numbers of scans.
  Default is `FALSE`.

## Value

Numeric vector of noise estimates (length `nrow(X)`).

## Details

In NMR spectroscopy, noise scales predictably with the number of scans
(\\NS\\). For otherwise identical acquisition settings:

- Signal increases approximately proportional to \\NS\\.

- Noise increases approximately proportional to \\\sqrt{NS}\\.

- Consequently, signal-to-noise ratio (SNR) increases proportional to
  \\\sqrt{NS}\\.

Equivalently, the noise standard deviation scales as:

\$\$\sigma(NS) \propto \frac{1}{\sqrt{NS}},\$\$

assuming a fixed underlying signal scale and comparable acquisition
conditions.

To compare noise levels across datasets acquired with different numbers
of scans, a scan-normalised noise estimate may be used:

\$\$\sigma\_{\mathrm{norm}} = \sigma \cdot \sqrt{NS}.\$\$

Under stable receiver gain and processing conditions, this normalised
noise should be approximately constant across runs.

## Examples

``` r
data(hiit_raw)
X <- hiit_raw$X
ppm <- hiit_raw$ppm
sigma <- noise_sd(X, ppm, where = c(10,11))
plot(hiit_raw$meta$a_NS, sigma,
  xlab = "Number of scans (NS)",
  ylab = expression(sigma~"(noise estimate)"),
   pch = 16)
 lines(lowess(hiit_raw$meta$a_NS, sigma), col = "red", lwd = 2)
```
