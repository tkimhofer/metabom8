# Subset Optimisation by Reference Matching (STORM)

Selects an optimal subset of spectra that best match a specified target
signal region, improving downstream correlation-based structural
analysis such as STOCSY.

## Usage

``` r
storm(X, ppm, b = 30, q = 0.05, idx.refSpec, shift)
```

## Arguments

- X:

  Numeric matrix (or data.frame) of NMR spectra with samples in rows and
  spectral variables in columns.

- ppm:

  Numeric vector of chemical shift values corresponding to the columns
  of `X`.

- b:

  Integer. Half-window size expressed as number of spectral variables
  (data points). The effective window width therefore depends on the ppm
  spacing.

- q:

  Numeric. P-value threshold for including spectral variables when
  updating the reference region.

- idx.refSpec:

  Integer. Row index of `X` defining the initial reference spectrum.

- shift:

  Numeric vector of length 2 giving the ppm range (minimum, maximum)
  defining the initial target signal region.

## Value

Integer vector of row indices of `X` defining the selected spectral
subset.

## Details

STORM iteratively refines a subset of spectra exhibiting consistent
signal position and multiplicity within a specified ppm region.

Starting from an initial reference spectrum and ppm window:

1.  Spectra positively correlated with the current reference signal are
    retained.

2.  A driver peak (maximum intensity within the reference window) is
    identified.

3.  Correlation and covariance are evaluated within a local window of
    size `b` around the driver peak.

4.  The reference region is updated using variables that satisfy the
    p-value threshold (`q`) and show positive correlation.

The procedure continues until the selected subset stabilises. The
resulting row indices define spectra that most consistently represent
the structural pattern of the target signal.

STORM does not perform metabolite identification directly. Instead, it
refines the dataset to enhance structural coherence prior to
correlation-based interpretation methods such as STOCSY.

## References

Posma, J. M., et al. (2012). Subset optimisation by reference matching
(STORM): an optimised statistical approach for recovery of metabolic
biomarker structural information from \\^1H\\ NMR spectra of biofluids.
*Analytical Chemistry*, 84(24), 10694–10701.

## See also

Other structural_annotation:
[`stocsy()`](https://tkimhofer.github.io/metabom8/reference/stocsy.md)

## Examples

``` r
## Simulated example with three Gaussian signals
set.seed(123)

n <- 100          # number of spectra
S <- 1000         # number of spectral variables
ppm <- seq(10, 0, length.out = S)

gauss <- function(x, centre, width, height) {
  height * exp(-((x - centre)^2) / (2 * width^2))
}

## Generate base signals
sig1 <- gauss(ppm, 7,   0.05, 10)
sig2 <- gauss(ppm, 3.5, 0.10,  8)
sig3 <- gauss(ppm, 1,   0.07,  5)

spectra <- matrix(0, n, S)

for (i in seq_len(n)) {
  spectra[i, ] <- sig1 + sig2 + rnorm(S, 0, 0.1)
  if (i <= 25) {
    spectra[i, ] <- spectra[i, ] + sig3
  }
}

## Apply STORM to refine spectra containing the 1 ppm signal
idx <- storm(
  X = spectra,
  ppm = ppm,
  b = 30,
  q = 0.05,
  idx.refSpec = 1,
  shift = c(0.75, 1.25)
)

length(idx)   # number of spectra retained
#> [1] 25
```
