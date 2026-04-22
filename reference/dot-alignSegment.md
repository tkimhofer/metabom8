# Align NMR Spectra via Cross-Correlation

Aligns rows of an NMR spectral matrix to a reference spectrum by
maximizing cross-correlation. Optionally uses the row-wise median
spectrum as reference. Useful for minor spectral misalignments.

## Usage

``` r
.alignSegment(
  seg,
  idx_ref = 1,
  clim = 0.7,
  med = TRUE,
  norm = FALSE,
  lag.max = 20
)
```

## Arguments

- seg:

  Numeric matrix. Each row is a 1D NMR spectrum segment.

- idx_ref:

  Integer. Row index to use as reference spectrum. Ignored if
  `med = TRUE`.

- clim:

  Numeric. Minimum cross-correlation threshold. Segments with lower
  similarity are not shifted.

- med:

  Logical. If `TRUE`, the row-wise median spectrum is used as the
  reference.

- norm:

  Logical. If `TRUE`, rows are z-scaled prior to alignment.

## Value

A numeric matrix of aligned spectra with the same dimensions as `seg`.
If `med = TRUE`, the output excludes the median row.

## Details

Each spectrum is aligned to the reference by computing the
cross-correlation and applying a lag-based linear interpolated shift.
The shift is applied only if the cross-correlation exceeds `clim`. This
function is intended for fine-tuning alignment across short ppm windows
(i.e., segments).
