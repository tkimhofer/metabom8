# Dynamic Interval Segmentation (RSPA-style)

Defines peak-system intervals from the median spectrum for recursive
segment-wise alignment (RSPA).

## Usage

``` r
.dynamicIntervalsMedian(X, ppm, half_win_ppm = 0.02, min_snr = 8)
```

## Arguments

- X:

  Numeric matrix (spectra in rows).

- ppm:

  Numeric ppm vector.

- half_win_ppm:

  Numeric. Half window size around detected peaks.

## Value

List of index vectors defining alignment intervals.

## See also

Other NMR:
[`get_idx()`](https://tkimhofer.github.io/metabom8/reference/get_idx.md)
