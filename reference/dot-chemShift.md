# Calculate Chemical Shift Axis

Computes a 1D NMR chemical shift axis based on sweep width, offset, and
size.

## Usage

``` r
.chemShift(swidth, offset, si)
```

## Arguments

- swidth:

  Numeric. Sweep width (Hz).

- offset:

  Numeric. Frequency offset (ppm).

- si:

  Integer. Number of data points.

## Value

Numeric vector of chemical shift values (ppm).
