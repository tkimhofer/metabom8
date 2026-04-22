# Linear Interpolated Shift (Internal)

Applies a linear shift to a numeric vector using interpolation. Values
outside the original domain are extended using the boundary values (no
zero padding).

## Usage

``` r
.shift_interp(x, lag, ppm = NULL)
```

## Arguments

- x:

  Numeric vector.

- lag:

  Numeric scalar. Shift in index units. Positive values shift the signal
  to the right, negative values shift to the left.

## Value

Numeric vector of the same length as `x`.
