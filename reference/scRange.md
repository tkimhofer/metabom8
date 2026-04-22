# Min-Max Scaling to Arbitrary Range

Rescales a numeric vector to a specified range using min-max scaling.
This is a generalized form of min-max normalization allowing any output
range.

## Usage

``` r
scRange(x, ra)
```

## Arguments

- x:

  Numeric vector. Input values to be scaled.

- ra:

  Numeric vector of length 2. Desired output range (e.g., `c(5, 10)`).

## Value

A numeric vector of the same length as `x`, scaled to the range `ra`.

## Details

The scaled values are computed as: \$\$x\_{scaled} = \frac{x -
\min(x)}{\max(x) - \min(x)} \cdot (r\_{max} - r\_{min}) + r\_{min}\$\$

## See also

[`minmax()`](https://tkimhofer.github.io/metabom8/reference/minmax.md)

## Examples

``` r
x <- rnorm(20)
plot(x, type = 'l'); abline(h = range(x), lty = 2)
points(scRange(x, ra = c(5, 10)), type = 'l', col = 'red');
abline(h = c(5, 10), col = 'red', lty = 2)

```
