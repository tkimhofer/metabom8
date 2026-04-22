# Min-Max Scaling to \\\[0,1\]\\

Scales a numeric vector to the range \\\[0,1\]\\ using min-max
normalization. This is a special case of
[`scRange`](https://tkimhofer.github.io/metabom8/reference/scRange.md).

## Usage

``` r
minmax(x, na.rm = FALSE)
```

## Arguments

- x:

  Numeric vector. Input values to be scaled.

- na.rm:

  Logical; if `TRUE`, ignore `NA`s when computing the range.

## Value

A numeric vector of the same length as `x`, scaled to the range
\\\[0,1\]\\.

## Details

The scaled values are computed as: \$\$x\_{scaled} = \frac{x -
\min(x)}{\max(x) - \min(x)}\$\$

Equivalent to `scRange(x, ra = c(0, 1))`.

## See also

[`scRange()`](https://tkimhofer.github.io/metabom8/reference/scRange.md)
for flexible output ranges.

## Examples

``` r
x <- rnorm(20)
plot(x, type = 'l'); abline(h = range(x), lty = 2)
points(minmax(x), type = 'l', col = 'blue')
abline(h = c(0, 1), col = 'blue', lty = 2)

```
