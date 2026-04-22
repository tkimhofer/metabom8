# Plot STOCSY result

Generates a STOCSY plot (covariance trace coloured by absolute
correlation).

## Usage

``` r
plotStocsy(stoc_mod, shift = c(0, 10), title = NULL)
```

## Arguments

- stoc_mod:

  An object of class `m8_stocsy1d` returned by
  [`stocsy()`](https://tkimhofer.github.io/metabom8/reference/stocsy.md).

- shift:

  Numeric vector of length 2 specifying the chemical shift range (ppm).

- title:

  Optional character plot title.

## Value

A `ggplot2` object.

## Examples

``` r
# st <- stocsy(X, ppm, driver = 5.233, plotting = FALSE)
# plotStocsy(st, shift = c(5.15, 5.30), title = "Glucose")

data(covid)
cs = 5.233 # anomeric H of gluc
s1 <- stocsy(covid$X, driver=cs, plotting = FALSE)
plotStocsy(s1)

```
