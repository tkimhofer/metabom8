# Cliff's Delta Effect Size

Calculates Cliff's delta (`Cd`) as a non-parametric effect size for
comparing two numeric vectors. `Cd` quantifies the directional
difference between a reference (`ref`) and comparison (`comp`)
distribution based on pairwise comparisons of their values.

`Cd` ranges from -1 to 1, where:

- -1: values in `comp` are consistently larger than those in `ref`

- 0: both distributions are similar, with no systematic difference

- 1: values in `ref` are consistently larger than those in `comp`

## Usage

``` r
cliffs_d(ref, comp)

es_cdelta(ref, comp)
```

## Arguments

- ref:

  Numeric vector representing the reference group.

- comp:

  Numeric vector representing the comparator group.

## Value

A single numeric value: Cliff's delta effect size.

## Details

The effect size is calculated as the scaled difference in dominance
between groups. Missing or infinite values are removed with a message.
Synonyms are

## References

Cliff, N. (1993). Dominance statistics: Ordinal analyses to answer
ordinal questions. *Psychological Bulletin*, 114(3), 494–509.
[doi:10.1037/0033-2909.114.3.494](https://doi.org/10.1037/0033-2909.114.3.494)

## Examples

``` r
ref <- rnorm(100, mean = 0)
comp <- rnorm(100, mean = 1)

hist(ref,  col=rgb(0,0,1,0.3), breaks=50, xlim=c(-4,4))
hist(comp, col=rgb(1,0,0,0.3), breaks=50, add=TRUE)
legend("topright", legend=c("ref","comp"), fill=c("blue","red"))


cliffs_d(ref, comp)
#> [1] -0.532
cliffs_d(comp, ref)
#> [1] 0.532

comp1 <- rnorm(100, mean = 10)
cliffs_d(ref, comp1)
#> [1] -1

cliffs_d(ref, ref)
#> [1] 0
```
