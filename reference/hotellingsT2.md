# Hotelling T^2 Statistic

Computes the Hotelling T^2 statistic defining a multivariate confidence
region for a set of observations. The region corresponds to an ellipsoid
in p-dimensional space and is commonly visualised as an ellipse in
two-dimensional (OPLS) score plots.

## Usage

``` r
hotellingsT2(X, alpha = 0.95)
```

## Arguments

- X:

  Numeric matrix with observations in rows and variables (dimensions) in
  columns.

- alpha:

  Numeric scalar. Confidence level for the T^2 region. Default is 0.95.

## Value

A named `list` with elements:

- center:

  Numeric vector of column means.

- cov:

  Sample covariance matrix.

- c2:

  Squared radius of the Hotelling T^2 region.

## Details

The Hotelling T^2 region is defined as \$\$(x - \mu)^T S^{-1} (x - \mu)
\le c^2\$\$ where \\c^2 = (p (n - 1) / (n - p)) F\_{p, n-p}(\alpha)\\.

## Examples

``` r
set.seed(1)
X <- matrix(rnorm(200), ncol = 2)
t2 <- hotellingsT2(X)
ell <- ellipse2d(t2)
plot(X, asp = 1)
lines(ell$x, ell$y, col = "red", lwd = 2)
```
