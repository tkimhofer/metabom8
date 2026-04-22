# Compute Variable Importance in Projection (VIP)

Calculates VIP scores for a PLS/OPLS model using predictive components
only.

## Usage

``` r
.vip(Tx, W, C)
```

## Arguments

- Tx:

  Numeric matrix (n x A). Score matrix.

- W:

  Numeric matrix (A x p). Weight matrix. Rows correspond to components,
  columns to variables.

- C:

  Numeric vector or matrix of length A. Y-loadings.

## Value

Numeric vector of length p containing VIP scores.

## Details

VIP is computed as:

\$\$ VIP_j = \sqrt{ p \frac{\sum_a SSY_a w\_{aj}^2}{\sum_a SSY_a} } \$\$

where SSY_a = (t_a^T t_a) \* c_a^2

By construction, mean(VIP^2) = 1.
