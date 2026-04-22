# Monte-Carlo cross-validation strategy

Monte-Carlo cross-validation strategy

## Usage

``` r
mc(k, split)
```

## Arguments

- k:

  Integer. Number of repeated random splits.

- split:

  Numeric. Fraction of samples assigned to the training set (e.g.
  `2/3`).

## Value

A named `list` with elements:

- train:

  List of integer vectors containing training set indices for each
  resampling iteration.

- strategy:

  Character string indicating the resampling strategy.

- n:

  Integer. Number of samples in the dataset.

- seed:

  Integer. Random seed used to generate the resampling splits, ensuring
  reproducibility.

## Details

Monte-Carlo cross-validation generates `k` random train/test splits
without replacement. No stratification is applied; samples are drawn
uniformly at random.

## See also

Other resampling strategies:
[`balanced_boot()`](https://tkimhofer.github.io/metabom8/reference/balanced_boot.md),
[`balanced_mc()`](https://tkimhofer.github.io/metabom8/reference/balanced_mc.md),
[`kfold()`](https://tkimhofer.github.io/metabom8/reference/kfold.md),
[`stratified_kfold()`](https://tkimhofer.github.io/metabom8/reference/stratified_kfold.md)

## Examples

``` r
n <- 100
# bivariate outcome
thr <- 1.5
Y <- c(rnorm(80, thr-3, 0.3), rnorm(20, thr+3, 0.3))  # unbalanced low/high outcome
mean(Y>thr)
#> [1] 0.2

cv_mc <- mc(k = 10, split = 2/3)
mc_inst <- metabom8:::.arg_check_cv(cv_pars=cv_mc, model_type='R', n=n, Y_prepped=cbind(Y))
sapply(mc_inst$train, function(i) length(i))
#>  [1] 66 66 66 66 66 66 66 66 66 66
```
