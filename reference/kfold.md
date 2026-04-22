# K-fold cross-validation strategy

K-fold cross-validation strategy

## Usage

``` r
kfold(k)
```

## Arguments

- k:

  Integer number of folds.

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

Partitions the data into `k` folds. Each fold is used once as a test
set, with the remaining folds used for training. No stratification is
applied; folds are created by random partitioning.

## See also

Other resampling strategies:
[`balanced_boot()`](https://tkimhofer.github.io/metabom8/reference/balanced_boot.md),
[`balanced_mc()`](https://tkimhofer.github.io/metabom8/reference/balanced_mc.md),
[`mc()`](https://tkimhofer.github.io/metabom8/reference/mc.md),
[`stratified_kfold()`](https://tkimhofer.github.io/metabom8/reference/stratified_kfold.md)

## Examples

``` r
n <- 100
thr <- 1.5
Y <- c(rnorm(80, thr - 3, 0.3), rnorm(20, thr + 3, 0.3))  # unbalanced outcome
mean(Y > thr)
#> [1] 0.2
cv_k  <- kfold(k = 10)
k_inst <- metabom8:::.arg_check_cv(cv_pars=cv_k, model_type='R', n=n, Y_prepped=cbind(Y))
sapply(k_inst$train, function(i) length(i))
#>  [1] 90 90 90 90 90 90 90 90 90 90
```
