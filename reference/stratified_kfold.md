# Y-stratified k-fold cross-validation strategy

Y-stratified k-fold cross-validation strategy

## Usage

``` r
stratified_kfold(k, type = c("DA", "R"), probs = NULL)
```

## Arguments

- k:

  Integer. Number of folds.

- type:

  Character. Either `"DA"` (classification) or `"R"` (regression).

- probs:

  Numeric vector of quantile probabilities used to stratify continuous
  `Y` when `type = "R"`.

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

For classification (`type = "DA"`), folds are generated such that class
proportions are approximately preserved in each fold.

For regression (`type = "R"`), `Y` is discretised into bins defined by
`probs`, and folds are generated to approximately preserve the bin
distribution.

## See also

Other resampling strategies:
[`balanced_boot()`](https://tkimhofer.github.io/metabom8/reference/balanced_boot.md),
[`balanced_mc()`](https://tkimhofer.github.io/metabom8/reference/balanced_mc.md),
[`kfold()`](https://tkimhofer.github.io/metabom8/reference/kfold.md),
[`mc()`](https://tkimhofer.github.io/metabom8/reference/mc.md)

## Examples

``` r
set.seed(1)
n <- 100
thr <- 1.5
Y <- c(rnorm(80, thr - 3, 0.3), rnorm(20, thr + 3, 0.3))  # unbalanced outcome
mean(Y > thr)
#> [1] 0.2

q80 <- quantile(Y, 0.8)  # defines the rare "high" stratum (top 20%)

cv_k  <- kfold(k = 10)
cv_sk  <- stratified_kfold(k = 10, type = "R", probs = c(0, 0.8, 1))

k_inst <-  metabom8:::.arg_check_cv(cv_pars=cv_k, model_type='R', n=n, Y_prepped=cbind(Y))
sk_inst <-  metabom8:::.arg_check_cv(cv_pars=cv_sk, model_type='R', n=n, Y_prepped=cbind(Y))

round(sapply(k_inst$train,  function(i) mean(Y[i] > q80)), 2)  # reflects imbalance
#>  [1] 0.19 0.18 0.20 0.22 0.21 0.21 0.20 0.20 0.20 0.19
round(sapply(sk_inst$train, function(i) mean(Y[i] > q80)), 2)  # balanced across strata
#>  [1] 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2
```
