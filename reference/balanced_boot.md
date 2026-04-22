# Balanced bootstrap resampling strategy

Balanced bootstrap resampling strategy

## Usage

``` r
balanced_boot(k, split, type = c("DA", "R"), probs = NULL)
```

## Arguments

- k:

  Integer. Number of bootstrap resamples.

- split:

  Numeric. Fraction of samples drawn for the training set (e.g. `2/3`).
  Sampling is performed with replacement.

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

Generates `k` bootstrap samples (training sets sampled with
replacement). The remaining samples (the out-of-bag set) can be used as
a test set.

Balancing ensures equal representation of strata in the training data:

- `type = "DA"`:

  Class labels define the strata, and sampling is balanced across
  classes.

- `type = "R"`:

  The response is discretised into bins using quantiles defined by
  `probs`, and each bin contributes equally to the training set.

## See also

Other resampling strategies:
[`balanced_mc()`](https://tkimhofer.github.io/metabom8/reference/balanced_mc.md),
[`kfold()`](https://tkimhofer.github.io/metabom8/reference/kfold.md),
[`mc()`](https://tkimhofer.github.io/metabom8/reference/mc.md),
[`stratified_kfold()`](https://tkimhofer.github.io/metabom8/reference/stratified_kfold.md)

## Examples

``` r
n <- 100
# bivariate outcome
thr <- 1.5
Y <- c(rnorm(80, thr-3, 0.3), rnorm(20, thr+3, 0.3))  # unbalanced low/high outcome
mean(Y>thr)
#> [1] 0.2

cv_k <- kfold(k = 10)
cv_boot <- balanced_boot(k = 10, split = 2/3, type = "R", probs = c(0, 0.8, 1))

k_inst <- metabom8:::.arg_check_cv(cv_pars=cv_k, model_type='R', n=n, Y_prepped=cbind(Y))
b_inst <- metabom8:::.arg_check_cv(cv_pars=cv_boot, model_type='R', n=n, Y_prepped=cbind(Y))

# balanced splits: proportion above global median stays ~0.5
q80 <- quantile(Y, 0.8)
round(sapply(k_inst$train, function(i) mean(Y[i] > q80)), 2) # resembles original Y distr.
#>  [1] 0.20 0.22 0.21 0.22 0.17 0.18 0.21 0.20 0.20 0.19
round(sapply(b_inst$train, function(i) mean(Y[i] > q80)), 2) # balanced strata (low/high)
#>  [1] 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5
```
