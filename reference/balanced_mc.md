# Balanced Monte-Carlo resampling strategy

Balanced Monte-Carlo resampling strategy

## Usage

``` r
balanced_mc(k, split, type = c("DA", "R"), probs = NULL)
```

## Arguments

- k:

  Integer. Number of repeated random splits.

- split:

  Numeric. Fraction of samples assigned to the training set (e.g.
  `2/3`).

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

Generates `k` Monte-Carlo resampling splits by randomly partitioning the
data into training and test sets without replacement.

Balancing ensures equal representation of strata in the training data:

- `type = "DA"`:

  Class labels define the strata, and sampling is balanced across
  classes.

- `type = "R"`:

  The response is discretised into bins using quantiles defined by
  `probs`, and each bin contributes equally to the training set.

This strategy can improve robustness of model evaluation in settings
with limited samples size and imbalanced or unevenly distributed outcome
variables.

## See also

Other resampling strategies:
[`balanced_boot()`](https://tkimhofer.github.io/metabom8/reference/balanced_boot.md),
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
cv_mc <- balanced_mc(k = 10, split = 2/3, type = "R", probs = c(0, 0.8, 1))

k_inst <- metabom8:::.arg_check_cv(cv_pars=cv_k, model_type='R', n=n, Y_prepped=cbind(Y))
mc_inst <- metabom8:::.arg_check_cv(cv_pars=cv_mc, model_type='R', n=n, Y_prepped=cbind(Y))

# balanced splits: proportion above global median stays ~0.5
q80 <- quantile(Y, 0.8)
round(sapply(k_inst$train, function(i) mean(Y[i] > q80)), 2) # resembles original Y distr.
#>  [1] 0.19 0.21 0.21 0.20 0.20 0.19 0.18 0.22 0.19 0.21
round(sapply(mc_inst$train, function(i) mean(Y[i] > q80)), 2) # balanced strata (low/high)
#>  [1] 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5
```
