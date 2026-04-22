# Class-balanced resampling (with replacement)

Generates k training index sets. Total train size per set is floor(N \*
split), where N is the total number of samples. Each stratum contributes
(approximately) equally, sampling WITH replacement within each stratum.

## Usage

``` r
.balanced_bootstrap_mc(
  k,
  split,
  stratified,
  remainder = c("distribute", "drop")
)
```

## Arguments

- k:

  Integer. Number of resampling sets.

- split:

  Numeric in (0,1). Fraction of total N used as training size.

- stratified:

  List of length 3:

  1.  type: "DA" for classification or "R" for regression

  2.  Y: matrix with one column (n x 1)

  3.  probs: numeric vector of quantile probs for regression binning
      (e.g. c(0, .33, .66, 1))

- remainder:

  Character. How to handle n_train %% G remainder:

  - "drop": keep perfect balance, total size becomes G \*
    floor(n_train/G)

  - "distribute": distribute leftover 1-by-1 to random strata (still
    near-balanced)

## Value

List of length k. Each element is an integer vector of training indices.
