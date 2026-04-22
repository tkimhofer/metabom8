# Class-balanced Monte Carlo Cross-Validation (MCCV) splits

Generates class/group-balanced training set indices for each MCCV round.

## Usage

``` r
.mcBalanced(k, split, stratified)
```

## Arguments

- k:

  Integer. Number of MCCV training sets to generate.

- split:

  Numeric. Fraction of observations per class to include in each
  training set. Must be between 0 and 1.

- stratified:

  List of three elements:

  - 1: Character. Outcome type: 'R' (regression) or 'DA' (discriminant
    analysis), optionally with suffix '-mY' for multi-column Y.

  - 2: Matrix. Response matrix `Y`, with `ncol(Y) == 1`.

  - 3: Numeric vector of probabilities used to stratify numeric Y if
    regression.

## Value

A list of length `k`, each element containing a vector of row indices
for the training set.

## Details

Each round samples a class-balanced subset without replacement. Useful
for high-variance, imbalanced-class modeling.
