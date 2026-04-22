# k-fold cross-validation index generator

Creates a list of training set indices for k-fold cross-validation (CV),
without considering class balance in Y. Each element of the list
represents the row indices used for training in one CV fold.

## Usage

``` r
.kFold(k, Y)
```

## Arguments

- k:

  Integer. Number of CV folds. If `k` is not valid or too high relative
  to number of rows in Y, it defaults to leave-one-out CV (LOO-CV).

- Y:

  Matrix. Outcome matrix (observations x variables).

## Value

A list of length `k`, where each element contains the training indices
(integers) for one CV fold.
