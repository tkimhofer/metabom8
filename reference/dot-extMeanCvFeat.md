# Extract Cross-Validated OPLS Features

Extracts a specific data structure (e.g., scores or predictions) from
cross-validation output returned by
[`.oplsComponentCv()`](https://tkimhofer.github.io/metabom8/reference/dot-oplsComponentCv.md).
Depending on the cross-validation type and model structure, the output
is aggregated across folds by computing the mean, standard deviation
(SD), and coverage.

Here, **Coverage** refers to the proportion of cross-validation folds in
which a given sample's value was observed (i.e., not missing), serving
as an estimate of how frequently that value was evaluated across folds.

## Usage

``` r
.extMeanCvFeat(cv_obj, feat = "t_xp")
```

## Arguments

- cv_obj:

  Named list. Output of
  [`.oplsComponentCv()`](https://tkimhofer.github.io/metabom8/reference/dot-oplsComponentCv.md)
  for a specific OPLS component.

- feat:

  Character. Feature name to extract from each CV fold object (e.g.,
  `"t_xp"` or `"y_pred_test"`).

## Value

A numeric matrix or array of the requested feature:

- For Monte Carlo CV and single-column Y: matrix with rows `"mean"`,
  `"sd"`, and `"coverage"` x samples.

- For Monte Carlo CV and multi-column Y: 3D array with shape
  `[3, samples, outcomes]`.

- For k-fold CV: matrix or array containing values for each sample (no
  aggregation, just fold values).
