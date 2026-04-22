# Evaluate model fit progression based on cross-validated performance

Determines whether an additional component should be fitted based on
cross-validated generalization performance. For regression models
(`type = "R"`), the decision is based on the \\Q^2\\ statistic. For
discriminant analysis models (`type = "DA"`), the decision is based on
the cross-validated AUROC (area under the ROC curve).

## Usage

``` r
.evalFit(
  type,
  q2s,
  cv_auc,
  pc_max = 5,
  min_delta = 0.05,
  sat_q2 = 0.98,
  sat_auc = 0.97,
  min_q2 = 0.05,
  min_auc = 0.6
)
```

## Arguments

- type:

  Character. Model type: `"R"` for regression or `"DA"` for discriminant
  analysis. The value may optionally include a suffix such as `"-mY"` to
  indicate multi-response models.

- q2s:

  Numeric vector of \\Q^2\\ values for the fitted components (used for
  regression models).

- cv_auc:

  Numeric vector of cross-validated AUROC values for the fitted
  components (used for classification models).

- pc_max:

  Integer. Maximum number of components allowed.

- min_delta:

  Numeric. Minimum improvement in the cross-validated performance metric
  required to justify fitting an additional component.

- sat_q2:

  Numeric. Saturation threshold for \\Q^2\\ in regression models.

- sat_auc:

  Numeric. Saturation threshold for cross-validated AUROC in
  discriminant analysis models.

- min_q2:

  Numeric. Minimum acceptable \\Q^2\\ value for the first component.

- min_auc:

  Numeric. Minimum acceptable cross-validated AUROC for the first
  component.

## Value

A list containing:

- stop:

  Logical indicating whether model fitting should stop.

- reason:

  Character string describing the stopping condition (e.g., `"pc_max"`,
  `"cv_decreased"`, `"cv_improvement_negligible"`).

- metric:

  The current cross-validated performance value.

- delta:

  Improvement relative to the previous component.

## Details

AUROC is computed on cross-validation test folds to ensure that it
reflects out-of-sample classification performance, analogous to the use
of \\Q^2\\ in regression.
