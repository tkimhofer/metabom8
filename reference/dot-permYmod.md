# OPLS Y-permutation Modeling

Performs OPLS modeling on permuted Y to estimate cross-validated model
performance under the null hypothesis. Extracts predictive performance
(R2, Q2, AUROC) from cross-validation for either regression or
classification.

## Usage

``` r
.permYmod(XcsTot, Y, cv, type, nc_o)
```

## Arguments

- XcsTot:

  Numeric matrix. Input data matrix (samples x features).

- Y:

  Numeric matrix. Response variable (numeric or dummy-coded).

- cv:

  List. Cross-validation parameters, including `method` and `cv_sets`.

- type:

  Character. Model type: `"DA"`, `"R"`, possibly with `"-mY"` suffix for
  multi-column Y.

- nc_o:

  Integer. Number of orthogonal components to be fitted.

## Value

List with elements:

- `r2_comp`: Numeric, R2 statistic for test data.

- `q2_comp`: Numeric, Q2 statistic from cross-validation.

- `aucs_tr`: Numeric, AUROC on training data (for classification only).

- `aucs_te`: Numeric, AUROC on test data (for classification only).

## Details

When `nc_o > 1`, the function fits additional orthogonal components
sequentially, updating the model object. Classification performance is
computed using AUROC (via `pROC`). Regression performance uses R2 and
Q2.
