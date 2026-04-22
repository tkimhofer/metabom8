# OPLS Component Estimation via Cross-Validation

Performs orthogonal projections to latent structures (OPLS) modeling on
training folds of a cross-validation (CV) scheme. The function fits one
orthogonal and one predictive component per fold, tracking scores,
predictions, and residuals. For the first orthogonal component
(`nc = 1`), the full model structure is initialized. For later
components (`nc > 1`), modeling proceeds on residual matrices carried
over from previous iterations.

## Usage

``` r
.oplsComponentCv(X, Ycs_fold, cv.set, nc, mod.cv, acc)
```

## Arguments

- X:

  Numeric matrix. Input data matrix (samples x features); only required
  for `nc = 1`.

- Ycs_fold:

  List containing the response values (`Y`) for each cross-validation
  split.

- cv.set:

  List of integer vectors. Each element contains training indices for
  one CV round. These indices are adjusted internally by -1 due to
  0-based indexing in the underlying C++ (Rcpp) routines.

- nc:

  Integer. Component number currently being estimated (`nc = 1` for
  first orthogonal component).

- mod.cv:

  List. Model state to be updated with results from each CV iteration.

## Value

A list of the same length as `cv.set`. Each list element contains:

- `t_xo`: Orthogonal component scores (matrix, samples x components).

- `t_xp`: Predictive component scores (matrix, samples x 1).

- `y_pred_train`: Predicted responses for training samples.

- `y_pred_test`: Predicted responses for held-out test samples.

- `x_res`: Residual matrix after filtering out orthogonal structure.

## Note

Training indices in `cv.set` are internally adjusted by -1 before being
passed to Rcpp routines. This is required because C++ uses zero-based
indexing.
