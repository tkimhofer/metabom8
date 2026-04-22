# \\R^2\\ and \\Q^2\\ Calculation for OPLS Models

Computes the coefficient of determination (\\R^2\\ or \\Q^2\\) for OPLS
regression models using the formula: \$\$R^2 = 1 -
\frac{PRESS}{TSS}\$\$, where PRESS is the prediction error sum of
squares, and TSS is the total sum of squares. If `ytss` is not provided,
it is computed directly from `Y`.

This function supports matrix inputs (e.g., for multi-class outcomes)
and averages across columns to return a single summary \\R^2\\/\\Q^2\\
value.

## Usage

``` r
.r2(Y, Yhat, ytss = NULL)
```

## Arguments

- Y:

  Numeric matrix. True response values. For classification, this should
  be a dummy matrix.

- Yhat:

  Numeric matrix. Predicted response values from the model.

- ytss:

  Numeric (optional). Total sum of squares of `Y`. If not provided, it
  will be computed internally.

## Value

A single numeric value representing \\R^2\\ or \\Q^2\\.
