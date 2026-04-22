# Distance to the Model in X-Space (DModX)

Calculates the orthogonal distance of each observation to an OPLS model
in X-space. DModX can be used for identifying outliers.

## Usage

``` r
dmodx(mod, plot = TRUE)
```

## Arguments

- mod:

  An OPLS model object of class `m8_model` (engine = "opls").

- plot:

  Logical. If TRUE, a plot of DModX values with an approximate cutoff is
  shown.

## Value

A data frame with columns `ID`, `DmodX`, and `passed`.

## Details

DModX is computed from the X-residual matrix as a scaled RMSE of
residuals. The empirical cutoff (dashed line in plot) uses a t-based
confidence interval and assumes approximate normality of DModX values.
This assumption may not be satisfied in all datasets, so the resulting
threshold should be regarded as a pragmatic heuristic for outlier
detection.

## See also

Other model_validation:
[`cv_anova()`](https://tkimhofer.github.io/metabom8/reference/cv_anova.md),
[`opls_perm()`](https://tkimhofer.github.io/metabom8/reference/opls_perm.md)

## Examples

``` r
data(covid)
cv <- balanced_mc(k=5, split=2/3)
scaling <- uv_scaling(center=TRUE)
model <-opls(X=covid$X, Y=covid$an$type, scaling, cv)
#> Performing discriminant analysis.
#> An O-PLS-DA model with 1 predictive and 1 orthogonal components was fitted.
dX <- dmodx(model)

print(dX[[1]])
#>  [1]  1  2  3  4  5  6  7  8  9 10
df <- dX[[2]]; head(df)
#> [1] 1.0456047 0.9667600 0.3198258 0.8996049 0.9106332 0.9668600
```
