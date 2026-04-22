# Compute X residual matrix Returns the residual matrix (E) of an OPLS model.

Compute X residual matrix Returns the residual matrix (E) of an OPLS
model.

## Usage

``` r
xres(object)

# S4 method for class 'm8_model'
xres(object)
```

## Arguments

- object:

  An object.

## Value

Numeric matrix containing the X residuals.

## Examples

``` r
data(covid)
cv <- balanced_mc(k = 5, split = 2/3)
scaling <- uv_scaling(center = TRUE)
model <- opls(X = covid$X, Y = covid$an$type, scaling, cv)
#> Performing discriminant analysis.
#> An O-PLS-DA model with 1 predictive and 1 orthogonal components was fitted.
X_res <- xres(model)
dim(X_res) == dim(covid$X)
#> [1] TRUE TRUE
```
