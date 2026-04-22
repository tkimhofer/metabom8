# m8_model class Model object returned by `pca()`, `pls()`, and `opls()`.

m8_model class Model object returned by
[`pca()`](https://tkimhofer.github.io/metabom8/reference/pca.md),
[`pls()`](https://tkimhofer.github.io/metabom8/reference/pls.md), and
[`opls()`](https://tkimhofer.github.io/metabom8/reference/opls.md).

## Usage

``` r
# S4 method for class 'm8_model'
summary(object)

# S4 method for class 'm8_model'
show(object)
```

## Arguments

- object:

  An object of class `m8_model`.

## Value

An object of class `m8_model`.

## Functions

- `summary(m8_model)`: Summarise model performance and component
  selection.

- `show(m8_model)`: Show a compact model header.

## Slots

- `engine`:

  Character. Model engine ("pca", "pls", "opls").

- `ctrl`:

  List. Engine-specific control and performance information.

- `fit`:

  List. Fitted data (engine specific).

- `cv`:

  Resampling instance (may be NULL if not used).

- `prep`:

  Scaling and centering information

- `provenance`:

  Preprocessing attributes of spectral matrix X

- `session`:

  R-session information

- `call`:

  Function call

- `dims`:

  List with `n` and `p`.

## Examples

``` r
data(covid)
cv <- balanced_mc(k=5, split=2/3)
scaling <- uv_scaling(center=TRUE)
model <-opls(X=covid$X, Y=covid$an$type, scaling, cv)
#> Performing discriminant analysis.
#> An O-PLS-DA model with 1 predictive and 1 orthogonal components was fitted.
class(model)
#> [1] "m8_model"
#> attr(,"package")
#> [1] "metabom8"
show(model)
#> 
#> m8_model <opls>
#> ----------------------------------------
#> Dimensions  : 10 samples x 27819 variables
#> Mode        : classification 
#> Preprocess  : center | UV 
#> Components  : 2 (3 tested)
#> Validation  : BalancedMonteCarlo (k = 5) 
#> Stop rule   : cv_improvement_negligible 
#> ----------------------------------------
#> Use summary() for performance metrics.
#> 
```
