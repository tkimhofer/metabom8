# Check Validity of Selected Principal or Orthogonal Component

Internal consistency check for principal (or orthogonal) component
selection used in plotting or extraction functions.

## Usage

``` r
.check_pc_model(pc, mod, le = 1, type = "p")
```

## Arguments

- pc:

  Character or numeric. Selected component, e.g., `1` or `"o1"`.

- mod:

  An object of class `PCA_metabom8` or `OPLS_metabom8`.

- le:

  Integer. Expected maximum length of `pc`; usually 1.

- type:

  Character. Type of component to check: `"p"` for projection or `"t"`
  for scores. Currently only used for logic branching.

## Value

Nothing. Function is used for its side effect (i.e., error checking).
