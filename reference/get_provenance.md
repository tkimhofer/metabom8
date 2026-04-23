# Retrieve metabom8 provenance metadata

Extracts preprocessing provenance stored in the `"m8_prep"` attribute.
Allows access to the full processing log, a specific step, or a specific
parameter within a step.

## Usage

``` r
get_provenance(x, step = NULL, param = NULL)
```

## Arguments

- x:

  A metabom8 object (named list with element `X`) or a numeric matrix
  containing metabom8 provenance metadata.

- step:

  Optional. Either:

  - Numeric index of the preprocessing step

  - Character string matching the recorded step name

  If `NULL`, the full provenance log is returned.

- param:

  Optional character string specifying a parameter name within the
  selected step. If provided, only this parameter value is returned.

## Value

Depending on the arguments:

- Full provenance list (if `step = NULL`)

- A single preprocessing step (if `step` specified)

- A single parameter value (if `param` specified)

## Details

Provenance metadata are recorded automatically by metabom8 preprocessing
functions and stored as structured attributes on the spectral matrix.
This function provides programmatic access to these records.

## See also

Other provenance:
[`add_note()`](https://tkimhofer.github.io/metabom8/reference/add_note.md),
[`print_provenance()`](https://tkimhofer.github.io/metabom8/reference/print_provenance.md)

## Examples

``` r
data(hiit_raw)

hiit_proc <- hiit_raw |>
  calibrate(type = "tsp") |>
  excise()

# Retrieve full log
log <- get_provenance(hiit_proc)

# Retrieve specific step
get_provenance(hiit_proc, step = 2)
#> $step
#> [1] "excise1d"
#> 
#> $params
#> $params$regions
#> $params$regions$upfield_noise
#> [1] -5.217943  0.250000
#> 
#> $params$regions$water
#> [1] 4.5 5.2
#> 
#> $params$regions$urea
#> [1] 5.5 6.0
#> 
#> $params$regions$downfield_noise
#> [1]  9.70000 14.80367
#> 
#> 
#> 
#> $notes
#> [1] "Specified chemical shift regions removed using get_idx()."
#> 
#> $time
#> [1] "2026-04-23 10:05:41.679568"
#> 
#> $pkg
#> [1] "metabom8"
#> 
#> $pkg_version
#> [1] "0.99.10"
#> 

# Retrieve parameter from a named step
get_provenance(hiit_proc, step = "calibrate", param = "target")
#> $mode
#> [1] "singlet"
#> 
#> $window
#> [1] -0.2  0.2
#> 
#> $centre
#> [1] 0
#> 
```
