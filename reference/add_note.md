# Add user note to metabom8 provenance Appends a user annotation to the `"m8_prep"` attribute. The step title is formatted as `"note {username}"`. The timestamp is stored in `params`, and the user message is stored in `notes`.

Add user note to metabom8 provenance Appends a user annotation to the
`"m8_prep"` attribute. The step title is formatted as
`"note {username}"`. The timestamp is stored in `params`, and the user
message is stored in `notes`.

## Usage

``` r
add_note(x, note, params = NULL)
```

## Arguments

- x:

  A metabom8 object or matrix with `"m8_prep"` metadata.

- note:

  Character string describing the annotation.

- params:

  Named list providing paramter key-value pairs.

## Value

The input object with updated provenance metadata.

## See also

Other provenance:
[`get_provenance()`](https://tkimhofer.github.io/metabom8/reference/get_provenance.md),
[`print_provenance()`](https://tkimhofer.github.io/metabom8/reference/print_provenance.md)

## Examples

``` r
params <- list(
  runtime  = "docker",
  image    = Sys.getenv("IMAGE", "docker-image-dummy"),
  workflow = Sys.getenv("M8_WORKFLOW", "std_prof-urine"),
  agent    = paste0("snakemake/", Sys.getenv("SNAKEMAKE_VERSION", "v?")),
  run_id   = Sys.getenv("M8_RUN_ID", "m8-2605-001")
)

data(hiit_raw)
print_provenance(hiit_raw)
#> No metabom8 preprocessing metadata found.

hiit_proc <- hiit_raw |>
  calibrate(type = "tsp") |>
  excise() |>
  add_note('dilution-adaptive acquisition mode -> verify snr after normalisation',
    params)

print_provenance(hiit_proc)
#> metabom8 processing pipeline:
#> ==============================
#> [1] calibrate
#>      params:
#>         target :
#>             mode   : singlet
#>             window : [-0.2, 0.2]
#>             centre : 0
#>      notes: Chemical shift calibration applied.
#> 
#> [2] excise1d
#>      params:
#>         regions :
#>             upfield_noise   : [-5.21794, 0.25]
#>             water           : [4.5, 5.2]
#>             urea            : [5.5, 6]
#>             downfield_noise : [9.7, 14.8037]
#>      notes: Specified chemical shift regions removed using get_idx().
#> 
#> [3] user note
#>      params:
#>         runtime  : docker
#>         image    : docker-image-dummy
#>         workflow : std_prof-urine
#>         agent    : snakemake/v?
#>         run_id   : m8-2605-001
#>      notes: dilution-adaptive acquisition mode -> verify snr after normalisation
#> 
```
