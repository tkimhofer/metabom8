# Print metabom8 preprocessing pipeline

Displays the recorded preprocessing history attached to a metabom8
spectral matrix. The function reads the `"m8_prep"` attribute and prints
each processing step in chronological order, including parameters and
notes.

## Usage

``` r
print_provenance(x, detail = FALSE, max_items = 8)
```

## Arguments

- x:

  A metabom8 object or spectral matrix with attached `"m8_prep"`
  provenance metadata.

- detail:

  Prints full information trail, incl. timestamp / versioning

- max_items:

  Limits individual list entries (e.g., parameter) to specified number

## Value

Invisibly returns `NULL`. This function is called for its side effect of
printing pipeline information.

## Details

The input can be either:

- A metabom8-style list containing `$X`, or

- A matrix with metabom8 provenance attributes attached.

If no preprocessing metadata is found, a message is printed.

metabom8 records preprocessing steps as an ordered list of
transformation descriptors stored in the `"m8_prep"` attribute. Each
step typically contains:

- `step`: Name of the preprocessing operation

- `params`: Parameters used

- `notes`: Optional description

- `time`: Timestamp

- `pkg`: Package name and version

This function provides a compact audit trail of the processing workflow,
facilitating reproducibility and provenance inspection.

## See also

Other provenance:
[`add_note()`](https://tkimhofer.github.io/metabom8/reference/add_note.md),
[`get_provenance()`](https://tkimhofer.github.io/metabom8/reference/get_provenance.md)

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
