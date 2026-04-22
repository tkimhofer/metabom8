# Import 1D NMR spectra (TopSpin processed)

Imports TopSpin-processed 1D NMR spectra together with spectrometer
acquisition and TopSpin processing parameters (`acqus` and `procs`,
respectively).

## Usage

``` r
read1d(
  path,
  exp_type = list(pulprog = "noesygppr1d"),
  n_max = 1000,
  filter = TRUE,
  recursive = TRUE,
  verbose = 1,
  to_global = FALSE
)

read1d_proc(
  path,
  exp_type = list(pulprog = "noesygppr1d"),
  n_max = 1000,
  filter = TRUE,
  recursive = TRUE,
  verbose = 1,
  to_global = FALSE
)
```

## Arguments

- path:

  Character. Directory path containing Bruker NMR experiments.

- exp_type:

  Named list. Optional filtering specification based on acquisition or
  processing metadata. Each list element must correspond to a metadata
  field (e.g. `pulprog`, `ns`, `rg`). Filtering supports:

  - Exact match: `list(pulprog = "noesygppr1d")`

  - Membership: `list(pulprog = c("zg30", "noesygppr1d"))`

  - Numeric range: `list(ns = list(range = c(16, 128)))`

  - Generic comparison: `list(ns = list(op = ">=", value = 32))`

  Multiple fields are combined using logical AND.

- n_max:

  Integer. Maximum number of spectra to import. Default: 1000.

- filter:

  Logical. Filter out experiments with incomplete file systems.

- recursive:

  Logical. Search `path` recursively. Default: TRUE.

- verbose:

  Logical or numeric. Verbosity level.

- to_global:

  Logical. If `TRUE`, the returned objects are additionally assigned to
  the global environment.

## Value

A named list with three elements:

- X:

  A numeric matrix of spectra (rows = samples, columns = ppm values).

- ppm:

  A numeric vector of chemical shift values (ppm).

- meta:

  A data frame of acquisition and processing metadata, row-aligned with
  `X`.

If `to_global = TRUE`, objects with the same names in the global
environment will be overwritten.

## Examples

``` r
path <- system.file("extdata", package = "metabom8")

read1d_proc(path, exp_type = list(pulprog = "noesygppr1d"), n_max = 2)
#> Imported 2 spectra.
```
