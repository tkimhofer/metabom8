# Check for Intact File Systems - Helper Function for read1d

Scans the specified directory recursively to find intact sets of Bruker
NMR files: 'procs', 'acqus', and '1r' files. Optionally filters
incomplete experiments.

## Usage

``` r
.detect1d_procs(
  datapath,
  n_max = 10,
  filter = TRUE,
  recursive = TRUE,
  verbose = 1
)
```

## Arguments

- datapath:

  character. Path to directory containing spectra.

- n_max:

  integer. Maximum number of spectra to read (default 10).

- filter:

  logical. Whether to filter out incomplete file sets (default TRUE).

- recursive:

  logical. Whether to search directories recursively.

- verbose:

  integer. Verbosity level for messaging.

## Value

A list with elements:

- path:

  Character vector of experiment folder paths (intact).

- exp_no:

  Character vector of experiment folder IDs.

- f_procs:

  Character vector of paths to 'procs' files.

- f_acqus:

  Character vector of paths to 'acqus' files.

- f_1r:

  Character vector of paths to '1r' files.
