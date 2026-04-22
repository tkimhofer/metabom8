# Check for Intact Bruker NMR File Structures

Verifies that each NMR experiment has a matching `acqus` and `fid` file.

## Usage

``` r
.check1d_files_fid(datapath, n_max = 10, filter = TRUE, recursive, verbose)
```

## Arguments

- datapath:

  Character. Directory containing NMR spectra.

- n_max:

  Integer. Maximum number of experiments to process.

- filter:

  Logical. Whether to filter out incomplete experiments.

- recursive:

  Logical. Whether to search directories recursively.

- verbose:

  Integer. Controls verbosity level (0 = silent, 1 = basic, 2 =
  detailed).

## Value

A list with paths to intact experiments and matching metadata.

## See also

[`.extract_acq_pars1d`](https://tkimhofer.github.io/metabom8/reference/dot-extract_acq_pars1d.md),
[`.filterExp_files`](https://tkimhofer.github.io/metabom8/reference/dot-filterExp_files.md)
