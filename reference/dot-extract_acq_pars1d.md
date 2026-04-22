# Extract Bruker NMR Acquisition Parameters

Parses `acqus` files to extract acquisition parameters for each
experiment.

## Usage

``` r
.extract_acq_pars1d(f_list)
```

## Arguments

- f_list:

  List. Contains paths to NMR experiment folders (e.g., as returned by
  [`.check1d_files_fid()`](https://tkimhofer.github.io/metabom8/reference/dot-check1d_files_fid.md)).

## Value

A data frame containing extracted acquisition parameters for each
experiment.

## See also

[`.filterExp_files`](https://tkimhofer.github.io/metabom8/reference/dot-filterExp_files.md),
[`.check1d_files_fid`](https://tkimhofer.github.io/metabom8/reference/dot-check1d_files_fid.md)
