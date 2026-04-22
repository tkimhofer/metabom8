# Read Bruker NMR Parameter Files

Helper function used by
[`read1d()`](https://tkimhofer.github.io/metabom8/reference/read1d.md)
to extract acquisition (`acqus`) and processing (`procs`) parameters
from Bruker-formatted 1D NMR experiments.

## Usage

``` r
.extract_pars1d(f_list)
```

## Arguments

- f_list:

  List of file paths. Output from
  [`.detect1d_procs()`](https://tkimhofer.github.io/metabom8/reference/dot-detect1d_procs.md).

## Value

A data frame of extracted acquisition and processing metadata. Row names
correspond to spectrum filenames.
