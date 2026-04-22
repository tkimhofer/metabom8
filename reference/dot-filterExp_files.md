# Filter Bruker NMR Experiments

Filters Bruker 1D experiments based on acquisition metadata and
optionally limits the number of returned entries. Character parameters
are matched by value, while numeric parameters support equality, set
inclusion, range filtering (list(range = c(min, max))), and generic
operators (list(op = "\>", value = x)).

## Usage

``` r
.filterExp_files(pars, exp_type, f_list, n_max)
```

## Arguments

- pars:

  Data frame. Parsed acquisition and processing parameters.

- exp_type:

  Named list. Filtering conditions for parameters. Elements may be
  character vectors, numeric values/vectors, or operator/range lists for
  numeric fields.

- f_list:

  List. File paths of spectra and associated metadata (e.g., f_fid,
  f_1r).

- n_max:

  Integer. Maximum number of experiments to retain.

## Value

A list containing filtered and ordered `f_list` and `pars`.

## See also

[`.extract_acq_pars1d`](https://tkimhofer.github.io/metabom8/reference/dot-extract_acq_pars1d.md),
[`.check1d_files_fid`](https://tkimhofer.github.io/metabom8/reference/dot-check1d_files_fid.md)
