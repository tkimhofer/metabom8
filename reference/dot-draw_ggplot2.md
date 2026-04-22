# Draw NMR Spectra Using ggplot2

Internal backend renderer using `ggplot2`. Converts spectra to long
format prior to plotting.

## Usage

``` r
.draw_ggplot2(dat, ...)
```

## Arguments

- dat:

  List returned by
  [`.prep_spec()`](https://tkimhofer.github.io/metabom8/reference/dot-prep_spec.md).

- ...:

  Additional arguments for `ggplot2`.
