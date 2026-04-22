# High-intensity interval training (HIIT) 1H NMR urine dataset

Urine samples collected from a single individual performing a
\\VO_2max\\-type exercise protocol over a time period of 3h.

## Format

A list with four elements:

- X:

  Numeric matrix of spectral intensities (samples x variables).

- ppm:

  Numeric vector of chemical shift values corresponding to columns of
  `X`.

- meta:

  Data frame containing acquisition metadata.

- df:

  Data frame containing sample annotations.

## Source

Example dataset bundled with the package.

## Details

The spectra were acquired on a 600 MHz Bruker NMR spectrometer. The
dataset is included for demonstration of preprocessing and modelling
workflows in metabom8.

## Examples

``` r
data(hiit_raw)
```
