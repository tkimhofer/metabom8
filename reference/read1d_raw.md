# Read raw FIDs and process to spectra

Reads Bruker 1D NMR FIDs, corrects the digital filter (group delay),
applies apodisation (windowing), optional zero-filling, FFT, phasing and
ppm calibration.

Returns either absorption-, dispersion-, or magnitude-mode spectra.

If `to_global = TRUE`, objects are assigned to the global environment.

## Usage

``` r
read1d_raw(
  path,
  exp_type = list(exp = "PROF_PLASMA_CPMG128_3mm", pulprog = "noesygppr1d"),
  apodisation = list(fun = "exponential", lb = 0.2),
  zerofill = 1L,
  mode = c("absorption", "dispersion", "magnitude"),
  verbose = 1,
  recursive = TRUE,
  n_max = 1000,
  filter = TRUE,
  to_global = FALSE
)
```

## Arguments

- path:

  Character. Path to the root directory containing Bruker experiment
  folders.

- exp_type:

  Named list. Acquisition-parameter filter to select experiments (e.g.,
  `list(PULPROG = "noesygppr1d")`).

- apodisation:

  Named list. Apodisation function and parameters. `fun` must be one of
  `"exponential"`, `"cosine"`, `"sine"`, `"sem"`.

- zerofill:

  Integer. Zero-filling exponent (`1` doubles points, `2` quadruples,
  …).

- mode:

  Character. Spectrum type to return: `"absorption"`, `"dispersion"`, or
  `"magnitude"`.

- verbose:

  Integer. Verbosity level: `0` = silent, `1` = summary (default), `2` =
  detailed, `3` = debug.

- recursive:

  Logical. Recursively search subdirectories for FIDs.

- n_max:

  Integer. Maximum number of experiments to process.

- filter:

  Logical. Remove experiments with incomplete file structures.

- to_global:

  Logical. If TRUE, objects are also assigned to the global environment;
  otherwise, only an invisible list is returned.

## Value

A named list with three elements:

- X:

  A numeric matrix of spectra (rows = samples, columns = ppm values).

- ppm:

  A numeric vector of chemical-shift axis (ppm).

- meta:

  A data frame of acquisition metadata, row-aligned with `X`.

If `to_global = TRUE`, these objects are also assigned to the global
environment. In that case, any existing objects with the same names will
be overwritten.

## Details

FIDs are read from Bruker acquisition folders and processed by the
following pipeline:

1.  Digital-filter (group-delay) correction: the initial *n* complex
    points are invalid due to the causal DSP decimation filter and are
    discarded; *n* equals `GRPDLY` when present, or is looked up from
    Bruker tables indexed by `DECIM` and `DSPFVS` on older systems.

2.  Apodisation (windowing).

3.  Zero-filling (optional).

4.  FFT to the frequency domain.

5.  Phase correction.

6.  PPM calibration (e.g., to TSP).

A common ppm scale is then interpolated across spectra.

**Digital filter note:** On newer systems `GRPDLY` is written in
`acqus`/`acqu2s` and should be used directly. For older data sets
(`GRPDLY < 0` or missing), the group delay is derived from `DECIM` and
`DSPFVS` via an internal look-up table.

**Apodisation functions**:

- `"uniform"`

- `"cosine"`,

- `"sine"`,

- `"exponential"` (parameter: `lb`)

- `"sem"` (sine \* exponential, parameter: `lb`)

- `"gauss"` (parameter: `lb`, `'gb'`, and `'para'`)

- `"expGaus_resyG"` (parameter: `lb`, `'gb'`, and `'aq_t'`)

## See also

[`read1d_proc`](https://tkimhofer.github.io/metabom8/reference/read1d.md)
for importing TopSpin-processed spectra

## Examples

``` r
path <- system.file("extdata", package = "metabom8")
read1d_raw(
  path,
  exp_type    = list(PULPROG = "noesygppr1d"),
  apodisation = list(fun = "exponential", lb = 0.2),
  zerofill    = 1,
  n_max       = 3
)
#> Processing 2 experiments.
```
