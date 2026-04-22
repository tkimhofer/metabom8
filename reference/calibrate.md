# Chemical Shift Calibration

Aligns 1D \\^1\\H NMR spectra to a reference signal on the existing ppm
grid. Supports singlet (e.g. TSP or custom range) and predefined doublet
references (glucose, alanine).

## Usage

``` r
calibrate(X, ppm = NULL, type = "tsp")
```

## Arguments

- X:

  Numeric matrix/vector of spectra or a metabom8 data list.

- ppm:

  Numeric chemical shift vector. If `X` is a metabom8 list, this is
  taken from `X$ppm`.

- type:

  Character (`"tsp"`, `"glucose"`, `"alanine"`), or list.

## Value

Calibrated spectra in the same structure as input.

## Details

In addition to predefined references, custom calibration targets can be
supplied as a list with elements:

- `mode`: `"singlet"` or `"doublet"`

- `window`: numeric vector of length 2 defining the ppm search region

- `centre`: target ppm position (optional; defaults to mean(window))

- `j`: numeric vector of length 2 specifying expected J-coupling
  (required for doublet calibration)

Custom doublet calibration requires the expected J-coupling range (`j`)
in ppm to distinguish the two peaks of the multiplet.

For example:

    calibrate(
      X, ppm,
      list(
        mode   = "doublet",
        window = c(1.2, 1.35),
        j      = c(0.007, 0.009)
      )
    )

## See also

Other preprocessing:
[`align_segment()`](https://tkimhofer.github.io/metabom8/reference/align_segment.md),
[`align_spectra()`](https://tkimhofer.github.io/metabom8/reference/align_spectra.md),
[`binning()`](https://tkimhofer.github.io/metabom8/reference/binning.md),
[`correct_baseline()`](https://tkimhofer.github.io/metabom8/reference/correct_baseline.md),
[`correct_lw()`](https://tkimhofer.github.io/metabom8/reference/correct_lw.md),
[`pqn()`](https://tkimhofer.github.io/metabom8/reference/pqn.md),
[`print_preprocessing()`](https://tkimhofer.github.io/metabom8/reference/print_preprocessing.md)

## Examples

``` r
data("covid_raw")
X=covid_raw$X
ppm=covid_raw$ppm
X_tsp <- calibrate(X, ppm, type = "tsp")
X_glu <- calibrate(X, ppm, type = "glucose")
X_custom <- calibrate(X, ppm, type = c(1.9, 2.1))
```
