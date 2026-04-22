# COVID-19 blood plasma proton NMR spectra (raw)

1D proton NMR spectra from SARS-CoV-2–positive patients (n = 10) and
healthy controls (n = 13), collected in a research study in Perth,
Western Australia. Spectra are raw and require processing before
statistical analysis (see
[`?covid`](https://tkimhofer.github.io/metabom8/reference/covid.md) for
processed spectra).

## Format

A numeric matrix/data frame with 23 rows and 29,782 columns:

- rows:

  Spectra (samples)

- columns:

  Chemical shift variables in parts per million (ppm)

## Source

Australian National Phenome Centre (ANPC), Murdoch University.

## Details

FIDs were acquired using a Carr–Purcell–Meiboom–Gill (CPMG) pulse
sequence on a 600 MHz Bruker Avance II spectrometer using IVDr methods
for blood plasma (300 K, 32 scans). The spectrometer was equipped with a
double-resonance broadband (BBI) probe and a refrigerated autosampler at
4\\^{\circ}\\C.

## References

Kimhofer et al. (2020)
[doi:10.1021/acs.jproteome.0c00280](https://doi.org/10.1021/acs.jproteome.0c00280)

## Examples

``` r
data(covid_raw)
```
