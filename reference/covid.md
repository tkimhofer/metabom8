# COVID-19 blood plasma proton NMR spectra (processed)

1D proton NMR spectra from SARS-CoV-2–positive patients (n = 10) and
healthy controls (n = 13), collected in Perth, Western Australia.
Spectra were pre-processed (residual water and signal-free regions
excised, baseline corrected, and normalized to account for line-width
differences).

## Format

A numeric matrix/data frame with 23 rows (samples) and 27,819 columns
(chemical shift variables in parts per million, ppm).

## Source

Australian National Phenome Centre (ANPC), Murdoch University.

## Details

FIDs were acquired with a standard 90\\^{\circ}\\ RF pulse sequence on a
600 MHz Bruker Avance II spectrometer using IVDr methods for blood
plasma (300 K, 32 scans). The spectrometer was equipped with a
double-resonance broadband (BBI) probe and a refrigerated autosampler.

## References

Kimhofer et al. (2020)
[doi:10.1021/acs.jproteome.0c00280](https://doi.org/10.1021/acs.jproteome.0c00280)

## Examples

``` r
data(covid)
```
