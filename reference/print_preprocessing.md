# List available preprocessing functions Returns the preprocessing utilities provided by metabom8.

List available preprocessing functions Returns the preprocessing
utilities provided by metabom8.

## Usage

``` r
print_preprocessing()
```

## Value

A named character vector describing preprocessing functions.

## See also

Other preprocessing:
[`align_segment()`](https://tkimhofer.github.io/metabom8/reference/align_segment.md),
[`align_spectra()`](https://tkimhofer.github.io/metabom8/reference/align_spectra.md),
[`binning()`](https://tkimhofer.github.io/metabom8/reference/binning.md),
[`calibrate()`](https://tkimhofer.github.io/metabom8/reference/calibrate.md),
[`correct_baseline()`](https://tkimhofer.github.io/metabom8/reference/correct_baseline.md),
[`correct_lw()`](https://tkimhofer.github.io/metabom8/reference/correct_lw.md),
[`pqn()`](https://tkimhofer.github.io/metabom8/reference/pqn.md)

## Examples

``` r
list_preprocessing()
#>                                               calibrate 
#>        "Chemical shift calibration to reference signal" 
#>                                                  excise 
#>  "Remove spectral regions (e.g., residual water, urea)" 
#>                                        baseline_correct 
#>                                   "Baseline correction" 
#>                                            norm_erectic 
#>              "Normalise spectra based on ERETIC signal" 
#>                                                     pqn 
#>      "Apply Probabilistic Quantile Normalisation (PQN)" 
#>                                                 binning 
#>                                           "Bin spectra" 
#>                                              correct_lw 
#>            "Apply line width correction (experimental)" 
#>                                           align_segment 
#> "Spectral alignment of a single chemical shift segment" 
#>                                           align_spectra 
#>        "Automated segment-wise full spectrum alignment" 
```
