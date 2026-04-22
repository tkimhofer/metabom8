# List available preprocessing steps

Returns a named character vector describing the preprocessing operations
implemented in metabom8. For more information on individual
functionalities please refer to the function help pages.

## Usage

``` r
list_preprocessing()
```

## Value

A named character vector with preprocessing function identifiers as
names and short descriptions as values.

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
