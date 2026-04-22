# OPLS Model Validation via Y-Permutation

Performs Y-permutation tests to assess the robustness of OPLS models by
comparing model performance statistics on real vs. permuted response
labels.

## Usage

``` r
opls_perm(smod, n = 10, plot = TRUE, mc = FALSE)
```

## Arguments

- smod:

  An OPLS model object of class `OPLS_metabom8`.

- n:

  Integer. Number of permutations to perform.

- plot:

  Logical. If `TRUE`, generates a visual summary of permutation
  statistics.

- mc:

  Logical. If `TRUE`, enables multicore processing (currently not
  implemented).

## Value

A `data.frame` with model metrics (e.g., R2, Q2, AUROC) from both
permuted and original models.

## Details

Each permutation shuffles the response labels and fits a new OPLS model.
The function captures model statistics (R2, Q2, AUC) to compare against
the non-permuted model. This helps determine whether the original model
performance is better than expected by chance.

## References

Wiklund, S. et al. (2008). A Tool for Improved Validation of OPLS
Models. *Journal of Chemometrics*, 22(11–12), 594–600.

## See also

Other model_validation:
[`cv_anova()`](https://tkimhofer.github.io/metabom8/reference/cv_anova.md),
[`dmodx()`](https://tkimhofer.github.io/metabom8/reference/dmodx.md)
