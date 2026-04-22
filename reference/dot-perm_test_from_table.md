# Permutation-test summary from an opls_perm out_df table

Takes the `out_df` you showed (perm rows + one "non-permuted" row) and
returns observed values, permutation distributions, and permutation
p-values.

## Usage

``` r
.perm_test_from_table(
  out_df,
  observed_label = "non-permuted",
  alternative = c("greater", "less"),
  add_one = TRUE,
  na_rm = TRUE
)
```

## Arguments

- out_df:

  data.frame with columns like q2_comp, r2_comp, aucs_te, aucs_tr, model
  (with one row "non-permuted"), and optionally r/r_abs.

- observed_label:

  character. Label used in `model` for the observed row.

- alternative:

  character. "greater" (default) tests obs \> perm; "less" tests obs \<
  perm.

- add_one:

  logical. If TRUE, uses (1 + count)/(B + 1) p-value correction.

- na_rm:

  logical. Drop NA values before computing p-values.

## Value

A list with:

- `observed`: named numeric vector of observed metrics

- `perm`: list of numeric vectors for each metric

- `p_value`: named numeric vector of permutation p-values

- `B`: number of permutations used
