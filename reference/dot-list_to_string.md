# Convert named list to compact string representation

Internal utility that converts a named list into a single character
string of the form `name=value`. Nested sublists are represented
recursively using braces.

## Usage

``` r
.list_to_string(lst)
```

## Arguments

- lst:

  A named list.

## Value

A character string.

## Details

This function is primarily used for compact parameter logging in
provenance metadata and should not be called directly by users.
