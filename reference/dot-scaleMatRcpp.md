# Column-wise matrix scaling (C++ backend)

Column-wise matrix scaling (C++ backend)

## Usage

``` r
.scaleMatRcpp(X, idc, center, scale_type)
```

## Arguments

- X:

  num matrix

- idc:

  int row indices of X

- center:

  logical mean centering

- scale_type:

  int 0: no scaling, 1: SD scaling, 2: Pareto scaling

## Value

list: 1. scale X matrix, 2. mean (sd vec), 3: sd (num vec)
