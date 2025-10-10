# perpare data for package
library(metabom8)

load("CoV_rawCPMG.Rdata") # imports R objects: Xc, ppm, meta, an
### loading selected spectra from 1D CPMG experiments of human plasma samples
### data represent a selection of NMR spectra reported in this publication:
### https://doi.org/10.1021/acs.jproteome.1c00273

## following pre-processing was applied

matspec(Xc, ppm, shift = c(5.15, 5.3))
Xc[is.na(Xc)] <- 0

# calibrate to TSP (use examples to calibrate to glucose)
Xc <- calibrate(Xc, ppm, type = "tsp")
matspec(Xc, ppm, shift = c(-0.1, 0.1))

idx <- which(an$ind_timepoint == 1 | an$sample_lab == "Healthy")
X <- Xc[idx, ]
meta <- meta[idx, ]

rownames(X) <- 1:nrow(X)
rownames(meta) <- rownames(X)

# keep only 3 raw spectra (vignette uses nmrdata lib)
idx <- get.idx(c(-0.1, 10), ppm)
ppm <- ppm[idx]
X <- X[1:2, idx]
meta <- meta[1:2, ]

save(X, ppm, meta, file = "data/covid_raw.rda")


