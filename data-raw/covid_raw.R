# perpare data for package
usethis::use_data("covid_raw")

load("/Users/TKimhofer/Rproj/cov19/CoV_rawCPMG_calibrated.Rdata")
matspec(Xc, ppm, shift = c(5.15, 5.3))
Xc[is.na(Xc)] <- 0

# calibrate to TSP (use examples to calibrate to glucose)
Xc <- calibrate(Xc, ppm, type = "tsp")
matspec(Xc, ppm, shift = c(-0.1, 0.1))

idx <- which(an$ind_timepoint == 1 | an$sample_lab == "Healthy")
# an=an[idx,colnames(an) %in% c('type', 'hospital')]
X <- Xc[idx, ]
meta <- meta[idx, ]

# rownames(an)=1:nrow(an)
rownames(X) <- 1:nrow(X)
rownames(meta) <- rownames(X)

# generate raw spectra, these need only be two, stats will be performed with
# processed data
idx <- get.idx(c(-0.1, 10), ppm)
ppm <- ppm[idx]
X <- X[1:2, idx]
meta <- meta[1:2, ]

save(X, ppm, meta, file = "data/covid_raw.rda")


