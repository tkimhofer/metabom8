##### provenance hiit_raw.rda data

library(metabom8)

exp_type <- list(exp="PROF_PLASMA_CPMG")

path = '../HIIT7/'
exp_type = list(pulprog = c("noesygppr1d"), ns= list(op = "!=", value = 32))

hiit_raw <- read1d_proc(path, exp_type)
# dim(hiit_raw$X) # 8 131072

usethis::use_data(hiit_raw,  overwrite = TRUE, compress = "xz")

# For information on sample collection & processing and nmr data acquisition,
# please contact the package author tkimhofer@gmail.com
