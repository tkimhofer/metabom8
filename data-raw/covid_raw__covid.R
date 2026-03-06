##### provenance covid_raw.rda / covid.rda

library(metabom8)
exp_type <- list(exp="PROF_PLASMA_CPMG")


##### RAW covid dataset
path <- "../rack_1"
covid_raw  <- read1d_proc(path, exp_type, n_max=2)

save(covid_raw, file='data/covid_raw.rda')


##### PROCESSED / ANNOTATED covid dataset
path <- "../f2"
covid_raw_f2 <- read1d_proc(path, exp_type)

covid_f2 <-
  covid_raw_f2 |>
  calibrate(type = "glucose") |>
  excise(
    regions = list(
      uf=c(min(ppm), 0.25),
      rw=c(4,55, 4.8),
      df=c(9.0, max(ppm)))
  ) |>
  correct_baseline(method='asls')

an <- read.csv('../tk/annotation_covid.csv')
covid <- list(X=covid_f2$X, ppm=covid_f2$ppm, an=an)

save(covid, file='data/covid.rda')


# For information on sample collection & processing and nmr data acquisition,
# please refer to the original publication doi:10.1021/acs.jproteome.1c00273
