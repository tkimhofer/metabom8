#' COVID-19 blood plasma proton NMR spectra (processed).
#'
#' 1D NMR spectra from SARS-CoV 2 positive patients (n=10) and healthy controls (n=13), collected within a research study framework performed in Perth, Western Australia. Spectra were pre-processed which included excision of residual water and signal-free areas, baseline correction and normalisatoin to account for different line widths.
#' FID's were acquired using a standard 90˚ RF pulse sequence on a 600 MHz Avance II spectrometer using Bruker IVDr methods for blood plasms (300 K, 32 scans). The spectrometer was equipped with a double resonance broadban probe (BBI) and a refrigerated autosampler that was operated at 4˚C. Further information on the study setup and NMR pulse sequencuence and FID processing can be found in Kimhofer et al., 'Integrative Modelling of Quantitative Plasma Lipoprotein, Metabolic and Amino Acid Data Reveals a Multi-organ Pathological Signature of SARS-CoV-2 Infection'. Journal of Proteome Research (Aug 2020).
#'
#' @format A data frame with 11 rows and 131072 colums:
#' \describe{
#' \item{rows}{represent spectra}
#' \item{colums}{represent chemical shift variables in part per million (ppm)}
#' }
#' #' ppm num array, chemical shift position, its lengths matches the number of columns in X
#' meta, data frame, containing spectrometer metadata.
#'
#' Spectral data have not been pre-processed for statistical analysis.
#' @source Australian National Phenome Centre (ANPC) @ Murdoch University \url[https://www.murdoch.edu.au/research/institutes-centres/health-futures-institute/australian-national-phenome-centre/covid-19-critical-research-programme]{ANPC}
#' @examples
#' data(covid)
#'covid'




#' COVID-19 blood plasma proton NMR spectra (raw)
#'
#' 1D proton NMR spectra from SARS-CoV 2 positive patients (n=10) and healthy controls (n=13), collected within a research study framework performed in Perth, Western Australia. Spectra are raw and require processig before statistical analysis (see \code{?data(covid) for processed spectra).}).
#' FID's were acquired using a Carr-Purcell-Meiboom-Gill pulse sequence (spin echo) pulse sequence on a 600 MHz Avance II spectrometer using Bruker in vitro diagnostics research (IVDr) methods for blood plasma (300 K, 32 scans). The spectrometer was equipped with a double resonance broadban probe (BBI) and a refrigerated autosampler that was operated at 4˚C. Detailed information on the study setup, incl. NMR pulse sequencuence and FID processing can be found in Kimhofer et al., 'Integrative Modelling of Quantitative Plasma Lipoprotein, Metabolic and Amino Acid Data Reveals a Multi-organ Pathological Signature of SARS-CoV-2 Infection'. Journal of Proteome Research (Aug 2020).
#'
#' @format A data frame with 11 rows and 131072 colums:
#' \describe{
#' \item{rows}{represent spectra}
#' \item{colums}{represent chemical shift variables in part per million (ppm)}
#' }
#'
#' Spectral data have not been pre-processed for statistical analysis.
#'
#' @source Australian National Phenome Centre (ANPC) @ Murdoch University \url[https://www.murdoch.edu.au/research/institutes-centres/health-futures-institute/australian-national-phenome-centre/covid-19-critical-research-programme]{ANPC}
#' @examples
#' data(covid_raw)
#'covid_raw'
