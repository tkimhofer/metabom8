#' COVID-19 blood plasma proton NMR spectra (processed)
#'
#' 1D NMR spectra from SARS-CoV-2 positive patients (n=10) and healthy controls (n=13), collected in Perth, Western Australia.
#' Spectra were pre-processed including excision of residual water and signal-free areas, baseline correction, and normalization to account for different line widths.
#' FIDs were acquired using a standard 90˚ RF pulse sequence on a 600 MHz Avance II spectrometer using Bruker IVDr methods for blood plasma (300 K, 32 scans).
#' The spectrometer was equipped with a double resonance broadband probe (BBI) and a refrigerated autosampler o


#' COVID-19 blood plasma proton NMR spectra (raw)
#'
#' 1D proton NMR spectra from SARS-CoV-2 positive patients (n=10) and healthy controls (n=13), collected in a research study in Perth, Western Australia. Spectra are raw and require processing before statistical analysis (see \code{?covid} for processed spectra).
#' FIDs were acquired using a Carr-Purcell-Meiboom-Gill (CPMG) pulse sequence on a 600 MHz Avance II spectrometer using Bruker IVDr methods for blood plasma (300 K, 32 scans). The spectrometer was equipped with a double resonance broadband probe (BBI) and a refrigerated autosampler at 4˚C.
#' Detailed information is available in Kimhofer et al., *Integrative Modelling of Quantitative Plasma Lipoprotein, Metabolic and Amino Acid Data Reveals a Multi-organ Pathological Signature of SARS-CoV-2 Infection*, Journal of Proteome Research (Aug 2020).
#'
#' @format A data frame with 11 rows and 131,072 columns:
#' \describe{
#'   \item{rows}{Spectra (samples)}
#'   \item{columns}{Chemical shift variables in parts per million (ppm)}
#' }
#'
#' @source Australian National Phenome Centre (ANPC), Murdoch University
#'
#' @examples
#' data(covid_raw)
#' head(covid_raw)
#'
#' @docType data
#' @keywords datasets
"covid_raw"
