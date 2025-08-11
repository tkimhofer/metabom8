#' COVID-19 blood plasma proton NMR spectra (processed)
#'
#' 1D proton NMR spectra from SARS-CoV-2–positive patients (n = 10) and healthy controls (n = 13),
#' collected in Perth, Western Australia. Spectra were pre-processed (residual water and signal-free
#' regions excised, baseline corrected, and normalized to account for line-width differences).
#'
#' FIDs were acquired with a standard 90\eqn{^{\circ}} RF pulse sequence on a 600 MHz Bruker Avance II
#' spectrometer using IVDr methods for blood plasma (300 K, 32 scans). The spectrometer was equipped
#' with a double-resonance broadband (BBI) probe and a refrigerated autosampler.
#'
#' @format A numeric matrix/data frame with 23 rows (samples) and 27,819 columns
#'   (chemical shift variables in parts per million, ppm).
#' @source Australian National Phenome Centre (ANPC), Murdoch University.
#' @references Kimhofer et al. (2020) \doi{10.1021/acs.jproteome.0c00280}
#' @usage data(covid)
#' @examples
#' data(covid)
#' dim(covid)
#' @docType data
#' @keywords datasets
#' @name covid
"covid"

#' COVID-19 blood plasma proton NMR spectra (raw)
#'
#' 1D proton NMR spectra from SARS-CoV-2–positive patients (n = 10) and healthy controls (n = 13),
#' collected in a research study in Perth, Western Australia. Spectra are raw and require processing
#' before statistical analysis (see \code{?covid} for processed spectra).
#'
#' FIDs were acquired using a Carr–Purcell–Meiboom–Gill (CPMG) pulse sequence on a 600 MHz Bruker
#' Avance II spectrometer using IVDr methods for blood plasma (300 K, 32 scans). The spectrometer was
#' equipped with a double-resonance broadband (BBI) probe and a refrigerated autosampler at
#' 4\eqn{^{\circ}}C.
#'
#' @format A numeric matrix/data frame with 23 rows and 29,782 columns:
#' \describe{
#'   \item{rows}{Spectra (samples)}
#'   \item{columns}{Chemical shift variables in parts per million (ppm)}
#' }
#'
#' @source Australian National Phenome Centre (ANPC), Murdoch University.
#' @references Kimhofer et al. (2020) \doi{10.1021/acs.jproteome.0c00280}
#' @usage data(covid_raw)
#' @examples
#' data(covid_raw)
#' dim(covid_raw)
#' @docType data
#' @keywords datasets
#' @name covid_raw
"covid_raw"
