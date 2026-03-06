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
#'
#' @examples
#' data(covid)
#' @docType data
#' @keywords datasets
#' @name covid
NULL

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
#'
#' @examples
#' data(covid_raw)
#' @docType data
#' @keywords datasets
#' @name covid_raw
NULL

#' High-intensity interval training (HIIT) 1H NMR urine dataset
#'
#' Urine samples collected from a single individual performing a \eqn{VO_2max}-type
#' exercise protocol over a time period of 3h.
#'
#' @docType data
#'
#' @format A list with four elements:
#' \describe{
#'   \item{X}{Numeric matrix of spectral intensities (samples x variables).}
#'   \item{ppm}{Numeric vector of chemical shift values corresponding to columns of \code{X}.}
#'   \item{meta}{Data frame containing acquisition metadata.}
#'   \item{df}{Data frame containing sample annotations.}
#' }
#'
#' @details
#' The spectra were acquired on a 600 MHz Bruker NMR spectrometer.
#' The dataset is included for demonstration of preprocessing and
#' modelling workflows in \pkg{metabom8}.
#'
#' @source Example dataset bundled with the package.
#' @examples
#' data(hiit_raw)
#'
#' @keywords datasets
#' @name hiit_raw
NULL
