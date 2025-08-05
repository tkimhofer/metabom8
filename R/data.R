#' Bariatric surgery metabolic profiling data
#'
#' Metabolic profiling NMR dataset of murine urine samples collected before and after bariatric surgery (Roux-en-Y gastric bypass).
#'
#' One dimensional proton NMR spectra of murine urine samples were acquired on a Bruker Avance III spectrometer (600.13 MHz, 300 K). A standard NMR pulse sequence (\emph{recycle delay - 90° - t1 - 90° - tm - 90° - acquisition}) was applied, where \emph{t1} was 3 μs and \emph{tm} (mixing time) was 100 ms. Water suppression was achieved using selective irradiation during a recycle delay of 2 s and \emph{tm}. A 90° pulse of 10 μs was used. A total of 128 scans were collected into 64k data points with a spectral width of 20 ppm.
#'
#' For further information on sample collection and processing see Li et al. (2011).
#'
#' @format A list with three elements:
#' \describe{
#'   \item{X.pqn}{Matrix of pre-processed proton NMR spectra (rows = spectra)}
#'   \item{ppm}{Chemical shift vector (ppm) with length equal to number of columns in \code{X.pqn}}
#'   \item{meta}{Spectrometer metadata}
#' }
#'
#' @references Li, Jia V., et al. (2011) Metabolic surgery profoundly influences gut microbial-host metabolic cross-talk. \emph{Gut}. 60(9), 1214-1223. \url{https://www.ncbi.nlm.nih.gov/pubmed/21572120}
#'
#' @usage data(bariatric)
#' @docType data
#' @keywords datasets
"bariatric"


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
