# metabom8 1.2.0 (2025-08-11)

## Added
- GitHub Actions: R CMD check matrix and pkgdown site deployment.

## Changed
- Minimum R version is now 4.1.
- Switched docs/labels to plotmath/UTF-8 escapes for Greek letters (e.g., `alpha`) to avoid locale issues.git 
- Standardized plotting text to be Unicode-safe across OSes.

## Fixed
- Non-portable filenames in pkgdown outputs (`docs/` excluded via `.Rbuildignore`).
- Vignette install error on macOS (ensure `inst/doc` not ignored).
- Rd cross-reference warning caused by `[0, 1]` being parsed as links.

## Internal
- Build uses C++17 via `src/Makevars` (`CXX_STD = CXX17`); removed non-portable compiler flags.
- Added unit tests (`testthat`) and improved example coverage.


# metabom8 1.1.0 (prior to 2025)

## Added
- Datasets: `covid` (processed) and `covid_raw` (raw) with documentation.
- Expanded preprocessing and modeling workflows (PCA/OPLS) 
- Apdated examples for ANPC release.




