## XCMS-GNPS Integration tools
## Taken from https://github.com/jorainer/xcms-gnps-tools

#' @title Format MS2 spectra for export in GNPS-MGF format
#'
#' @description
#'
#' Re-format MS2 spectrum information for export of the data in Mgf format
#' supported by GNPS. In detail, the function replaces the acquisition number
#' of each spectrum with the feature ID (expected to be present in the
#' `"feature_id"` column of `mcols(x)`) converted to an integer by removing
#' the ID's leading `"FT"`.
#'
#' @param x `Spectra`.
#'
#' @return `Spectra` with the acquisition number replaced.
#'
#' @author Johannes Rainer
#'
#' @noRd
formatSpectraForGNPS <- function(x) {
  fids <- S4Vectors::mcols(x)$feature_id
  if (!length(fids))
    stop("No column named 'feature_id' present in 'mcols(x)'")
  fids <- as.integer(sub("^FT", "", fids))
  S4Vectors::mendoapply(x, fids, FUN = function(z, id) {
    z@acquisitionNum <- id
    z
  })
}
