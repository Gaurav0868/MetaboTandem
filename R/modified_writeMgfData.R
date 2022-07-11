##
## MODIFIED the original writeMgfData
## from https://github.com/lgatto/MSnbase/blob/master/R/readWriteMgfData.R
##
## by Nathalia Graf Grachet
##
## Now, this function write the correct polarity of the features
## (original automatically added positive), and
## the title of the features are now the actual feature identification instead of
## the long description

# Based on the code contributed by Guangchuang Yu <guangchuangyu@gmail.com>
# Modified by Sebastian Gibb <mail@sebastiangibb.de>
# setMethod("mod_writeMgfData",
#           signature = signature("Spectrum"),
#           function(object,
#                    con = "spectrum.mgf",
#                    COM = NULL,
#                    TITLE = NULL) {
#               mod_writeMgfDataFile(list(object), con = con, COM = COM, TITLE = TITLE,
#                                verbose = isMSnbaseVerbose())
#           })
#
# setMethod("mod_writeMgfData",
#           signature = signature("MSnExp"),
#           function(object,
#                    con = "experiment.mgf",
#                    COM = NULL,
#                    verbose = isMSnbaseVerbose()) {
#               mod_writeMgfDataFile(spectra(object), con = con, COM = COM,
#                                verbose = verbose)
#           })

#' @param addFields `data.frame` or `matrix` with optional additional
#'     fields to be added to each spectrum (each column one field).
#'
#' @noRd
mod_writeMgfDataFile <- function(splist, con, COM = NULL, TITLE = NULL,
                                 verbose = MSnbase::isMSnbaseVerbose(),
                                 addFields = NULL) {

  if (class(con) == "character" && file.exists(con)) {
    message("Overwriting ", con, "!")
    unlink(con)
  }

  if (length(addFields)) {
    if (length(dim(addFields)) != 2)
      stop("'addFields' has to be a matrix or data.frame.")
    if (!is.matrix(addFields))
      addFields <- do.call(cbind, lapply(addFields, as.character))
    if (is.null(colnames(addFields)))
      stop("Column names required on 'addFields'.")
    if (nrow(addFields) != length(splist))
      stop("nrow of 'addFields' has to match length of 'splist'")
  }

  if (class(con)[1] == "character") {
    con <- file(description = con, open = "at")
    on.exit(close(con))
  }

  if (is.null(COM)) {
    COM <- paste0(ifelse(length(splist) <= 1, "Spectrum", "Experiment"),
                  "exported by MSnbase on ", date())
  }
  cat(paste0("COM=",COM), file = con, sep = "")

  verbose <- verbose & length(splist) > 1

  if (verbose)
    pb <- utils::txtProgressBar(min = 0, max = length(splist), style = 3)

  for (i in seq(along=splist)) {
    if (verbose)
      utils::setTxtProgressBar(pb, i)

    mod_writeMgfContent(splist[[i]],
                        TITLE = splist@elementMetadata@listData[["feature_id"]][[i]],
                        con = con,
                        addFields = addFields[i, ])
  }
  if (verbose)
    close(pb)
}

mod_writeMgfContent <- function(sp, TITLE = NULL, con, addFields = NULL) {
  .cat <- function(..., file = con, sep = "", append = TRUE) {
    cat(..., file = file, sep = sep, append = append)
  }

  .cat("\nBEGIN IONS\n",
       "SCANS=", MSnbase::acquisitionNum(sp))

  if (is.null(TITLE)) {
    .cat("\nTITLE=msLevel ", MSnbase::msLevel(sp),
         "; retentionTime ", MSnbase::rtime(sp),
         "; scanNum ", MSnbase::acquisitionNum(sp))

    if (length(MSnbase::scanIndex(sp))) {
      .cat("; scanIndex ", MSnbase::scanIndex(sp))
    }

    if (MSnbase::msLevel(sp) > 1) {
      .cat("; precMz ", MSnbase::precursorMz(sp),
           "; precCharge ", MSnbase::precursorCharge(sp))
    }
  } else {
    .cat("\nTITLE=", TITLE)
  }

  signal = ""

  if(MSnbase::polarity(sp) == 1){
    signal = "+"
  } else signal = "-"

  if (MSnbase::msLevel(sp) > 1) {
    .cat("\nRTINSECONDS=", MSnbase::rtime(sp),
         "\nPEPMASS=", MSnbase::precursorMz(sp))
    if (length(MSnbase::precursorCharge(sp)) &&
        !is.na(MSnbase::precursorCharge(sp))) {
      .cat("\nCHARGE=", MSnbase::precursorCharge(sp), signal)
    }
  } else .cat("\nRTINSECONDS=", MSnbase::rtime(sp))

  if (length(addFields) && !is.null(names(addFields)))
    .cat("\n", paste(toupper(names(addFields)),
                     addFields, sep = "=", collapse = "\n"))

  .cat("\n", paste(MSnbase::mz(sp), MSnbase::intensity(sp), collapse = "\n"))
  .cat("\nEND IONS\n")
}
