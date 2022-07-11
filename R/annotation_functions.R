#' Adduct annotation with CAMERA
#'
#' Function to apply peak picking parameters on a data subset
#'
#'
#' @param data An [MSnExp-class] object with peaks picked
#' @param group_by column
#'
#' @export
#'
#'
annotate_CAMERA <- function(data, group_by){
  ms1_data <- xcms::filterMsLevel(data, msLevel = 1L)
  xs <- as(ms1_data, 'xcmsSet')

  # Setting up sampclass
  sampclass(xs) <- MSnbase::pData(data)[[group_by]]
  sampnames(xs) <- MSnbase::pData(data)$SampleID

  an <- CAMERA::xsAnnotate(xs)
  anF <- CAMERA::groupFWHM(an, perfwhm = 0.6)
  anI <- CAMERA::findIsotopes(anF, mzabs = 0.01)
  anIC <- CAMERA::groupCorr(anI, cor_eic_th = 0.75)
  anFA <- CAMERA::findAdducts(anIC, polarity="positive")


  peaklist <- CAMERA::getPeaklist(anFA)

  return(peaklist)
}

#' Annotation with GNPS
#'
#'
