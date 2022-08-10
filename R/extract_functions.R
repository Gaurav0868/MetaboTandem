#' Extract number of peaks
#'
#' Function to apply peak picking parameters on a data subset
#'
#'
#' @param data An [MSnExp-class] object with peaks picked
#' @param save_table Path to save the extracted table or `NULL`
#' @param save_figure Path to save the figure or `NULL`
#'
#' @export
extract_number_peaks <- function(data,
                                 save_table = NULL,
                                 save_figure = NULL){

  chrom_peaks_df <- as.data.frame(xcms::chromPeaks(data)) %>%
    dplyr::count(sample) %>%
    dplyr::rename('sample_index' = 'sample', 'total_peaks_detected' = 'n') %>%
    dplyr::left_join(tibble::rowid_to_column(Biobase::pData(data)),
                     by = c('sample_index' = 'rowid')) %>%
    dplyr::select(-FileName)

  par(mar=c(5,6,1,1))
  xcms::plotChromPeakImage(data,
                     binSize=100,
                     xlab="Retention time (sec)",
                     cex.sub = 0.5,
                     yaxt = "n",
                     main = 'Numbers of peaks detected')
  axis(2,
       at = seq(0,1, length.out = length(Biobase::pData(data)$SampleID)),
       labels = Biobase::pData(data)$SampleID,
       cex.axis = 1,
       las = 2)

  if(is.null(save_table)){
    print('Table not saved')
  } else {
    readr::write_csv(chrom_peaks_df, save_table)
  }

  if(is.null(save_figure)){
    print('Figure not saved')
  } else {
    png(save_figure, width = 3000, height = 2000, res = 300)
    colfunc <- grDevices::colorRampPalette(c("lightblue", "blue4"))
    par(mar=c(5,6,1,1))
    xcms::plotChromPeakImage(data,
                       binSize=100,
                       xlab="Retention time (sec)",
                       cex.sub = 0.5,
                       yaxt = "n",
                       main = 'Numbers of peaks detected')
    axis(2,
         at = seq(0,1, length.out = length(xcms::pData(data)$SampleID)),
         labels = xcms::pData(data)$SampleID,
         cex.axis = 1,
         las = 2)
    dev.off()
  }

  return(chrom_peaks_df)
}


#' Extract features
#'
#' Function to apply peak picking parameters on a data subset
#'
#' @param data An [MSnExp-class] object with peaks picked.
#' @param metadata Sample information data.frame.
#' @param save_table Path to save the extracted table or `NULL`
#'
#' @export
extract_features <- function(data){

  ## extract feature values after filling in
  feature_ab_matrix <- as.data.frame(xcms::featureValues(data,
                                                          value="into",
                                                          method="maxint")) %>%
    tibble::rownames_to_column(var = 'FeatureID')
  ## replace NA with zero
  feature_ab_matrix[is.na(feature_ab_matrix)] <- 0
  ## replace file name with sample name
  colnames(feature_ab_matrix)[-1] <- paste0(MSnbase::pData(data)$SampleID)

  return(feature_ab_matrix)

}


#' Extract feature definitions
#'
#' Function to apply peak picking parameters on a data subset
#'
#' @param data An [MSnExp-class] object with peaks picked.
#' @param feature_abundance Data.frame with feature abundances per sample.
#' @param save_table Path to save the extracted table or `NULL`
#'
#' @export
extract_feature_definition <- function(data,
                                       feature_abundance,
                                       save_table = NULL){
  ## Get feature definitions and intensities
  feature_definition <- as.data.frame(xcms::featureDefinitions(data)) %>%
    tibble::rownames_to_column(var = 'FeatureID') %>%
    dplyr::select(-peakidx) %>%
    dplyr::left_join(feature_abundance, by = 'FeatureID')

  if(is.null(save_table)){
    print('Table not saved')
  } else {
    readr::write_csv(feature_definition, save_table)
  }

  return(feature_definition)
}


#' Extract spectra table
#'
#' Function to apply peak picking parameters on a data subset
#'
#' @param data An [MSnExp-class] object with peaks picked.
#' @param save_table Path to save the extracted table or `NULL`
#'
#' @export
extract_spectra_table <- function(data, save_table = NULL){

  spectra_table <- MSnbase::fData(data) %>%
    tibble::rownames_to_column(var = 'SpectraID')

  if(is.null(save_table)){
    print('Table not saved')
  } else {
    readr::write_csv(spectra_table, save_table)
  }

  return(spectra_table)
}


#' Extract MS2 data
#'
#' Function to apply peak picking parameters on a data subset
#'
#' @param data An [MSnExp-class] object with peaks picked.
#' @param save_table Path to save the extracted table or `NULL`
#'
#' @export
extract_MS2_consensus <- function(data){

  # This function uses custom functions for database compatibility from the
  # Github Repo:
  # https://github.com/jorainer/xcms-gnps-tools

  print('Starting extraction')
  BiocParallel::register(BiocParallel::SerialParam())

  filtMs2Spectra <- xcms::featureSpectra(data, return.type = "MSpectra",
                                         msLevel = 2)
  filtMs2Spectra <- MSnbase::clean(filtMs2Spectra,
                                   all = TRUE)
  filtMs2Spectra <- formatSpectraForGNPS(filtMs2Spectra)
  Ms2_consensus <- MSnbase::combineSpectra(filtMs2Spectra,
                                           fcol = "feature_id",
                                           method = MSnbase::consensusSpectrum,
                                           mzd = 0,
                                           minProp = 0.5,
                                           ppm = 25,
                                           intensityFun = median,
                                           mzFun = median)

  print('Extraction finished')

  return(Ms2_consensus)
}

