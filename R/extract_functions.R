#' Extract number of peaks
#'
#' Function to apply peak picking parameters on a data subset
#'
#'
#' @param data An [MSnExp-class] object with peaks picked
#'
#' @export
extract_number_peaks <- function(data,
                                 save.table = FALSE,
                                 save.figure = FALSE){

  chrom_peaks_df <- as.data.frame(chromPeaks(data_centroid)) %>%
    dplyr::count(sample) %>%
    dplyr::rename('sample_index' = 'sample', 'total_peaks_detected' = 'n') %>%
    dplyr::left_join(rowid_to_column(Biobase::pData(data_centroid)),
                     by = c('sample_index' = 'rowid')) %>%
    dplyr::select(-FileName)

  par(mar=c(5,6,1,1))
  plotChromPeakImage(data_centroid,
                     binSize=100,
                     xlab="Retention time (sec)",
                     cex.sub = 0.5,
                     yaxt = "n",
                     main = 'Numbers of peaks detected')
  axis(2,
       at = seq(0,1, length.out = length(pData(data_centroid)$SampleID)),
       labels = pData(data_centroid)$SampleID,
       cex.axis = 1,
       las = 2)

  if(save.table == FALSE){
    print('Table not saved')
  } else {
    write_csv(chrom_peaks_df, save.table)
  }

  if(save.figure == FALSE){
    print('Figure not saved')
  } else {
    png(save.figure, width = 3000, height = 2000, res = 300)
    colfunc <- colorRampPalette(c("lightblue", "blue4"))
    par(mar=c(5,6,1,1))
    plotChromPeakImage(data_centroid,
                       binSize=100,
                       xlab="Retention time (sec)",
                       cex.sub = 0.5,
                       yaxt = "n",
                       main = 'Numbers of peaks detected')
    axis(2,
         at = seq(0,1, length.out = length(pData(data_centroid)$SampleID)),
         labels = pData(data_centroid)$SampleID,
         cex.axis = 1,
         las = 2)
    dev.off()
  }

  return(chrom_peaks_df)
}

