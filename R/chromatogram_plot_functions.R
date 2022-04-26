#' Calculate TIC
#'
#' Calculate and plot Total Ion Chromatogram
#'
#'
#' @param data An [MSnExp-class] object.
#' @param metadata Sample information file.
#' @param show_plot Wheter to show the plot.
#' @param save_figure Path to save the figure or `NULL`.
#' @param color_by Column to color output.
#'
#' @export
calculate_tic <- function(data,
                          metadata,
                          show_plot = TRUE,
                          save_figure = NULL,
                          color_by = 'AUTO'){

  tic <- xcms::chromatogram(data, aggregationFun="sum")

  color_vector <- create_color_vector(metadata, color_by = color_by)

  tic_plot <- plot(tic,
                   col = color_vector,
                   ylab = "Intensity",
                   xlab = "Retention Time (sec)",
                   main = 'TIC')

  if(show_plot){
    tic_plot
    legend("topright",
           legend = unique(names(color_vector)),
           col = unique(color_vector),
           lty=1)
  }

  if(!is.null(save_figure)){
    print('TIC plot will not be saved')
  } else{
    png(save_figure, width = 3000, height = 2000, res = 300)
    plot(tic,
         col = color_vector,
         ylab = "Intensity",
         xlab = "Retention Time (sec)",
         main = 'TIC')
    legend("topright",
           legend = unique(names(color_vector)),
           col = unique(color_vector),
           lty=1)
    dev.off()
  }

  res <- list(tic = tic,
              tic_plot = tic_plot)

  return(res)
}

#' Calculate BPC
#'
#' Calculate and plot Base Peak Chromatogram
#'
#'
#' @param data An [MSnExp-class] object.
#' @param metadata Sample information file.
#' @param show_plot Whether to show the plot.
#' @param save_figure Path to save the figure or `NULL`.
#' @param color_by Column to color output.
#'
#' @export
calculate_bpc <- function(data,
                          metadata,
                          show_plot = TRUE,
                          save_figure = NULL,
                          color_by = 'AUTO'){
  #### Calculate, plot and save total ion chromatogram

  bpc <- xcms::chromatogram(data, aggregationFun="max")

  color_vector <- create_color_vector(metadata, color_by = color_by)

  bpc_plot <- plot(bpc,
                   col = color_vector,
                   ylab = "Intensity",
                   xlab = "Retention Time (sec)",
                   main = 'BPC')

  if(show_plot){
    bpc_plot
    legend("topright",
           legend = unique(names(color_vector)),
           col = unique(color_vector),
           lty=1)
  }

  if(is.null(save_figure)){
    print('BPC plot will not be saved')
  } else{
    png(save_figure, width = 3000, height = 2000, res = 300)
    plot(bpc,
         col = color_vector,
         ylab = "Intensity",
         xlab = "Retention Time (sec)",
         main = 'BPC')

    legend("topright",
           legend = unique(names(color_vector)),
           col = unique(color_vector),
           lty=1)
    dev.off()
  }

  res <- list(bpc = bpc,
              bpc_plot = bpc_plot)

  return(res)
}
