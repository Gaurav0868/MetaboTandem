#' Gap fill UI
#'
#' @param id character used to specify namespace, see [`shiny::NS`][shiny::NS()]
#'
#' @importFrom xcms featureValues
#'
#' @return a \code{shiny::\link[shiny]{tagList}} containing UI elements

gapFillingUI <- function(id){
  ns <- NS(id)
  tagList(
    h2('Before Gap filling'),
    plotOutput(ns('gap_plot')),
    actionButton(ns('apply_gap'), 'Apply Gap Filling',
                 class = 'btn-lg btn-success'),
    h2('After Gap Filling'),
    plotOutput(ns('gap_filled_plot')),
    tableOutput(ns('num_filled')),
    verbatimTextOutput(ns('is_gap_filled'))
  )
}

gapFillingServer <- function(id, data_grouped){
  moduleServer(id, function(input, output, session){

    features <- reactive({
      xcms::featureValues(data_grouped()) %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var = 'FeatureID') %>%
        tidyr::pivot_longer(!FeatureID,
                            names_to = 'sampleid',
                            values_to = 'abundance') %>%
        tidyr::drop_na(abundance) %>%
        dplyr::count(sampleid, name = 'Initial_peaks')
    }) %>%
      bindEvent(data_grouped())

    output$gap_plot <- renderPlot(

      features() %>%
        dplyr::mutate(peak = 'Initial_peaks') %>%
        ggplot() +
        geom_col(aes(x = sampleid,
                     y = Initial_peaks,
                     fill = peak),
                 color = 'black') +
        labs(y = 'Number of features') +
        scale_fill_manual(values = c('Initial_peaks' = 'blue3')) +
        theme_bw() +
        scale_y_continuous(expand = c(0,0,0.05,0),
                           breaks = scales::pretty_breaks()) +
        theme(axis.title.x = element_blank(),
              legend.position = 'bottom',
              legend.title = element_blank())
    )

    data_gap_filled <- reactive({
      notid <- showNotification('Gap filling...',
                                duration = NULL, closeButton = FALSE)
      on.exit(removeNotification(notid), add = TRUE)

      apply_gap_filling(data_grouped())
    }) %>%
      bindEvent(input$apply_gap)

    features_gap_filled <- reactive({
      xcms::featureValues(data_gap_filled()) %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var = 'FeatureID') %>%
        tidyr::pivot_longer(!FeatureID,
                            names_to = 'sampleid',
                            values_to = 'abundance') %>%
        tidyr::drop_na(abundance) %>%
        dplyr::count(sampleid) %>%
        dplyr::left_join(features(), by = 'sampleid') %>%
        dplyr::mutate(Gap_filled = n - Initial_peaks) %>%
        tidyr::pivot_longer(cols = c(Initial_peaks, Gap_filled),
                            names_to = 'peaks',
                            values_to = 'counts')
    }) %>%
      bindEvent(data_gap_filled())

    output$gap_filled_plot <- renderPlot(

      ggplot(features_gap_filled()) +
        geom_col(aes(x = sampleid,
                     y = counts,
                     fill = peaks),
                 color = 'black') +
        labs(y = 'Number of features') +
        scale_fill_manual(values = c('Initial_peaks' = 'blue3',
                                     'Gap_filled' = 'firebrick')) +
        theme_bw() +
        scale_y_continuous(expand = c(0,0,0,1),
                           breaks = scales::pretty_breaks()) +
        theme(axis.title.x = element_blank(),
              legend.position = 'bottom',
              legend.title = element_blank())
    )

    return(data_gap_filled)

  })
}
