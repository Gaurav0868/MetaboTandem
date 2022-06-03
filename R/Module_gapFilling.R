#' Gap fill UI
#'
#' @param id character used to specify namespace, see [`shiny::NS`][shiny::NS()]
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
    plotOutput(ns('no_gap_plot')),
    tableOutput(ns('num_filled')),
    verbatimTextOutput(ns('is_gap_filled'))
  )
}

gapFillingServer <- function(id, data_grouped){
  moduleServer(id, function(input, output, session){

    features <- featureValues(data_grouped())


    output$no_gap_plot <- renderPlot(

    )

  })

}
