#' Gap fill UI
#'
#' @param id character used to specify namespace, see [`shiny::NS`][shiny::NS()]
#'
#' @return a \code{shiny::\link[shiny]{tagList}} containing UI elements

gapFillingUI <- function(id){
  ns <- NS(id)
  tagList(
    h2('Gap filling'),
    fluidRow(
      column(6,
             actionButton(ns('gap_fill'), 'Apply Alignment and Correspondence',
                          class = 'btn-lg btn-success')
      ),
      column(6,
             tableOutput(ns('no_gap_table'))
      )
    ),
    verbatimTextOutput(ns('is_gap_filled'))
  )
}

gapFillingServer <- function(id, data_grouped){
  moduleServer(id, function(input, output, session){

    output$no_gap_table <- renderTable(
      featureValues(data_grouped()),
      striped = TRUE,
      bordered = TRUE
    )

  })

}
