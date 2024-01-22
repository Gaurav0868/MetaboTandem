#' Autotuner UI
#'`
#' @param id character used to specify namespace, see [`shiny::NS`][shiny::NS()]
#'
#' @return a \code{shiny::\link[shiny]{tagList}} containing UI elements

autotunerUI <- function(id){
  ns <- NS(id)
  dashboardPage(
    app_header(),
    dashboardSidebar(
      fluidRow(
        column(6, align = 'center', offset = 3,
               shinyWidgets::actionBttn(inputId = 'goHome_autotuner',
                                        icon = icon('house'),
                                        style = 'material-circle',
                                        color = 'success',
                                        size = 'sm'))
      )
    ),
    dashboardBody(
      tagList(
        h2('Parameters for Autotuner'),
        fluidRow(
          box(
            metadataUI(ns('metadata_autotuner'))
          ),
          column(6,
                 textInput(ns('group'), 'Group'),
                 numericInput(ns('lag'), 'Lag', value = 25, step = 1),
                 numericInput(ns('threshold'), 'Threshold', value = 3.1, step = 0.1),
                 numericInput(ns('influence'), 'Influence', value = 0.1, step = 0.1),
                 shinyWidgets::actionBttn(ns('get_signals'), 'Get Signals',
                                          style = 'material-flat',
                                          color = 'primary',
                                          size = 'md')

          )
        ),
        plotOutput(ns('signals_plot')),
        uiOutput(ns('finish_autotuner')),
        tableOutput(ns('final_parameters'))
      )
    )
  )


}

autotunerServer <- function(id){
  moduleServer(id, function(input, output, session){
    ns <- NS(id)
    metadata <- metadataServer('metadata_autotuner')

    autotuner_obj <- reactive({
      start_autotuner(metadata(),
                      group = input$group,
                      lag = input$lag,
                      threshold = input$threshold,
                      influence = input$influence,
                      plot = TRUE)
    }) %>%
      bindEvent(input$get_signals)

    output$signals_plot <- renderPlot(
      start_autotuner(metadata(),
                      group = input$group,
                      lag = input$lag,
                      threshold = input$threshold,
                      influence = input$influence,
                      plot = TRUE)
    ) %>%
      bindEvent(input$get_signals)

    output$finish_autotuner <- renderUI(
      tagList(
        shinyWidgets::actionBttn(ns('get_params'), 'Get Parameters',
                                 style = 'material-flat',
                                 color = 'success',
                                 size = 'md')
      )
    ) %>%
      bindEvent(autotuner_obj)

    output$final_parameters <- renderTable(
      extract_autotuner(autotuner_obj = autotuner_obj(),
                        massThr = 0.005),
      striped = TRUE,
      bordered = TRUE
    ) %>%
      bindEvent(input$get_params)


  })
}
