#' Annotation using SIRIUS UI
#'`
#' @param id character used to specify namespace, see [`shiny::NS`][shiny::NS()]
#'
#' @return a \code{shiny::\link[shiny]{tagList}} containing UI elements

siriusAnnotationUI <- function(id){
  ns <- NS(id)
  tagList(
    box(title = 'Important',
        solidHeader = TRUE,
        status = 'primary',
        width = 12,
        p('Since version 5, SIRIUS requires to have a ',
          strong('licence'), 'and ', strong('account'),
          'registered to use their ',
          'webservices. Please see more information at their ',
          'website: https://bio.informatik.uni-jena.de/software/sirius/'),
        p('If you do not have an account SIRIUS can still be used for ',
          'molecular formula annotation using fragmentation trees')),
    box(title = 'Setup',
        solidHeader = TRUE,
        status = 'primary',
        width = 12,
        shinyWidgets::materialSwitch(ns('account'),
                                     label = 'Do you have a SIRIUS account?',
                                     value = FALSE,
                                     status = 'success'),
        uiOutput(ns('account_info')),
        br(),
        actionButton(ns('annotate'), 'Annotate with SIRIUS')),

    dataTableOutput(ns('results'))

  )
}

#' Database annotation server
#'
#' @param id character used to specify namespace, see [`shiny::NS`][shiny::NS()]
#' @param ms1_data A [data.frame] with columns: 'name', 'mz' and 'rt'
#'

siriusAnnotationServer <- function(id, data_proc){
  moduleServer(id, function(input, output, session){
    ns <- NS(id)

    output$account_info <- renderUI({
      if(input$account){
        tagList(
          textInput(ns('user'), label = 'Email'),
          textInput(ns('pass'), label = 'Password')
        )
      }
    })

    annotation_df <- reactive({
      notid <- showNotification('Annotation with SIRIUS',
                                duration = NULL, closeButton = FALSE)
      on.exit(removeNotification(notid), add = TRUE)
      ms2.data <- suppressWarnings(extract_MS2_consensus(data_proc()))
      annotate_SIRIUS(ms2.data)
    }) %>%
      bindEvent(input$annotate)

    output$results <- renderDataTable({
      annotation_df()
    })

  })
}
