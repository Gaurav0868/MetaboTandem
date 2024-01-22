#' Align spectra UI
#'
#' @param id character used to specify namespace, see [`shiny::NS`][shiny::NS()]
#'
#' @return a \code{shiny::\link[shiny]{tagList}} containing UI elements

alignSpectraUI <- function(id){
  ns <- NS(id)
  tagList(
    h2('Spectra alignment parameters'),
    p('Select parameters for aligment and apply using the ',
      strong('Align Spectra'),
      'button.'),
    fluidRow(
      column(6,
             selectizeInput(ns('group'),
                            label = 'Choose grouping variable for alignment',
                            choices = NULL),
             br(),
             numericInput(ns('min_frac'),
                          'Minimum fraction of samples in a group a peak must,
                          be present',
                          value = 0.5,
                          max = 1),
             br(),
             numericInput(ns('min_samples'),
                          'Minimum number of samples in a group a peak must,
                          be present',
                          value = 1),
             br(),
             actionButton(ns('align'), 'Apply Alignment and Correspondence',
                          class = 'btn-lg btn-success')
      ),
      column(6,
             plotOutput(ns('aling_plot'))
      )
    ),
    verbatimTextOutput(ns('is_aligned')),
    br(),
    verbatimTextOutput(ns('is_grouped')),
    br(),
    plotOutput(ns('align_plot')),
    uiOutput(ns('next_buttonSA')),
    actionButton(inputId = 'back_buttonSA',
                 label = 'Back',
                 icon = icon('arrow-left'))
  )
}


alignSpectraServer <- function(id, metadata, data_pp){
  moduleServer(id, function(input, output, session){

    observe({
      updateVarSelectizeInput(session, 'group', data = metadata(),
                              server = TRUE)
    })

    data_aligned <- reactive({
      notid <- showNotification('Aligning spectra...',
                                duration = NULL, closeButton = FALSE)
      on.exit(removeNotification(notid), add = TRUE)
      apply_alignment(data = data_pp(),
                      metadata = metadata(),
                      min_frac = input$min_frac,
                      min_samples = input$min_samples,
                      group_by = input$group,
                      plot = FALSE)
    }) %>%
      bindEvent(input$align)

    output$align_plot <- renderPlot({
      color_vector <- create_col_vector(metadata(), color_by = input$group)
      xcms::plotAdjustedRtime(data_aligned(),
                              col = color_vector,
                              xlab="Retention Time (sec)",
                              font.lab=2,
                              cex.lab=2,
                              cex.axis=2,
                              font.main=2,
                              cex.main=2,
                              lwd=2)
      legend("topright",
             legend = unique(names(color_vector)),
             col = unique(color_vector), lty=1)
    }) %>%
      bindEvent(data_aligned())

    data_grouped <- reactive({
      notid <- showNotification('Applying correspondence...',
                                duration = NULL, closeButton = FALSE)
      on.exit(removeNotification(notid), add = TRUE)
      apply_correspondence(data = data_aligned(),
                           metadata = metadata(),
                           min_frac = input$min_frac,
                           min_samples = input$min_samples,
                           group_by = input$group)
    }) %>%
      bindEvent(data_aligned())

    output$is_aligned <- renderText({
      if(is(data_aligned(), 'XCMSnExp')){
        if(xcms::hasAdjustedRtime(data_aligned())){
          'Data has been aligned'
        } else {
          'Please align data'
        }
      } else {
        'Please apply peak picking'
      }

    })

    output$is_grouped <- renderText({
      if(is(data_grouped(), 'XCMSnExp')){
        if(xcms::hasFeatures(data_grouped())){
          'Data has been grouped'
        } else {
          'Please align data'
        }
      } else {
        'Please apply peak picking'
      }

    })

    output$next_buttonSA <- renderUI({
      if(is(data_grouped(), 'XCMSnExp')){
        if(xcms::hasFeatures(data_grouped())){
          actionButton('next_buttonSA', 'Next', class = 'btn-lg btn-success')
        }
      }
    })

    output$back_buttonSA <- renderUI({
      tagList
      actionButton(inputId = 'back_buttonSA',
                   label = 'Back',
                   icon = icon('arrow-left'))

    })

    return(data_grouped)

  })
}
