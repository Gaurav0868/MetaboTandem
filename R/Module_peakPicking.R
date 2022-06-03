#' Peak picking spectra UI
#'`
#' @param id, character used to specify namespace, see [`shiny::NS`][shiny::NS()]
#'
#' @return a \code{shiny::\link[shiny]{tagList}} containing UI elements

peakPickingUI <- function(id){
  ns <- NS(id)
  tagList(
    h2('Peak picking parameters'),
    p('Test parameters for peak picking using the',
      strong('Test Parameters'),
      'button.'),
    p('Apply desired parameters to all data using the',
      strong('Apply peak picking'),
      'button.'),
    fluidRow(
      column(6,
             numericInput(ns('noise'), 'Noise threshold', value = 0),
             numericInput(ns('snt'), 'Signal-to-noise threshold', value = 0),
             sliderInput(ns('p.width'), 'Peak Width',
                         value = c(20, 50),
                         min = 0,
                         max = 100),
             sliderInput(ns('rt_range'), 'Retention time range for testing [s]',
                         value = c(0, 60),
                         min = 0,
                         max = 1200),
             sliderInput(ns('mz_range'), 'Range of m/z for testing',
                         value = c(100, 200),
                         min = 0,
                         max = 1200)),
      column(6,
             actionButton(ns('test'), 'Test Parameters'),
             br(),
             tableOutput(ns('params')),
             br(), br(),
             actionButton(ns('pick'), 'Apply peak picking',
                          class = 'btn-lg btn-success'),
             br(), br(),
             verbatimTextOutput(ns('has_peaks')))
    ),
    waiter::use_waiter(),
    plotOutput(ns('test_plot'))
  )
}

#' Peak picking server-side processing
#'
#' @param id character used to specify namespace, see [shiny::NS][shiny::NS()]
#' @param data dataframe with sample information
#'
#' @return
#' \describe{
#'   \item{data_cent_pp}{A [MSnExp-class] object with identified peaks}
#' }

peakPickingServer <- function(id, data){
  moduleServer(id, function(input, output, session){

    params_df <- reactive({
      data.frame(Parameter = c('Noise threshold',
                               'Signal-to-noise ration threshold',
                               'Peak width'),
                 value = c(input$noise,
                           input$snt,
                           paste0(input$p.width[1], '-', input$p.width[2])))
    }) %>%
      bindEvent(input$test)

    test_plot <- reactive({
      waiter::Waiter$new(id = 'test_plot',
                         html = waiter::spin_folding_cube())$show()
      test_peak_picking(data$data_cent(),
                        p.width = c(input$p.width[1], input$p.width[2]),
                        mz.range = c(input$mz_range[1], input$mz_range[2]),
                        rt.range = c(input$rt_range[1], input$rt_range[2]),
                        snt = input$snt,
                        noise = input$noise)
    }) %>%
      bindEvent(input$test)

    output$test_plot <- renderPlot(
      test_plot()
    )

    output$params <- renderTable(
      params_df(),
      striped = TRUE,
      bordered = TRUE
    )

    data_cent_pp <- reactive({
      notid <- showNotification('Detecting peaks...',
                                duration = NULL, closeButton = FALSE)
      on.exit(removeNotification(notid), add = TRUE)
      apply_peak_picking(data$data_cent(),
                         p.width = c(input$p.width[1], input$p.width[2]),
                         snt = input$snt,
                         noise = input$noise)
    }) %>%
      bindCache(input$p.width,
                input$snt,
                input$noise) %>%
      bindEvent(input$pick)

    output$has_peaks <- renderText({
      if(is(data_cent_pp(), 'XCMSnExp')){
        if(hasChromPeaks(data_cent_pp())){
          'Peaks have been identified'
        } else {
          'Please apply peak picking'
        }
      } else {
        'Please apply peak picking'
      }

    })

    return(data_cent_pp)

  })

}

