#' Peak picking spectra UI
#'`
#' @param id character used to specify namespace, see [`shiny::NS`][shiny::NS()]
#'
#' @return a \code{shiny::\link[shiny]{tagList}} containing UI elements

peakPickingUI <- function(id){
  ns <- NS(id)

  # Box to select PP method
  tagList(

    headerbox_factory(
      title = 'Peak Picking Method',
      status = 'info',
      width = 12,
      content = tagList(
        selectInput(ns('pp_method'), 'Select method to use:',
                    c('centWave' = 'cw',
                      'Matched Filter' = 'mf',
                      'Massifquant' = 'mq'))
      )
    ),

    # Parameters for selected PP method
    uiOutput(ns('pp_params')),

    headerbox_factory(
      title = '',
      status = 'success',
      width = 6,
      content = tagList(
        uiOutput(ns('plot_pos')),
        uiOutput(ns('table_pos')),
        fluidRow(
          column(6, align = 'center', offset = 3,
                 shinyWidgets::actionBttn(ns('test'),
                                          label = 'Test',
                                          style = 'material-flat',
                                          color = 'warning',
                                          size = 'sm')
                 )
          ),
        br(), br(),
        fluidRow(
          column(6, align = 'center', offset = 3,
                 shinyWidgets::actionBttn(ns('pick'),
                                          label = 'Pick peaks',
                                          style = 'material-flat',
                                          color = 'primary',
                                          size = 'sm')
          )
        ),
        fluidRow(
          verbatimTextOutput(ns('has_peaks'))
        )

        )
      )
  )
}

#' Peak picking server-side processing
#'
#' @param id character used to specify namespace, see [`shiny::NS`][shiny::NS()]
#' @param data dataframe with sample information
#'
#' @return
#' \describe{
#'   \item{data_cent_pp}{A [MSnExp-class] object with identified peaks}
#' }

peakPickingServer <- function(id, data){
  moduleServer(id, function(input, output, session){

    ns <- NS(id)
    cw_params <- tagList(
      h3('Method parameters'),
      fluidRow(
        column(6,
               numericInput(ns('ppm'), 'Ppm threshold', value = 25),
               numericInput(ns('snt'), 'Signal-to-noise threshold', value = 3),
               numericInput(ns('p_width_min'), 'Min. peak width', value = 20),
               numericInput(ns('pf_k'), 'Number of peaks for pre-filtering', value = 3, step = 1)),
        column(6,
               numericInput(ns('noise'), 'Noise threshold', value = 1e6),
               numericInput(ns('mz_diff'), 'Mass difference for overlay peaks', value = 0.01),
               numericInput(ns('p_width_max'), 'Max. peak width', value = 50),
               numericInput(ns('pf_i'), 'Min. Intensity for prefiltering', value = 100))
      )
    )

    mf_params <- tagList(
      h3('Method parameters'),

      fluidRow(
        column(6,
               numericInput(ns('bin'), 'Bin size', value = 0.1),
               numericInput(ns('sigma'), 'sigma', value = 12.72),
               numericInput(ns('steps'), 'Number of bins to be merged', value = 2, step = 1)),
        column(6,
               numericInput(ns('fwhm'), 'fwhm', value = 30),
               numericInput(ns('max'), 'Max. peaks per slice', value = 10, step = 1))
      )
    )

    subset_tags <- tagList(
      hr(),
      h3('Subsetting for testing parameters'),
      sliderInput(ns('rt_range'), 'Retention time range for testing [s]',
                  value = c(0, 240),
                  min = 0,
                  max = 1200),
      sliderInput(ns('mz_range'), 'Range of m/z for testing',
                  value = c(100, 300),
                  min = 0,
                  max = 1200)
    )

    output$pp_params <- renderUI({

      if(input$pp_method != 'mf'){
        cont <- c(cw_params,
                  subset_tags)
      } else {
        cont <- c(mf_params,
                  subset_tags)
      }

      headerbox_factory(
        title = 'Method Parameters',
        status = 'success',
        content = cont
      )
    }) %>%
      bindEvent(input$pp_method)



    params_df <- reactive({
      data.frame(Parameter = c('Noise threshold',
                               'Signal-to-noise ration threshold',
                               'Peak width'),
                 value = c(input$noise,
                           input$snt,
                           paste0(input$p_width_min, '-', input$p_width_max)))
    }) %>%
      bindEvent(input$test)

    test_plot <- reactive({
      test_peak_picking(data$data_cent(),
                        p_width = c(input$p_width_min, input$p_width_max),
                        mz_range = c(input$mz_range[1], input$mz_range[2]),
                        rt_range = c(input$rt_range[1], input$rt_range[2]),
                        snt = input$snt,
                        noise = input$noise)
    }) %>%
      bindEvent(input$test)

    output$test_plot <- renderPlot(
      test_plot()
    )

    output$plot_pos <- renderUI(
      plotOutput(ns('test_plot'))
    ) %>%
      bindEvent(input$test)

    output$table_pos <- renderUI(
      tableOutput(ns('params'))
    ) %>%
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
                         method = input$pp_method,
                         p_width = c(input$p_width_min, input$p_width_max),
                         snt = input$snt,
                         noise = input$noise,
                         ppm = input$ppm,
                         prefilter = c(input$pf_k, input$pf_i),
                         mz_diff = input$mz_diff,
                         bin = input$bin,
                         fwhm = input$fwhm,
                         sigma = input$sigma,
                         max = input$max,
                         steps = input$steps)
    }) %>%
      bindCache(input$p_width_min,
                input$p_width_max,
                input$snt,
                input$noise) %>%
      bindEvent(input$pick)

    output$has_peaks <- renderText({
      if(is(data_cent_pp(), 'XCMSnExp')){
        if(xcms::hasChromPeaks(data_cent_pp())){
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

