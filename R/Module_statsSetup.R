#' Set options for statistical analysis
#'
#' @param id character used to specify namespace, see [`shiny::NS`][shiny::NS()]
#'
#' @return a \code{shiny::\link[shiny]{tagList}} containing UI elements

statsSetupUI <- function(id){
  ns <- NS(id)
  choice_list <- c('global.norm',
                   'median.norm',
                   'mean.norm',
                   'vsn.norm',
                   'cycloes.norm',
                   'none')
  names(choice_list) <- c('Global sum normalization',
                          'Median normalization',
                          'Mean normalization',
                          'Variance stabilizing normalization',
                          'Cyclic LOESS normalization',
                          'No normalization')
  tagList(
    box(
      title = 'Test normalization methods available',
      solidHeader = TRUE,
      collapsible = TRUE,
      collapsed = TRUE,
      width = 12,
      status = 'primary',
      actionButton(ns('test_norm'), label = 'Test normalization methods'),
      plotOutput(ns('test_plot'))
    ),
    box(
      title = 'Select Normalization method',
      solidHeader = TRUE,
      status = 'primary',
      width = 12,
      shinyWidgets::awesomeRadio(ns('norm_method'),
                                 label = 'Normalization method',
                                 choices = choice_list,
                                 checkbox = TRUE),
      hr(),
      shinyWidgets::materialSwitch(ns('log_transform'),
                                   label = 'Tranform data',
                                   status = 'success'),
      hr(),
      actionButton(ns('apply_norm'), label = 'Apply normalization',
                   class = 'btn-lg btn-success'),
      br(),br(),
      verbatimTextOutput(ns('is_normalized'))
    )
  )
}

stastSetupServer <- function(id, data_proc){
  moduleServer(id, function(input, output, session){

   # Testing normalization
    test_norm_plot <- reactive({
      notid <- showNotification('Testing normalization methods...',
                                duration = NULL, closeButton = FALSE)
      on.exit(removeNotification(notid), add = TRUE)

      features_df <- extract_features(data_proc())
      features_ready <- features_df %>%
        tibble::column_to_rownames(var = 'FeatureID')
      normalize_by_all(features_ready)
    }) %>%
      bindEvent(input$test_norm)

    output$test_plot <- renderPlot(
      test_norm_plot()
    )

    # Applying normalization
    norm_df <- reactive({
      notid <- showNotification('Normalizing data...')
      features_ready <- features_df() %>%
        tibble::column_to_rownames(var = 'FeatureID')
      norm_methods_all(features_ready,
                       input$norm_method,
                       transform = input$log_transform)

    }) %>%
      bindEvent(input$apply_norm)

    output$is_normalized <- renderText({
      if(is.data.frame(norm_df())){
        'Data has been normalized'
      }

    })

    return(norm_df)

  })
}
