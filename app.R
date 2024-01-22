#MetaboTandem UI
#'
#' UI for MetaboTandem pipeline
#'
#' @return UI elements for MetaboTandem
#'
#' @param ... UI elements to be passed to MetaboTandem app
#' @param id id for Namespacing
#'
#' @import shiny
#' @import shinydashboard
#' @import shinyFiles
#' @import shinydashboardPlus
#'
#' @export
MetaboTandemApp <- function(){
  ui <- tagList(
    shinyjs::useShinyjs(),
    div(id = 'home',
        style = 'display:none',
        home_ui
    ),
    div(id = 'database_app',
        style = 'display:none',
        database_appUI
    ),
    div(id = 'main_pipeline',
        style = 'display:none',
        ui_main()
    ),
    div(id = 'autotuner_app',
        style = 'display:none',
        autotunerUI('use_autotuner')
    )
  )

  server <- function(input, output, session) {
    shinyjs::show('home')


    observeEvent({input$goHome_database; input$goDatabase}, {
      shinyjs::toggle('database_app')
      shinyjs::toggle('home')

    })

    observeEvent({input$goHome_main; input$goMain}, {
      shinyjs::toggle('main_pipeline')
      shinyjs::toggle('home')
    })

    observeEvent({input$goHome_autotuner; input$goAutotuner}, {
      shinyjs::toggle('autotuner_app')
      shinyjs::toggle('home')
    })

    # Autotuner server

    autotunerServer('use_autotuner')

    # Pre-processing modules

    data <- load_dataServer('load_data')
    data_cent <- peakPickingServer('p_pick', data)
    data_grouped <- alignSpectraServer('align', data$metadata, data_cent)
    data_gap_filled <- gapFillingServer('gap', data_grouped)

    # Optional color picker function

    user_colors <- colorPickerServer('side', data$metadata)

    # Result server for pre-processing
    download_ResultspreprocServer('dl_preproc', data_gap_filled)

    # Annotation
    dbAnnotationServer('annot_dbs', data_gap_filled)
    siriusAnnotationServer('annot_sirius', data_gap_filled)

    # Statistical analysis module
    norm_df <- stastSetupServer('st_setup', data_gap_filled)
    multivariateServer('multi', norm_df, data$metadata, user_colors)
    diffExpressionServer('diffexp', norm_df, data$metadata)

    # Summary results module
    summaryServer('summ', data_gap_filled, data_cent)

  }

  shinyApp(ui, server)
}


