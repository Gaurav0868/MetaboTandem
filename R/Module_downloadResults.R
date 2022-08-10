#' Download results (preprocessing) UI
#'
#' @param id character used to specify namespace, see [`shiny::NS`][shiny::NS()]
#'
#' @return a \code{shiny::\link[shiny]{tagList}} containing UI elements

download_ResultspreprocUI <- function(id){
  ns <- NS(id)
  tagList(
    box(
      title = 'Feature Table',
      solidHeader = TRUE,
      collapsible = TRUE,
      status = 'primary',
      downloadButton(ns('d_feature'), 'Download .tsv',
                     icon = icon('file-export'))

    ),
    box(
      title = 'Feature Definitions',
      solidHeader = TRUE,
      collapsible = TRUE,
      status = 'primary',
      downloadButton(ns('d_feat_def'), 'Download .tsv',
                     icon = icon('file-export'))

    ),
    br(),
    box(
      title = 'Spectra Table',
      solidHeader = TRUE,
      collapsible = TRUE,
      status = 'primary',
      downloadButton(ns('spectra'), 'Download .tsv',
                     icon = icon('file-export'))

    ),
    box(
      title = 'MS2 file',
      solidHeader = TRUE,
      collapsible = TRUE,
      status = 'primary',
      actionButton(ns('extract_ms2'), 'Extract MS2'),
      downloadButton(ns('ms2'), 'Download .mgf',
                     icon = icon('file-export'))

    )
  )
}


download_ResultspreprocServer <- function(id, data_gap_filled){
  moduleServer(id, function(input, output, session){

    bindEvent({
      features_df <- reactive({
        extract_features(data_gap_filled())
      })

      features_definitions_df <- reactive({
        extract_feature_definition(data_gap_filled(),
                                   feature_abundance = features_df())
      })

      spectra_df <- reactive({
        extract_spectra_table(data_gap_filled())
      })
    },
    data_gap_filled())

    ms2_consensus <- reactive({
      notid <- showNotification('Extracting MS2...',
                                duration = NULL, closeButton = FALSE)
      on.exit(removeNotification(notid), add = TRUE)
      extract_MS2_consensus(data_gap_filled())
    }) %>%
      bindEvent(input$extract_ms2)

    output$d_feature <- downloadHandler(
      filename = function(){
        paste('feature_table.tsv')
      },
      content = function(file){
        vroom::vroom_write(features_df(), file)
      }
    )

    output$d_feat_def <- downloadHandler(
      filename = function(){
        paste('feature_definitions.tsv')
      },
      content = function(file){
        vroom::vroom_write(features_definitions_df(), file)
      }
    )

    output$spectra <- downloadHandler(
      filename = function(){
        paste('spectra_table.tsv')
      },
      content = function(file){
        vroom::vroom_write(spectra_df(), file)
      }
    )

    output$ms2 <- downloadHandler(
      filename = function(){
        paste('spectra_consensus.mgf')
      },
      content = function(file){
        mod_writeMgfDataFile(ms2_consensus(), file)
      }
    )


  })
}
