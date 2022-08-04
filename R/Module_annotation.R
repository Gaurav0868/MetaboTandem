#' Annotation UI
#'`
#' @param id character used to specify namespace, see [`shiny::NS`][shiny::NS()]
#'
#' @return a \code{shiny::\link[shiny]{tagList}} containing UI elements

dbAnnotationUI <- function(id){
  ns <- NS(id)
  tagList(
    h2('Data for annotation'),
    shinyWidgets::radioGroupButtons(ns('data_annot'),
                                    label = 'Select dataset',
                                    choiceNames = c('Current data set',
                                                    'Load new data set'),
                                    choiceValues = c('current_data',
                                                     'new_data'),
                                    checkIcon = list(
                                      yes = tags$i(class = "fa fa-check-square",
                                                   style = "color: steelblue"),
                                      no = tags$i(class = "fa fa-square-o",
                                                  style = "color: steelblue"))),
    uiOutput(ns('load_new_data')),
    br(),
    h2('Annotation databases'),
    shinyDirButton(ns('db_dir'), "Directory with databases",
                   'Select directory'),
    verbatimTextOutput(ns('db_dir_path')),
    box(title = 'Spectra information',
        solidHeader = TRUE,
        status = 'primary',
        selectInput(ns('polarity'),
                    label = 'Polarity',
                    choices = c('positive', 'negative')),
        selectInput(ns('column'),
                    label = 'Column Type',
                    choices = c('hilic', 'rp'))),
    box(title = 'Select databases for annotation',
        solidHeader = TRUE,
        status = 'primary',
        shinyWidgets::materialSwitch(ns('use_custom'),
                                     label = 'Custom Database',
                                     value = FALSE,
                                     status = 'success',
                                     right = TRUE),
        shinyWidgets::materialSwitch(ns('use_gnps'),
                                     label = 'GNPS',
                                     value = FALSE,
                                     status = 'success',
                                     right = TRUE),
        shinyWidgets::materialSwitch(ns('use_masbank'),
                                     label = 'MassBank',
                                     value = FALSE,
                                     status = 'success',
                                     right = TRUE),
        shinyWidgets::materialSwitch(ns('use_mona'),
                                     label = 'MassBank of North America (MoNA)',
                                     value = FALSE,
                                     status = 'success',
                                     right = TRUE),
        shinyWidgets::materialSwitch(ns('use_hmdb'),
                                     label = 'Human Metabolome Database',
                                     value = FALSE,
                                     status = 'success',
                                     right = TRUE),
        br(),
        br(),
        actionButton(ns('annotate'), 'Start Annotation')
    ),

    uiOutput(ns('annot_box'))
  )
}

#' Database annotation server
#'
#' @param id character used to specify namespace, see [`shiny::NS`][shiny::NS()]
#' @param data dataframe with sample information
#'
#' @return
#' \describe{
#'   \item{data_cent_pp}{A [MSnExp-class] object with identified peaks}
#' }

dbAnnotationServer <- function(id, ms1_data, ms2_data){
  moduleServer(id, function(input, output, session){
    ns <- NS(id)

    # Check if new data is required

    output$load_new_data <- renderUI({
      if(input$data_annot == 'new_data'){
        box(title = 'Select new data',
            solidHeader = TRUE,
            width = 12,
            status = 'primary',
            shinyFilesButton(ns('new_features_df'),
                             'Choose feature table (MS1 data)',
                             'Please select file',
                             FALSE),
            br(),
            verbatimTextOutput(ns('new_features_df_file')),
            br(),
            br(),
            shinyFilesButton(ns('new_mgf'),
                             'Choose mgf file (MS2 spectra)',
                             'Please select file',
                             FALSE),
            br(),
            verbatimTextOutput(ns('new_mgf_file'))
        )

      }
    })

    shinyFileChoose(input, 'new_features_df',
                    roots = c('wd' = '.'), session = session)
    feature_file <- reactive({parseFilePaths(roots = c('wd' = '.'),
                                             input$new_features_df)})

    shinyFileChoose(input, 'new_mgf',
                    roots = c('wd' = '.'), session = session)
    mgf_file <- reactive({parseFilePaths(roots = c('wd' = '.'),
                                         input$new_mgf)})

    shinyDirChoose(input, 'db_dir',
                   roots = c('wd' = '.'), session = session)
    db_dir <- reactive({parseDirPath(roots = c('wd' = '.'),
                                     input$db_dir)})

    output$new_features_df_file <- renderPrint({
      paste('Selected MS1 data: ', feature_file()[1])
    })

    output$new_mgf_file <- renderPrint({
      paste('Selected MS2 data: ', mgf_file()[1])
    })

    output$db_dir_path <- renderPrint({
      paste('Database directory: ', db_dir())
    })

    # Annotate data using selected databases


    annotation_custom <- reactive({
      #if(input$use_custom){
        metid::identify_metabolites(
          ms1.data = feature_file()[1],
          ms2.data = mgf_file()[1],
          ms2.match.tol = 0.5,
          ce = "all",
          path = './test_data/outs/',
          ms1.match.ppm = 15,
          rt.match.tol = 30,
          column = "rp",
          candidate.num = 3,
          database = 'snyder_database_rplc0.0.3'
        )
      #}
    }) %>%
      bindEvent(input$annotate)

     if(input$use_gnps){
       param_gnps <- metid::identify_metabolites_params(
         ms1.match.ppm = 15,
         rt.match.tol = 15,
         polarity = input$polarity,
         ce = "all",
         column = input$column,
         total.score.tol = 0.5,
         candidate.num = 3,
         threads = 3,
         database = 'gnps_db'
       )
     }

     if(input$use_massbank){
       param_massbank <- metid::identify_metabolites_params(
         ms1.match.ppm = 15,
         rt.match.tol = 15,
         polarity = input$polarity,
         ce = "all",
         column = input$column,
         total.score.tol = 0.5,
         candidate.num = 3,
         threads = 3,
         database = 'massbank_db'
       )
    }

     if(input$use_mona){
       param_mona <- metid::identify_metabolites_params(
         ms1.match.ppm = 15,
         rt.match.tol = 15,
         polarity = input$polarity,
         ce = "all",
         column = input$column,
         total.score.tol = 0.5,
         candidate.num = 3,
         threads = 3,
         database = 'mona_db'
       )
     }

      params <- c(param_gnps, param_massbank, param_mona)
      params <- params[c(input$use_gnps, input$use_massbank, input$use_mona)]

      db_used <- c('GNPS', 'MassBank', 'MoNA')
      db_used <- db_used[c(input$use_gnps, input$use_massbank, input$use_mona)]

    annotation <- reactive({
      t <- metid::identify_metabolite_all(
        # ms1.data = feature_file()[1],
        # ms2.data = mgf_file()[1],
        ms1.data = 'feature_definitions_annot.csv',
        ms2.data = 'MS2_spectra_consensus.mgf',
        parameter.list = params,
        path = './test_data/outs/'
      )
    }) %>%
      bindEvent(input$annotate)

    annot_table_custom <- reactive({
      if(is(annotation_custom(), 'metIdentifyClass')){
        metid::get_identification_table(annotation_custom())
      }
    }) %>%
      bindEvent(annotation_custom())


    ## Display annotation results

    output$annot_box <- renderUI({
      box(title = 'Annotation Results',
          solidHeader = TRUE,
          status = 'primary',
          width = 12,
          tabsetPanel(
            tabPanel('Custom Database',
                     dataTableOutput(ns('custom_annot'))),
            tabPanel('GNPS',
                     dataTableOutput(ns('gnps_annot'))),
            tabPanel('MassBank',
                     dataTableOutput(ns('massbank_output'))),
            tabPanel('MoNA',
                     dataTableOutput(ns('mona_output')))
          ))
    }) %>%
      bindEvent(input$annotate)

    output$custom_annot <- renderDataTable({
      annot_table_custom()
    })

  })

}
