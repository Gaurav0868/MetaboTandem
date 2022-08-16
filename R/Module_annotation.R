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
        shinyWidgets::materialSwitch(ns('use_massbank'),
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
        actionButton(ns('set_dbs'), 'Set databases'),
        hr(),
        verbatimTextOutput(ns('n_dbs'))
    ),

    box(title = 'Annotation options',
        solidHeader = TRUE,
        status = 'primary',
        selectInput(ns('polarity'),
                    label = 'Polarity',
                    choices = c('positive', 'negative')),
        selectInput(ns('column'),
                    label = 'Column Type',
                    choices = c('hilic', 'rp')),
        numericInput(ns('candidates'),
                     label = 'Number of candidate annotations',
                     min = 1,
                     step = 1,
                     value = 1),
        br(),
        actionButton(ns('extract'), 'Extract MS2'),
        br(),
        verbatimTextOutput(ns('ms2_data_check')),
        br(),
        uiOutput(ns('annotation_button')),
        uiOutput(ns('tables_button'))),

    uiOutput(ns('annot_box')),
    verbatimTextOutput(ns('test_annot'))
  )
}

#' Database annotation server
#'
#' @param id character used to specify namespace, see [`shiny::NS`][shiny::NS()]
#' @param ms1_data A [data.frame] with columns: 'name', 'mz' and 'rt'
#'

dbAnnotationServer <- function(id, data_proc){
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

    ms1.data <- reactive({
      ms1.data <- extract_features(data_proc()) %>%
        extract_feature_definition(data_proc(), feature_abundance = .)
    }) %>%
      bindEvent(input$extract)

    ms2.data <- reactive({
      notid <- showNotification('Extracting MS2...',
                                duration = NULL, closeButton = FALSE)
      on.exit(removeNotification(notid), add = TRUE)
      suppressWarnings(extract_MS2_consensus(data_proc()))
    }) %>%
      bindEvent(input$extract)

    output$ms2_data_check <- renderPrint({
      if(is(ms2.data(), 'MSpectra')){
        'MS2 data extracted succesfully'
      }
    })

    output$annotation_button <- renderUI({
      actionButton(ns('annotate'), 'Start Annotation')
    }) %>%
      bindEvent(ms2.data())

    available_dbs <- c('custom_db', 'gnps_db', 'massbank_db', 'mona_db', 'hmdb_db')

    selected_dbs <- reactive({
      available_dbs[c(input$use_custom,
                      input$use_gnps,
                      input$use_massbank,
                      input$use_mona,
                      input$use_hmdb)]
    }) %>%
      bindEvent(input$set_dbs)

    output$n_dbs <- renderPrint({
      n <- length(selected_dbs())
      paste(n, 'databases selected')
    }) %>%
      bindEvent(input$set_dbs)

    parameter.list <- reactive({
      p <- purrr::map(selected_dbs(), function(x){
        metid::identify_metabolites_params(
          ms1.match.ppm = 15,
          rt.match.tol = 15,
          polarity = input$polarity,
          ce = "all",
          column = input$column,
          total.score.tol = 0.5,
          candidate.num = 3,
          threads = parallelly::availableCores(),
          database = x)
      })
      names(p) <- selected_dbs()
      return(p)
    }) %>%
      bindEvent(selected_dbs())

    annotation <- reactive({
      notid <- showNotification('Annotating with selected databases...',
                                duration = NULL, closeButton = FALSE)
      on.exit(removeNotification(notid), add = TRUE)
      mod_identify_all(ms1.data(), ms2.data(), db_dir(), parameter.list())
    }) %>%
      bindEvent(input$annotate)

    output$tables_button <- renderUI({

    }) %>%
      bindEvent(annotation())

    output$annot_box <- renderUI({
      box(title = 'Annotation Results',
          solidHeader = TRUE,
          status = 'primary',
          width = 12,
          tabsetPanel(id = ns('annot_tabs'),
                      tabPanel('...',
                               actionButton(ns('get_table'), 'Get annotation_tables')))
      )
    }) %>%
      bindEvent(annotation())


    observeEvent(input$get_table, {
      tabname <- stringr::str_to_title(stringr::str_replace(selected_dbs(), '_db', 'Database'))
      table_id <- stringr::str_replace(selected_dbs(), '_db', '_input')

      purrr::map2(tabname, table_id, function(x,y){
        insertTab(inputId = 'annot_tabs',
                  tabPanel(x,
                           dataTableOutput(ns(y))))
      })
    })


    # insertTab(inputId = 'annot_tabs',
    #           tab = tabPanel('dd', h2('ha')),
    #           target = 'Results')


    # ))
    # })

    # annot_tables <- reactive({
    #   custom_annot
    # })
    #
    # output$custom_output<- renderDataTable({
    #   if(input$use_custom){
    #     mod_get_identification_table(annotation()[['custom_db']])
    #   } else {
    #     tibble(Database = '', not_selected = '')
    #   }
    # }) %>%
    #   bindEvent(annotation())
    #
    # output$custom_output<- renderDataTable({
    #   if(input$use_custom){
    #     mod_get_identification_table(annotation()[['custom_db']])
    #   } else {
    #     tibble(Database = '', not_selected = '')
    #   }
    # }) %>%
    #   bindEvent(annotation())
    #
    # output$gnps_output <- renderDataTable({
    #   if(input$use_gnps){
    #     mod_get_identification_table(annotation()[['gnps_db']])
    #   } else {
    #     tibble(Database = '', not_selected = '')
    #   }
    # }) %>%
    #   bindEvent(annotation())
    #
    # output$mona_output <- renderDataTable({
    #   if(input$use_mona){
    #     mod_get_identification_table(annotation()[['mona_db']])
    #   } else {
    #     tibble(Database = '', not_selected = '')
    #   }
    # }) %>%
    #   bindEvent(annotation())
    # #
    # output$massbank_output <- renderDataTable({
    #   if(input$use_massbank){
    #     mod_get_identification_table(annotation()[['massbank_db']])
    #   } else {
    #     tibble(Database = '', not_selected = '')
    #   }
    # }) %>%
    #   bindEvent(annotation())
    #
    # output$hmdb_output <- renderDataTable({
    #   if(input$use_hmdb){
    #     mod_get_identification_table(annotation()[['hmdb']])
    #   } else {
    #     tibble(Database = '', not_selected = '')
    #   }
    # }) %>%
    #   bindEvent(annotation())

    output$test_annot <- renderPrint({
      names(annotation())
    })


  })
}
