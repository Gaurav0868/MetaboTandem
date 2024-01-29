#' Annotation UI
#'`
#' @param id character used to specify namespace, see [`shiny::NS`][shiny::NS()]
#'
#' @return a \code{shiny::\link[shiny]{tagList}} containing UI elements

dbAnnotationUI <- function(id){
  ns <- NS(id)
  tagList(
    uiOutput(ns('load_new_data')),
    br(),

    headerbox_factory(
      title = 'Select databases for annotation',
      status = 'primary',
      width = 6,
      content = tagList(
        shinyDirButton(ns('db_dir'), "Directory with databases",
                       'Select directory'),
        hr(),
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
        textOutput(ns('n_dbs'))
      )
    ),

    headerbox_factory(
      title = 'Annotation options',
      status = 'info',
      width = 6,
      content = tagList(
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
        uiOutput(ns('annotation_button'))
      )
    ),

    #uiOutput(ns('annot_order')),

    uiOutput(ns('annot_box'))
  )
}

#' Database annotation server
#'
#' @param id character used to specify namespace, see [`shiny::NS`][shiny::NS()]
#' @param ms1_data A [data.frame] with columns: 'name', 'mz' and 'rt'
#'

dbAnnotationServer <- function(id, data_proc, new_data = FALSE){
  moduleServer(id, function(input, output, session){
    ns <- NS(id)

    # Check if new data is required

    output$load_new_data <- renderUI({
      if(new_data == TRUE){
        tagList(
          h2('Annotation App'),
          headerbox_factory(title = 'Select new data',
                            width = 12,
                            status = 'primary',
                            content = tagList(
                              fileInput(ns('new_features_df'),
                                        'Upload MS1 data (feature table)',
                                        accept = c('.csv', '.tsv')),
                              hr(),
                              fileInput(ns('new_mgf_file'),
                                        'Upload MS2 data (mgf file)',
                                        accept = c('.mgf')),
                            )

          )
        )

      }
    })

    options(shiny.maxRequestSize=30*1024^2)
    if(new_data == TRUE){
      feature_table <- reactiveVal()

      observe({
        req(input$new_features_df)

        ext <- tools::file_ext(input$new_features_df$name)
        if(ext == 'csv'){
          feature_table(vroom::vroom(input$new_features_df$datapath, delim = ","))
        } else if(ext == 'tsv'){
          feature_table(vroom::vroom(input$new_features_df$datapath, delim = "\t"))
        } else{
          validate("Invalid file; Please upload a .csv or .tsv file")
        }
      })

      ms2_data <- reactiveVal()

      observe({
        req(input$new_mgf_file)

        ext <- tools::file_ext(input$new_mgf_file$name)
        if(ext == 'mgf'){
          ms2_data(read_mgf(input$new_mgf_file$datapath))
        } else{
          validate("Invalid file; Please upload a MGF file")
        }
      })
    }

    shinyDirChoose(input, 'db_dir',
                   roots = c('wd' = '.'), session = session)
    db_dir <- reactive({parseDirPath(roots = c('wd' = '.'),
                                     input$db_dir)})

    output$db_dir_path <- renderPrint({
      paste('Database directory: ', db_dir())
    })

    # Annotate data using selected databases

    # ms1.data <- reactive({
    #   ms1.data <- extract_features(data_proc()) %>%
    #     extract_feature_definition(data_proc(), feature_abundance = .)
    # }) %>%
    #   bindEvent(input$extract)
    #
    # ms2.data <- reactive({
    #   notid <- showNotification('Extracting MS2...',
    #                             duration = NULL, closeButton = FALSE)
    #   on.exit(removeNotification(notid), add = TRUE)
    #   suppressWarnings(extract_MS2_consensus(data_proc()))
    # }) %>%
    #   bindEvent(input$extract)
    #
    # output$ms2_data_check <- renderPrint({
    #   if(is(ms2.data(), 'MSpectra')){
    #     'MS2 data extracted succesfully'
    #   }
    # })

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

    output$annotation_button <- renderUI({
      actionButton(ns('annotate'), 'Start Annotation')
    }) %>%
      bindEvent(selected_dbs())

    annotation <- reactiveVal()

    observe({
      req(input$annotate)

      notid <- showNotification('Annotating with selected databases...',
                                duration = NULL, closeButton = FALSE)
      on.exit(removeNotification(notid), add = TRUE)
      annotation(mod_identify_all(feature_table(), ms2_data(), db_dir(), parameter.list()))
    })

    output$annot_box <- renderUI({
      box(title = 'Annotation Results',
          solidHeader = TRUE,
          status = 'primary',
          width = 12,
          tabsetPanel(id = ns('annot_tabs'),
                      tabPanel('...',
                               h3('Annotation successfull.\nCheck results using the button below'),
                               br(),
                               actionButton(ns('get_table'), 'Get annotation_tables')))
      )
    }) %>%
      bindEvent(annotation())


    observeEvent(input$get_table, {
      tabname <- stringr::str_to_title(stringr::str_replace(selected_dbs(), '_db', ' Database'))

      purrr::map2(tabname, selected_dbs(), function(x,y){
        insertTab(inputId = 'annot_tabs',
                  tabPanel(x, dataTableOutput(ns(y))))
      })

      annot_results <- reactive({
        res <- purrr::map(annotation(), function(x){
          mod_get_identification_table(x) %>%
            filter(!is.na(Compound.name))
        })
        return(res)
      })

      # TODO avoid cloning results

      for(res in selected_dbs()){
        output[[res]] <- renderDataTable(
          annot_results()[[res]],
          options = list(scrollX = TRUE,
                         dom = 'ltip')
        )
      }
    })

  })
}
