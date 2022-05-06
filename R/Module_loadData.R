#' Metadata uploading UI
#'
#' @param id, character used to specify namespace, see \code{shiny::\link[shiny]{NS}}
#'
#' @return a \code{shiny::\link[shiny]{tagList}} containing UI elements

metadataUI <- function(id){
  ns <- NS(id)
  tagList(
    h2('Select metadata file'),
    fileInput(ns('metadata_file'), 'Choose file', accept = c('.csv', '.tsv')),
    dataTableOutput(ns('metadata_table'))
  )
}

#' Metadata uploading  server-side processing
#'
#' @param id, character used to specify namespace, see \code{shiny::\link[shiny]{NS}}
#' @param input, output, session standard \code{shiny} boilerplate
#'
#' @return list with following components
#' \describe{
#'   \item{xvar}{reactive character string indicating x variable selection}
#'   \item{yvar}{reactive character string indicating y variable selection}
#' }

metadataServer <- function(id){
  moduleServer(id, function(input, output, session){

    metadata <- reactiveVal()

    observe({
      req(input$metadata_file)

      ext <- tools::file_ext(input$metadata_file$name)
      if(ext == 'csv'){
        metadata(vroom::vroom(input$metadata_file$datapath, delim = ","))
      } else if(ext == 'tsv'){
        metadata(vroom::vroom(input$metadata_file$datapath, delim = "\t"))
      } else{
        validate("Invalid file; Please upload a .csv or .tsv file")
      }
    })

    output$metadata_table <- renderDataTable(
      metadata()
    )
    return(metadata)
  })
}

spectraUI <- function(id){
  ns <- NS(id)
  tagList(
    h2('Select data directory'),
    shinyDirButton(ns('datadir'), 'Choose data folder',
                      'Please select folder with data',
                   FALSE),
    br(),
    verbatimTextOutput(ns('sel_directory')),
    br(),
    actionButton(ns('load'), 'Load Data'),
    br(),
    verbatimTextOutput(ns('is_loaded')),
    br()
  )
}

spectraServer <- function(id, metadata){
  moduleServer(id, function(input, output, session){

    #volumes <- getVolumes()

    shinyDirChoose(input, 'datadir', roots = c('wd' = '.'), session = session)

    data_dir <- reactive({parseDirPath(roots = c('wd' = '.'), input$datadir)})

    output$sel_directory <- renderPrint(
      paste('Directory selected: ', data_dir())
    )

    data <- eventReactive(input$load, {
      temp_met <- metadata()
      temp_data <- load_spectra_data(data_dir(),
                                     temp_met,
                                     input$format)
      centroid_check(temp_data,
                     transform = TRUE)
    })
    output$is_loaded <- renderText(
      if(is(data(), 'OnDiskMSnExp')){
        'Data loaded correctly'
      } else {
        'Please load your data'
      }
    )

    return(data)
  })
}

load_dataUI <- function(id){
  ns <- NS(id)
  fluidRow(
    column(6,
      metadataUI(ns('met1'))
    ),
    column(6,
      spectraUI(ns('spec1'))
    )
  )
}

load_dataServer <- function(id){
  moduleServer(id, function(input, output, session){
    metadata <- metadataServer('met1')
    data <- spectraServer('spec1', metadata)

    return(
      list(data = reactive({data}),
           metadata = reactive({metadata}))
    )
  })
}
