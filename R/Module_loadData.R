#' Metadata uploading UI
#'
#' @param id character used to specify namespace,
#' see [`shiny::NS`]([shiny::NS()])
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
#' @param id character used to specify namespace,
#' see \code{shiny::\link[shiny]{NS}}
#'
#' @return A [data.frame] with sample information
#'

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
      metadata(),
      options = list(scrollX = TRUE)
    )
    return(metadata)
  })
}

#' Spectra uploading UI
#'
#' @param id character used to specify namespace, see [`shiny::NS`][shiny::NS()]
#'
#' @return a \code{shiny::\link[shiny]{tagList}} containing UI elements

spectraUI <- function(id){
  ns <- NS(id)
  tagList(
    waiter::use_waiter(),
    h2('Select data directory'),
    shinyDirButton(ns('datadir'), 'Choose data folder',
                      'Please select folder with data',
                   FALSE),
    br(),
    verbatimTextOutput(ns('sel_directory')),
    br(),
    selectInput(ns('format'), 'Data format', c('mzML', 'mzXML')),
    br(),
    verbatimTextOutput(ns('is_loaded')),
    actionButton(ns('load'), 'Load Data', class = 'btn-primary')
  )
}

#' Spectra uploading server-side processing
#'
#' @param id character used to specify namespace, see [`shiny::NS`][shiny::NS()]
#' @param metadata dataframe with sample information
#'
#' @return A centroided [MSnExp-class] object

spectraServer <- function(id, metadata){
  moduleServer(id, function(input, output, session){

    #volumes <- getVolumes()

    shinyDirChoose(input, 'datadir', roots = c('wd' = '.'), session = session)

    data_dir <- reactive({parseDirPath(roots = c('wd' = '.'), input$datadir)})

    output$sel_directory <- renderPrint(
      paste('Directory selected: ', data_dir())
    )

    data_cent <- reactiveVal()

    observe({
      req(input$load)
      notid <- showNotification('Reading data...',
                                duration = NULL, closeButton = FALSE)
      on.exit(removeNotification(notid), add = TRUE)
      temp_metadata <- metadata()
      temp_data <- load_spectra_data(data_dir(),
                                     temp_metadata,
                                     input$format)
      data_cent(centroid_check(temp_data,
                               transform = TRUE))
    })

    output[['is_loaded']] <- renderText({
      if(is(data_cent(), 'OnDiskMSnExp')){
        'Data loaded correctly'
      } else {
        'Please load your data'
      }
    })

    return(data_cent)
  })
}

#' Load data UI
#'
#' @param id character used to specify namespace, see [`shiny::NS`][shiny::NS()]
#'
#' @return a \code{shiny::\link[shiny]{tagList}} containing UI elements

load_dataUI <- function(id){
  ns <- NS(id)
  # fluidRow(
  #   column(6,
  #     metadataUI(ns('metadata'))
  #   ),
  #   column(6,
  #     spectraUI(ns('spectra'))
  #   )
  # )
  fluidRow(
    box(
      metadataUI(ns('metadata'))
    ),
    box(
      spectraUI(ns('spectra'))
    )
  )

}

#' Load data server-side processing
#'
#' @param id character used to specify namespace, see [`shiny::NS`][shiny::NS()]
#'
#' @return list with following components
#'   \itemize{
#'     \item {`data_cent` A centroided [MSnExp-class] object}
#'     \item {`metadata` A dataframe with sample information}
#'   }

load_dataServer <- function(id){
  moduleServer(id, function(input, output, session){
    metadata <- metadataServer('metadata')
    data_cent <- spectraServer('spectra', metadata)

    return(
      list(data_cent = data_cent,
           metadata = metadata)
    )
  })
}
