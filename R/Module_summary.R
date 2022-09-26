summaryUI <- function(id){
  ns <- NS(id)
  tagList(
    box(title = 'Chromatogram',
        plotOutput(ns('chrom'))),
    shinyjqui::selectableTableOutput(ns('table'),
                                     selection_mode = 'row'),
  )
}

summaryServer <- function(id, data_proc){
  moduleServer(id, function(input, output, session){

    data <- reactive({
      data <- extract_features(data_proc()) %>%
        extract_feature_definition(data_proc(), feature_abundance = .)
    })

    output$table <- renderTable(data())

    output$chrom <- renderPlot({

      i <- input$table_selected
      data_proc() %>%
        MSnbase::filterRt(rt = c(data()$rtmin[i], data()$rtmax[i])) %>%
        MSnbase::filterMz(mz = c(data()$mzmin[i], data()$mzmax[i])) %>%
        MSnbase::chromatogram(., aggregationFun="max") %>%
        xcms::plot(., col = c('blue', 'blue', 'red', 'red'),
                   ylab = "Intensity",
                   xlab = "Retention Time (sec)")
    })
  }
  )
}

