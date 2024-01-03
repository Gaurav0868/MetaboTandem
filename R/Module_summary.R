summaryUI <- function(id){
  ns <- NS(id)
  tagList(
    box(title = 'Chromatogram',
        plotOutput(ns('chrom'))),
    shinyjqui::selectableTableOutput(ns('table'),
                                     selection_mode = 'row'),
  )
}

summaryServer <- function(id, data_proc, data_align){
  moduleServer(id, function(input, output, session){

    data <- reactive({
      data <- extract_features(data_proc()) %>%
        extract_feature_definition(data_proc(), feature_abundance = .)
    })

    output$table <- renderTable(data())

    output$chrom <- renderPlot({

      i <- input$table_selected

      mzr <- c(data()$mzmin[i] - 0.01, data()$mzmax[i] + 0.01)

      if(data()$rtmax[i] - data()$rtmin[i] < 5){
        rtr <- c(data()$rtmin[i]- 2.5, data()$rtmax[i] + 2.5)
      } else {
        rtr <- c(data()$rtmin[i], data()$rtmax[i])
      }


      data_align() %>%
        MSnbase::filterRt(rt = rtr) %>%
        MSnbase::filterMz(mz = mzr) %>%
        MSnbase::chromatogram(., aggregationFun="max") %>%
        xcms::plot(., col = c('blue', 'blue', 'red', 'red'),
                   ylab = "Intensity",
                   xlab = "Retention Time (sec)")
    })
  }
  )
}

