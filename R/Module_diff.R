#' Ui for differential expression analysis
#'
#' @param id character used to specify namespace, see [`shiny::NS`][shiny::NS()]
#'
#' @return a \code{shiny::\link[shiny]{tagList}} containing UI elements

diffExpressionUI <- function(id){
  ns <- NS(id)
  tagList(
    box(
      title = 'Treatment groups',
      solidHeader = TRUE,
      collapsible = TRUE,
      width = 12,
      status = 'primary',
      column(width = 4,
             selectizeInput(ns('group'),
                            label = 'Choose grouping variables for differential
                            analysis',
                            choices = NULL),
             actionButton(ns('refresh'), label = 'Refresh')),
      column(width = 4,
             selectizeInput(ns('control'),
                            label = 'Control treatment (baseline)',
                            choices = NULL)),
      column(width = 4,
             selectizeInput(ns('treatment'),
                            label = 'Treatment to compare',
                            choices = NULL)),
      hr(),
      actionButton(ns('calculate'), label = 'Calculate differential expression')
    ),
    box(
      title = 'Differential expression table',
      solidHeader = TRUE,
      width = 12,
      status = 'primary',
      dataTableOutput(ns('diff_table'))
    ),
    box(
      title = 'Volcano Plot',
      solidHeader = TRUE,
      width = 12,
      status = 'primary',
      numericInput(ns('lfc_t'), 'Log2 Fold-Change threshold',
                   min = 0, step = 0.5, value = 2),
      numericInput(ns('pval_t'), 'p-value threshold',
                   min = 0, step = 0.01, value = 0.5),
      shinyWidgets::materialSwitch(ns('p_adjust'), 'Adjust p-value',
                                   status = 'success'),
      hr(),
      plotOutput(ns('volcano'))
    ),

    box(
      title = 'Heatmap',
      solidHeader = TRUE,
      width = 12,
      status = 'primary',
      shinyWidgets::materialSwitch(ns('hmp_sig'), 'Use only significant features',
                                   status = 'success'),
      shinyWidgets::materialSwitch(ns('clus_rows'), 'Cluster Rows',
                                   status = 'success'),
      shinyWidgets::materialSwitch(ns('clus_cols'), 'Cluster Cols',
                                   status = 'success'),
      hr(),
      plotOutput(ns('hmp'))
    )
  )
}

diffExpressionServer <- function(id, norm_df, metadata){
  moduleServer(id, function(input, output, session){

    observe({
      updateVarSelectizeInput(session, 'group', data = metadata(),
                              server = TRUE)
    })

    treatment_list <- reactive({
      unique(metadata() %>%
               dplyr::select(dplyr::all_of(input$group)) %>%
               dplyr::pull())
    }) %>%
      bindEvent(input$refresh)

    #treatment_list <- c('CTR', 'WP')

    observe({
      updateSelectizeInput(session, 'control', choices = treatment_list(),
                           server = TRUE)
      updateSelectizeInput(session, 'treatment', choices = treatment_list(),
                           server = TRUE)
    }) %>%
      bindEvent(treatment_list())


    control_samples <- reactive({
      get_samples(metadata(), Treatment = input$group, value = input$control)
    }) %>%
      bindEvent(input$control)

    treatment_samples <- reactive({
      get_samples(metadata(), Treatment = input$group, value = input$treatment)
    }) %>%
      bindEvent(input$treatment)

    diff_table <- reactive({
      get_diff_table(norm_df(),
                     control.sample_list = control_samples(),
                     treatment.sample_list = treatment_samples(),
                     log2_transformed = TRUE)
    }) %>%
      bindEvent(input$calculate)

    output$diff_table <- renderDataTable(
      diff_table()
    )

    volcano <- reactive({

      if(isTRUE(input$p_adjust)){
        plot_volcano(diff_table(),
                     log2FC,
                     pval.adj,
                     log2FC.threshold = input$lfc_t,
                     pval.threshold = input$pval_t)
      } else{
        plot_volcano(diff_table(),
                     log2FC,
                     pval,
                     log2FC.threshold = input$lfc_t,
                     pval.threshold = input$pval_t)
      }
    }) %>%
      bindEvent(diff_table())

    output$volcano <- renderPlot(
      volcano()
    )

    hmp_matrix <- reactive({
      if(input$hmp_sig){
        norm_df()[diff_table$FeatureID(),]
      } else {
        norm_df()
      }
    })

    output$hmp <- renderPlot({
      mapcolor <- grDevices::colorRampPalette(
        RColorBrewer::brewer.pal(11, 'RdYlBu'))(100)[100:1]
      pheatmap::pheatmap(hmp_matrix,
                         cluster_rows = input$clus_rows,
                         cluster_cols = input$clus_cols,
                         color = mapcolor,
                         show_rownames = FALSE,
                         scale = 'row'
      )
    })

  })
}
