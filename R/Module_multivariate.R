#' Perform multivariate Analysis
#'
#' @param id character used to specify namespace, see [`shiny::NS`][shiny::NS()]
#'
#' @return a \code{shiny::\link[shiny]{tagList}} containing UI elements

multivariateUI <- function(id){
  ns <- NS(id)
  tagList(
    box(
      title = 'Options',
      solidHeader = TRUE,
      collapsible = TRUE,
      width = 12,
      status = 'primary',
      column(width = 6,
             selectizeInput(ns('group'),
                            label = 'Choose grouping variables for plot',
                            choices = NULL)),
      column(width = 6,
             shinyWidgets::materialSwitch(ns('label'),
                                          label = 'Toogle sample names',
                                          status = 'primary')),
      hr(),
      actionButton(ns('calculate'), label = 'Calculate ordination')
    ),
    box(
      title = 'NMDS',
      solidHeader = TRUE,
      status = 'primary',
      plotOutput(ns('stress')),
      plotOutput(ns('nmds'))
    ),
    box(
      title = 'PCA',
      solidHeader = TRUE,
      status = 'primary',
      plotOutput(ns('scree')),
      plotOutput(ns('PCA'))
    ),
    box(
      title = 'PERMANOVA',
      solidHeader = TRUE,
      width = 12,
      status = 'primary',
      tableOutput(ns('permanova'))
    )
  )
}

multivariateServer <- function(id, norm_df, metadata){
  moduleServer(id, function(input, output, session){

    observe({
      updateVarSelectizeInput(session, 'group', data = metadata(),
                              server = TRUE)
    })


    # Calculate multivariate statistics

    nmds_ord <- reactive({
      nmds_ordination(norm_df(), metadata(), mode = 'ra', color_by = input$group)
    }) %>%
      bindEvent(input$calculate)

    pca_ord <- reactive({
      pca_ordination(norm_df(), metadata(), color_by = input$group)
    }) %>%
      bindEvent(input$calculate)

    permanova_table <- reactive({
      calculate_permanova(norm_df(), metadata(), mode = 'ra',
                          group_by =  input$group)
    }) %>%
      bindEvent(input$calculate)

    # Calculate outputs

    output$stress <- renderPlot(
      vegan::stressplot(nmds_ord()$nmds)
    ) %>%
      bindEvent(nmds_ord())

    output$nmds <- renderPlot(
      if(input$label){
        nmds_ord()$nmds_plot +
          geom_label(aes(label = SampleID))
      } else {
        nmds_ord()$nmds_plot
      }
    ) %>%
      bindEvent(nmds_ord())

    output$scree <- renderPlot(
      pca_ord()$scree_plot
    ) %>%
      bindEvent(pca_ord())

    output$PCA <- renderPlot(
      if(isTRUE(input$label)){
        pca_ord()$pca_plot +
          geom_label(aes(label = SampleID))
      } else {
        pca_ord()$pca_plot
      }
    ) %>%
      bindEvent(pca_ord())

    output$permanova <- renderTable(
      permanova_table() %>%
        tibble::rownames_to_column(var = 'treatment'),
      striped = TRUE,
      bordered = TRUE
    )
  })
}

