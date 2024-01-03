#' MetaboTandem UI
#'
#' UI for MetaboTandem pipeline
#'
#' @return UI elements for MetaboTandem
#'
#' @param ... UI elements to be passed to MetaboTandem app
#' @param id id for Namespacing
#'
#' @import shiny
#' @import shinydashboard
#' @import shinyFiles
#' @import shinydashboardPlus
#'
#' @export

old_MetaboTandemApp <- function(){
  ui <- dashboardPage(
    dashboardHeader(
      title = tagList(
        span(class = "logo-lg", "MetaboTandem"),
        img(src = "logo.png")),
      dropdownMenu(
        type = 'notifications',
        icon = icon('question-circle'),
        headerText = 'Help',

        notificationItem('Github Repository',
                         icon = icon('github'),
                         href = 'https://github.com/Coayala/MetaboTandem'),
        notificationItem('User Guide',
                         icon = icon('file'),
                         href = 'https://github.com/Coayala/MetaboTandem')
      )
    ),

    # Sidebar content
    dashboardSidebar(
      sidebarMenu(
        menuItem('Data pre-processing',
                 tabName = 'preproc',
                 icon = icon('cogs'),
                 menuSubItem('Load Data',
                             tabName = 'load_data',
                             icon = icon('upload')),
                 menuSubItem('Peak picking',
                             tabName = 'p_pick',
                             icon = icon('check')),
                 menuSubItem('Alignment & Correspondence',
                             tabName = 'align',
                             icon = icon('align-center')),
                 menuSubItem('Gap Filling',
                             tabName = 'gap',
                             icon = icon('fill')),
                 startExpanded = TRUE),
        menuItem('Annotation',
                 tabName = 'annot',
                 icon = icon('tags'),
                 menuSubItem('Public or custom databases',
                             tabName = 'dbs_annot',
                             icon = icon('database')),
                 menuSubItem('Using SIRIUS',
                             tabName = 'sirius_annot',
                             icon = icon('computer')),
                 startExpanded = TRUE),
        menuItem('Statistical Analysis',
                 tabName = 'stat',
                 icon = icon('chart-line'),
                 menuSubItem('Setup',
                             tabName = 'stats-setup',
                             icon = icon('tasks')),
                 menuSubItem('Multivariate analysis',
                             tabName = 'stats-mult',
                             icon = icon('chart-bar')),
                 menuSubItem('Differential expression',
                             tabName = 'diff-exp',
                             icon = icon('expand-alt')),
                 startExpanded = TRUE),
        menuItem('Results summary',
                 tabName = 'res_summ',
                 icon = icon('table')),
        menuItem('Results Download',
                 tabName = 'res',
                 icon = icon('download'),
                 menuSubItem('Pre-processing Tables',
                             tabName = 'res_preproc',
                             icon = icon('file-download')))
      )
    ),

    # Body content
    dashboardBody(
      tags$head(
        tags$style(HTML("
      #sidebarItemExpanded > ul > :last-child {
        position: absolute;
        bottom: 0;
        width: 100%;
        background-color: steelblue;
      }

    "))),
      tabItems(
        # Load data tab

        ## Pre-processing tabs
        tabItem(tabName = 'load_data',
                h1('Load your data'),
                load_dataUI('load_data')),
        tabItem(tabName = 'p_pick',
                h1('Peak Picking'),
                peakPickingUI('p_pick')),
        tabItem(tabName = 'align',
                h1('Spectra alignment'),
                alignSpectraUI('align')),
        tabItem(tabName = 'gap',
                h1('Gap Filling'),
                gapFillingUI('gap')),

        ## Annotation module
        tabItem(tabName = 'dbs_annot',
                h1('Annotation using custom or public databases'),
                dbAnnotationUI('annot_dbs')),
        tabItem(tabName = 'sirius_annot',
                h1('Annotation using SIRIUS predictions'),
                siriusAnnotationUI('annot_sirius')),

        ## Statistical analysis
        tabItem(tabName = 'stats-setup',
                h1('Set options for statistical analysis'),
                statsSetupUI('st_setup')),
        tabItem(tabName = 'stats-mult',
                h1('Multivariate Analysis'),
                multivariateUI('multi')),
        tabItem(tabName = 'diff-exp',
                h1('Differential Analysis'),
                diffExpressionUI('diffexp')),

        ## Results tabs
        tabItem(tabName = 'res_summ',
                h1('Results summary'),
                summaryUI('summ')),

        tabItem(tabName = 'res_preproc',
                h1('Select data to download'),
                download_ResultspreprocUI('dl_preproc'))
      )
    ),

    # Sidebar contet
    dashboardControlbar(
      br(),
      box(
        title = 'Color palette selector',
        solidHeader = TRUE,
        width = 12,
        colorPickerUI('side')
      )
    )
  )

  server <- function(input, output, session) {

    # Pre-processing modules

    data <- load_dataServer('load_data')
    data_cent <- peakPickingServer('p_pick', data)
    data_grouped <- alignSpectraServer('align', data$metadata, data_cent)
    data_gap_filled <- gapFillingServer('gap', data_grouped)

    # Optional color picker function

    user_colors <- colorPickerServer('side', data$metadata)

    # Result server for pre-processing
    download_ResultspreprocServer('dl_preproc', data_gap_filled)

    # Annotation
    dbAnnotationServer('annot_dbs', data_gap_filled)
    siriusAnnotationServer('annot_sirius', data_gap_filled)

    # Statistical analysis module
    norm_df <- stastSetupServer('st_setup', data_gap_filled)
    multivariateServer('multi', norm_df, data$metadata, user_colors)
    diffExpressionServer('diffexp', norm_df, data$metadata)

    # Summary results module
    summaryServer('summ', data_gap_filled, data_cent)

  }

  shinyApp(ui, server)
}


