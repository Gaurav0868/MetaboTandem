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
#'
#' @export

MetaboTandemApp <- function(){
  ui <- dashboardPage(
    dashboardHeader(
      title = 'MetaboTandem',
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
        menuItem('Load Data',
                 tabName = 'load_data',
                 icon = icon('upload')),
        menuItem('Peak picking',
                 tabName = 'p_pick',
                 icon = icon('check')),
        menuItem('Alignment & Correspondence',
                 tabName = 'align',
                 icon = icon('align-center')),
        menuItem('Gap Filling',
                 tabName = 'gap',
                 icon = icon('fill')),
        menuItem('Results Download',
                 tabName = 'res',
                 icon = icon('file-download'))
      )
    ),

    # Body content
    dashboardBody(
      tabItems(
        # Load data tab
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
        tabItem(tabName = 'res',
                h1('Select data to download'))
        )
      )
    )

  server <- function(input, output, session) {
    data <- load_dataServer('load_data')
    data_cent_pp <- peakPickingServer('p_pick', data)
    data_grouped <- alignSpectraServer('align', data$metadata, data_cent_pp)
    gg <- gapFillingServer('gap', data_grouped)
  }

  shinyApp(ui, server)
}


