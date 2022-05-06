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
        menuItem('Load Data', tabName = 'load_data', icon = icon('upload')),
        menuItem('Spectra alignment', tabName = 'align', icon = icon('align-center') )
      )
    ),

    # Body content
    dashboardBody(
      tabItems(
        # Load data tab
        tabItem(tabName = 'load_data',
                h2('Load your data'),
                load_dataUI('dat1')),
        tabItem(tabName = 'align',
                h2('Spectra alignment'),
                alignTestUI('test1'))
        )
      )
    )

  server <- function(input, output, session) {
    data <- load_dataServer('dat1')
    alignTestServer('test1', data)
  }

  shinyApp(ui, server)
}


