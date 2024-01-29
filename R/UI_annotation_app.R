#' Solo annotation UI
#'`
#' @param id character used to specify namespace, see [`shiny::NS`][shiny::NS()]
#'
#' @return a \code{shiny::\link[shiny]{tagList}} containing UI elements
#'

annotationappUI <- function(id){
  dashboardPage(
    app_header(),
    dashboardSidebar(
      fluidRow(
        column(6, align = 'center', offset = 3,
               shinyWidgets::actionBttn(inputId = 'goHome_annot',
                                        icon = icon('house'),
                                        style = 'material-circle',
                                        color = 'success',
                                        size = 'sm'))
      )
    ),
    dashboardBody(
      dbAnnotationUI(id)
  )
)
}
