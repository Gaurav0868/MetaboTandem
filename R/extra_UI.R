sidebarUI <- function(id){
  ns <- NS(id)
  tagList(
    selectInput(ns('group'), label = 'Grouping Variables',
                choices = c('treatment'), selected = 'treatment'),
    hr(),
    colourpicker::colourInput(ns('col1'), label = 'CTR'),
    colourpicker::colourInput(ns('col2'), label = 'WP')
  )
}
