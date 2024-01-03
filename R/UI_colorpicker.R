#' ColorPickerUI
#'
#' UI for color picker
#'
#' @return UI elements for colorpicker widget
#'
#' @param ... UI elements to be passed to MetaboTandem app
#' @param id id for Namespacing
#'
#'
#' @export
colorPickerUI <- function(id){
  ns <- NS(id)
  tagList(
    shinyWidgets::radioGroupButtons(ns('set_colors'),
                                    'Color palette',
                                    choices = c('Default',
                                                'User defined'),
                                    individual = TRUE,
                                    checkIcon = list(
                                      yes = tags$i(class = "fa fa-circle",
                                                   style = "color: steelblue"),
                                      no = tags$i(class = "fa fa-circle-o",
                                                  style = "color: steelblue"))),

    uiOutput(ns('color_picker')),

    textOutput(ns('message'))
  )
}


colorPickerServer <- function(id, metadata){
  moduleServer(id, function(input, output, session){

    ns <- NS(id)

    # Get columns of the metadata for selectizeInput

    output$color_picker <- renderUI({
      if(input$set_colors == 'User defined'){
        list(
          varSelectizeInput(ns('color_group'),
                            label = 'Choose grouping variables for plot',
                            data = metadata()),

          uiOutput(ns('pickers'))
        )
      }
    })

    treatments <- reactive({
      paste0('color_', unique(dplyr::pull(metadata(), input$color_group)))
    })

    output$pickers <- renderUI({

      purrr::map(treatments(), function(x){
        label <- stringr::str_remove(x, 'color_')
        colourpicker::colourInput(ns(x), label = label)
      })
    }) %>%
      bindEvent(treatments())

    user_colors <- reactive({
      if(input$set_colors == 'User defined'){
        colors <- purrr::map_chr(treatments(), ~purrr::`%||%`(input[[.x]], ""))
        treatments <- stringr::str_remove(treatments(), 'color_')
        names(colors) <- treatments
        colors
      } else {
        'Default'
      }
    })

    output$message <- renderText({
      if(all(user_colors() == 'Default')){
        'Default colors selected'
      } else {
        paste('User colors selected based on:', input$group)
      }
    })

    return(user_colors)

  })


}

