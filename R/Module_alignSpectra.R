alignTestUI <- function(id){
  ns <- NS(id)
  tagList(
    numericInput(ns('noise'), 'Noise threshold', value = 0),
    numericInput(ns('snt'), 'Signal-to-noise threshold', value = 0),
    sliderInput(ns('rt_range'), 'Range of retention time for testing [s]', value = c(0, 60), max = 1200, min = 0),
    sliderInput(ns('mz_range'), 'Range of m/z for testing', value = c(100, 200), max = 1200, min = 0),
    actionButton(ns('test'), 'Test Parameters'),
    #plotOutput(ns('test_plot'))
    verbatimTextOutput(ns('test_table'))
  )
}

alignTestServer <- function(id, data){
  moduleServer(id, function(input, output, session){

    # data_cent <- eventReactive(input$test,{
    #   data()$data
    # })
    # output$test_plot <- renderPlot({
    #   test_peak_picking(data_cent(),
    #                     p.width = c(0,100),
    #                     mz.range = c(input$mz_range[1], input$mz_range[2]),
    #                     rt.range = c(input$rt_range[1], input$rt_range[2]),
    #                     snt = input$snt,
    #                     noise = input$noise)
    # })

    data_cent <- NULL

    observe({
      req(input$test)

      data_cent <- data$data()
    })

    output$test_table <- renderPrint(
      print(data_cent)
    )
  })

}

