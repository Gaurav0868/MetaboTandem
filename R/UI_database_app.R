# TODO finish database app
database_appUI <- dashboardPage(
  dashboardHeader(),
  dashboardSidebar(),
  dashboardBody(
    shinyWidgets::actionBttn(
      inputId = 'goHome_database',
      label = 'Return Home',
      style = 'material-flat',
      color = 'success'
    )
  )
)
