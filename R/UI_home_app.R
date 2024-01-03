home_ui <- dashboardPage(
  #md = TRUE,
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
    ),
    leftUi = tagList(
      h4('Welcome to MetaboTandem', style = "color:white")
    )
  ),

  # Sidebar content
  dashboardSidebar(),

  # Body content
  dashboardBody(
    fluidRow(
      box(title = 'Estimate Parameters',
          headerBorder = FALSE,
          background = 'purple',
          collapsible = FALSE,
          closable = FALSE,
          gradient = TRUE,
          footer = p('Use Autotuner to estimate the best parameters for your run',
                     style = 'color:gray'),
          fluidRow(
            column(6, align = 'center', offset = 3,
                   shinyWidgets::actionBttn(inputId = 'goAutotuner',
                                            icon = icon('lightbulb'),
                                            style = 'jelly',
                                            size = 'lg'))
          )
      ),
      box(title = 'Full Data analysis',
          headerBorder = FALSE,
          background = 'light-blue',
          collapsible = FALSE,
          closable = FALSE,
          gradient = TRUE,
          footer = p('Run the full analysis starting from mas spectra data',
                     style = 'color:gray'),
          fluidRow(
            column(6, align = 'center', offset = 3,
                   shinyWidgets::actionBttn(inputId = 'goMain',
                                            icon = icon('chart-simple'),
                                            style = 'jelly',
                                            size = 'lg'))
          )
      ),
      box(title = 'Manage Databases',
          headerBorder = FALSE,
          background = 'olive',
          collapsible = FALSE,
          closable = FALSE,
          gradient = TRUE,
          footer = p('Create custom database or download public databases',
                     style = 'color:gray'),
          fluidRow(
            column(6, align = 'center', offset = 3,
                   shinyWidgets::actionBttn(inputId = 'goDatabase',
                                            icon = icon('database'),
                                            color = 'default',
                                            style = 'jelly',
                                            size = 'lg'))
          )
      ),
      box(title = 'Annotate MGF',
          headerBorder = FALSE,
          background = 'maroon',
          collapsible = FALSE,
          closable = FALSE,
          gradient = TRUE,
          footer = p('Annotate MS2 spectra using custom or public databases',
                     style = 'color:gray'),
          fluidRow(
            column(6, align = 'center', offset = 3,
                   shinyWidgets::actionBttn(inputId = 'goAnnotate',
                                            icon = icon('magnifying-glass'),
                                            style = 'jelly',
                                            size = 'lg'))
          )
      )
    )
  )
)


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

test_appUI <- dashboardPage(
  dashboardHeader(),
  dashboardSidebar(),
  dashboardBody(
    shinyWidgets::actionBttn(
      inputId = 'goHome_autotuner',
      label = 'Return Home',
      style = 'material-flat',
      color = 'success'
    )
  )
)

