

headerbox_factory <- function(title, status, content, width = 6){
  box(title = title,
      status = status,
      solidHeader = TRUE,
      collapsible = FALSE,
      closable = FALSE,
      width = width,
      content
      )
}
