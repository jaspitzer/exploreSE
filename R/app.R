build_app <- function(...) {
  ui <- .build_ui(...)
  server <- .build_server(...)
  app <- shinyApp(ui = ui, server = server)
}
