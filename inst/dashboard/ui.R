## ui.R
## Generates the html making the dashboard UI. Uses a custom solution made with Pug.js and Bootstrap.
## Author: Matt Piekenbrock
addResourcePath(prefix = "js", directoryPath = system.file(file.path("dashboard", "www", "js"), package = "Mapper"))
addResourcePath(prefix = "css", directoryPath = system.file(file.path("dashboard", "www", "css"), package = "Mapper"))
ui <- shinyUI(
  includeHTML(system.file(file.path("dashboard", "www", "index.html"), package = "Mapper"))
)
