#' Run the Shiny RNA-seq App
#'
#' @export
run_app <- function() {
options(shiny.maxRequestSize = 250 * 1024^2)

  shiny::runApp(system.file("app", package = "TranscriptoPathR"))
}

