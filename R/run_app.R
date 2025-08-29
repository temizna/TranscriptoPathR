#' Launch the TranscriptoPathR Shiny app
#'
#' @param ... Additional arguments passed to [shiny::runApp()]
#' @param width  Browser window width (for chrome app-mode launcher). Default 1440.
#' @param height Browser window height (for chrome app-mode launcher). Default 900.
#' @param app_mode If TRUE, launch Chrome in app-mode (chromeless window). Default TRUE.
#'
#' @return Runs the app; returns whatever [shiny::runApp()] returns.
#' @export
#' @importFrom shiny runApp
#' @importFrom thematic thematic_shiny
run_app <- function(..., width = 1440, height = 900, app_mode = TRUE) {
  options(shiny.maxRequestSize = 250 * 1024^2)

  thematic::thematic_shiny()

  # Use app theme option if already set; otherwise fall back to your iOS-like theme.
  bs <- getOption("TranscriptoPathR.bs", default = tpr_theme())
  options(TranscriptoPathR.bs = bs)

  app_dir <- system.file("app", package = "TranscriptoPathR")
  if (!nzchar(app_dir) || !dir.exists(app_dir)) {
    stop("App directory not found in inst/app. Did you install the package with inst/app included?")
  }

  # Allow override via options(TranscriptoPathR.launch_browser = ...)
  default_launcher <-
    if (exists("chrome_launcher", mode = "function")) {
      chrome_launcher(width, height, app_mode = app_mode)
    } else {
      TRUE  # fall back to default browser if chrome_launcher isn't available
    }

  launcher <- getOption("TranscriptoPathR.launch_browser", default = default_launcher)

  shiny::runApp(app_dir, launch.browser = launcher, display.mode = "normal", ...)
}

