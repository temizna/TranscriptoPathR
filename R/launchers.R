#' Launch Chrome at a given window size
#' @param width,height window size in pixels
#' @param app_mode use Chrome's app mode (no tabs)
#' @export
chrome_launcher <- function(width = 1440, height = 900, app_mode = TRUE) {
  function(url) {
    args <- c("--new-window", sprintf("--window-size=%d,%d", width, height))
    if (app_mode) args <- c(args, paste0("--app=", shQuote(url))) else args <- c(args, shQuote(url))
    
    # Try common Chrome paths cross-platform
    paths <- c(
      Sys.which("google-chrome"),
      Sys.which("chromium"),
      Sys.which("chrome"),
      "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
      "C:/Program Files/Google/Chrome/Application/chrome.exe",
      "C:/Program Files (x86)/Google/Chrome/Application/chrome.exe"
    )
    paths <- paths[nzchar(paths)]
    if (length(paths) && file.exists(paths[1])) {
      system2(paths[1], args, wait = FALSE)
    } else {
      # Fallback: whatever the system default browser is
      utils::browseURL(url)
    }
  }
}
