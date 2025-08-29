# R/mod_easy_compare_ui.R

#' Easy Comparison Builder UI
#'
#' Build Test vs Reference comparisons without contrasts.
#'
#' @param id module id
#' @importFrom shiny NS tabPanel sidebarLayout sidebarPanel mainPanel
#' @importFrom shiny selectInput radioButtons actionButton uiOutput tableOutput helpText br
#' @export
mod_easy_compare_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tabPanel(
    "Complex Comparison Builder",
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        shiny::helpText("Build Test vs Reference comparisons without writing contrasts."),
        shiny::selectInput(ns("group_var"), "Grouping column", choices = NULL),
        shiny::radioButtons(
          ns("mode"), "Comparison type",
          choices = c(
            "One level vs one level" = "simple",
            "Pool levels vs pool levels" = "pooled",
            "Paired (subject ID required)" = "paired"
          ),
          selected = "simple"
        ),
        shiny::uiOutput(ns("mode_controls")),
        shiny::uiOutput(ns("covariate_ui")),
        shiny::br(),
        shiny::actionButton(ns("apply"), "Apply comparison")
      ),
      shiny::mainPanel(
        shiny::uiOutput(ns("preview_counts_ui")),
        shiny::br(),
        shiny::tableOutput(ns("preview_samples")),
        shiny::br(),
        shiny::tableOutput(ns("summary"))
      )
    )
  )
}
