# R/mod_pathway_plots_ui.R
#' Pathway Plots UI (tabbed: heatmap / tree / upset)
#' @param id module id
#' @export
mod_pathway_plots_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    sidebarLayout(
      sidebarPanel(
        shiny::downloadButton(ns("download_pathheatmap_plot"), "Download Heatmap Plot"),
        shiny::downloadButton(ns("download_tree_plot"),        "Download Tree Plot"),
        shiny::downloadButton(ns("download_upset_plot"),       "Download Upset Plot")
      ),
      mainPanel(
        shiny::tabsetPanel(
          shiny::tabPanel("Heatmap", shiny::plotOutput(ns("pathheatmapPlot"))),
          shiny::tabPanel("Tree",    shiny::plotOutput(ns("treePlot"))),
          shiny::tabPanel("Upset",   shiny::plotOutput(ns("upsetPlot")))
        )
      )
    )
  )
}
