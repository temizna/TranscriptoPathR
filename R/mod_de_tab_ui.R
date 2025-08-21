#' Differential Expression tabPanel (UI module)
#' @param id module id
#' @export
mod_de_tab_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tabPanel(
    "Differential Expression",
    sidebarLayout(
      sidebarPanel(
        shiny::selectInput(ns("metadata_column"), "Variable to test:", choices = NULL),
        shiny::selectInput(ns("reference_condition"), "Reference Condition:", choices = NULL),
        shiny::selectInput(ns("test_condition"), "Test Condition:", choices = NULL),
        shiny::sliderInput(ns("lfc_threshold"), "Log2 Fold Change Threshold:",
                           min = 0, max = 4, value = 1, step = 0.25, ticks = TRUE),
        shiny::sliderInput(ns("padj_threshold"), "Adjusted P-Value Threshold:",
                           min = 0, max = 0.5, value = 0.05, step = 0.01, ticks = TRUE),
        shiny::actionButton(ns("run_de"), "Run Differential Expression"),
        shiny::numericInput(ns("num_genes"), "Number of Top Genes:",
                            value = 100, min = 10, max = 500, step = 10),
        shiny::checkboxInput(ns("cluster_columns"), "Cluster Columns", value = TRUE),
        shiny::downloadButton(ns("download_heatmap"), "Save as PDF")
      ),
      mainPanel(
        shiny::plotOutput(ns("heatmapPlot")),
        shiny::br(),
        DT::DTOutput(ns("deTable")),
        shiny::br(),
        shiny::downloadButton(ns("download_de_table"), "Download DE Table")
      )
    )
  )
}
