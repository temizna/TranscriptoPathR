# --- Module UI: Non-overlap Genes Pathway Analysis ---------------------------
#' Non-overlap Pathway tabPanel (UI)
#' @param id module id
#' @importFrom shiny NS tabPanel sidebarLayout sidebarPanel mainPanel
#' @importFrom shiny actionButton selectInput sliderInput numericInput
#' @importFrom shiny downloadButton plotOutput br
#' @importFrom DT DTOutput
#' @export
mod_nonoverlap_tab_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tabPanel(
    "Non-overlap Genes Pathway Analysis",
    sidebarLayout(
      sidebarPanel(
        shiny::selectInput(
          ns("pathway_db_nonOL"), "Pathway database:",
          choices = c("GO", "KEGG", "Reactome"), selected = "GO"
        ),
        shiny::sliderInput(ns("padj_threshold_nonOL"),
                           "Adjusted P-Value Threshold:",
                           min = 0, max = 0.5, value = 0.05, step = 0.01),
        shiny::sliderInput(ns("pathway_qval_nonOL"),
                           "Pathway Q-value threshold:",
                           min = 0, max = 0.5, value = 0.10, step = 0.01),
        shiny::numericInput(ns("showCategory_nonOL"),
                            "Show top N categories:", value = 10, min = 3, max = 50),
        shiny::actionButton(ns("run_non_overlap_pathway"),
                            "Run Non Overlap Pathway Analysis"),
        
        shiny::downloadButton(ns("download_nonOL_dot_plot"),     "Download Dot Plot"),
        shiny::downloadButton(ns("download_nonOL_heatmap_plot"), "Download Heatmap Plot"),
        shiny::downloadButton(ns("download_nonOL_tree_plot"),    "Download Tree Plot"),
        shiny::downloadButton(ns("download_nonOL_pathway_table"),"Download Pathway Table")
      ),
      mainPanel(
        shiny::tabsetPanel(
          shiny::tabPanel("Dotplot",  shiny::plotOutput(ns("dotPlot_nonOL"))),
          shiny::tabPanel("Heatmap",  shiny::plotOutput(ns("heatmapPlot_nonOL"))),
          shiny::tabPanel("Treeplot", shiny::plotOutput(ns("treePlot_nonOL"))),
          shiny::tabPanel("Table",    DT::DTOutput(ns("pathwayTable_nonOL")))
        )
      )
    )
  )
}
