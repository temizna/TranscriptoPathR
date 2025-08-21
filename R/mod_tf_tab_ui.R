#' Transcription Factor Enrichment UI (content-only)
#' @param id module id
#' @importFrom shiny NS sidebarLayout sidebarPanel mainPanel
#' @importFrom shiny selectInput sliderInput numericInput actionButton
#' @importFrom shiny downloadButton plotOutput tabsetPanel tabPanel br
#' @importFrom DT DTOutput
#' @export
mod_tf_tab_ui <- function(id) {
  ns <- shiny::NS(id)
  sidebarLayout(
    sidebarPanel(
      shiny::selectInput(
        ns("tf_data_source"), "Select Transcription Factor Dataset:",
        choices = c(
          "TRANSFAC and JASPAR PWMs",
          "ENCODE and ChEA Consensus TFs from ChIP-X",
          "TRRUST_Transcription_Factors_2019",
          "TF_Perturbations_Followed_by_Expression",
          "hTFtarget",
          "TFLink"
        ),
        selected = NULL
      ),
      shiny::selectInput(
        ns("enrichment_method"), "Select Enrichment Method:",
        choices = c("Over_representation", "GSEA")
      ),
      shiny::selectInput(
        ns("gene_direction"), "Direction:",
        choices = c("Up", "Down", "Both")
      ),
      shiny::sliderInput(
        ns("lfc_threshold"), "Log2 Fold Change Threshold:",
        min = 0, max = 4, value = 1, step = 0.25, ticks = TRUE
      ),
      shiny::sliderInput(
        ns("padj_threshold"), "Adjusted P-Value Threshold:",
        min = 0, max = 0.5, value = 0.05, step = 0.01, ticks = TRUE
      ),
      shiny::sliderInput(
        ns("tf.qval"), "TF Q-value threshold for enrichment:",
        min = 0, max = 0.5, value = 0.1, step = 0.01, ticks = TRUE
      ),
      shiny::actionButton(ns("run_tf_enrichment"), "Run TF Enrichment"),
      shiny::downloadButton(ns("download_tf_dotplot"),     "Download Dot Plot"),
      shiny::downloadButton(ns("download_tf_ridgeplot"),   "Download Ridge Plot"),
      shiny::downloadButton(ns("download_tf_results_table"), "Download Table")
    ),
    mainPanel(
      shiny::tabsetPanel(
        shiny::tabPanel("Dotplot",   shiny::plotOutput(ns("tf_dotplot"))),
        shiny::tabPanel("Ridgeplot", shiny::plotOutput(ns("tf_ridgeplot"))),
        shiny::tabPanel("Table",     DT::DTOutput(ns("tf_results_table")))
      )
    )
  )
}
