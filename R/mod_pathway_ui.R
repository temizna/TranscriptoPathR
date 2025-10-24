#' Pathway Analysis UI
#' @param id module id
#' @export
mod_pathway_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    sidebarLayout(
      sidebarPanel(
        shiny::checkboxInput(ns("pathway_use_custom_file"), "Use custom DE results file", value = FALSE),
        shiny::fileInput(ns("pathway_custom_file"), "Upload DE results (TSV or CSV)",
                         accept = c(".txt", ".tsv", ".csv")),
        shiny::selectInput(ns("pathway_db"), "Select Database:",
                           choices = c("GO", "KEGG", "Reactome")),
        shiny::selectInput(ns("pathway_direction"), "Direction:",
                           choices = c("All", "Up", "Down"), selected = "All"),
        shiny::numericInput(ns("lfc_threshold"), "Log2 Fold Change Threshold:", 1, 0, 8, 0.25),
        shiny::numericInput(ns("padj_threshold"), "Adjusted p-value Threshold:", 0.05, 0, 1, 0.01),
        shiny::numericInput(ns("pathway.qval"), "Q-value cutoff:", 0.2, 0, 1, 0.01),
        shiny::numericInput(ns("max_genes"), "Max genes to use:", 1000, 10, 10000, 100),
        shiny::selectInput(ns("circular_layout"), "Circular Plot Layout:",
                           choices = c("kk", "fr", "circle", "lgl", "dh", "gem")),
        shiny::actionButton(ns("run_pathway"), "Run Pathway Analysis"),
        shiny::hr(),
        shiny::downloadButton(ns("download_dot_plot"), "Download Dot Plot"),
        shiny::downloadButton(ns("download_emap_plot"), "Download Enrichment Map"),
        shiny::downloadButton(ns("download_cnet_plot"), "Download Cnet Plot"),
        shiny::downloadButton(ns("download_circular_plot"), "Download Circular Plot"),
        shiny::downloadButton(ns("download_pathway_table"), "Download Pathway Table")
      ),
      mainPanel(
        shiny::tabsetPanel(
          shiny::tabPanel("Dot Plot", shiny::plotOutput(ns("dotPlot"))),
          shiny::tabPanel("Enrichment Map", shiny::plotOutput(ns("emapPlot"))),
          shiny::tabPanel("Cnet Plot", shiny::plotOutput(ns("cnetPlot"))),
          shiny::tabPanel("Circular Plot", shiny::plotOutput(ns("circularPlot"))),
          shiny::tabPanel("Results Table", DT::DTOutput(ns("pathwayTable")))
        )
      )
    )
  )
}
