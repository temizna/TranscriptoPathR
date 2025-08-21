# R/mod_pathway_ui.R
#' Pathway Analysis UI (tabbed)
#' @param id module id
#' @export
mod_pathway_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    sidebarLayout(
      sidebarPanel(
        shiny::selectInput(ns("pathway_db"), "Select Pathway Database:",
                           choices = c("GO", "KEGG", "Reactome")),
        shiny::selectInput(ns("pathway_direction"), "Direction:", choices = c("Up", "Down", "Both")),
        shiny::selectInput(ns("circular_layout"), "Circular Plot Layout:",
                           choices = c("circle", "kk", "mds"), selected = "circle"),
        shiny::sliderInput(ns("lfc_threshold"),  "Log2 Fold Change Threshold:", min = 0, max = 4, value = 1, step = 0.25, ticks = TRUE),
        shiny::sliderInput(ns("padj_threshold"), "Adjusted P-Value Threshold:", min = 0, max = 0.5, value = 0.05, step = 0.01, ticks = TRUE),
        shiny::sliderInput(ns("pathway.qval"),   "Pathway Q-value threshold for enrichment:", min = 0, max = 0.5, value = 0.1, step = 0.01, ticks = TRUE),
        shiny::numericInput(ns("max_genes"), "Max Genes For Pathway Analysis:", value = 1000, min = 100, max = 1500, step = 100),
        shiny::actionButton(ns("run_pathway"), "Run Pathway Analysis"),
        shiny::hr(),
        shiny::downloadButton(ns("download_dot_plot"),      "Download Dot Plot"),
        shiny::downloadButton(ns("download_emap_plot"),     "Download Emap Plot"),
        shiny::downloadButton(ns("download_cnet_plot"),     "Download Cnet Plot"),
        shiny::downloadButton(ns("download_circular_plot"), "Download Circular Plot"),
        shiny::downloadButton(ns("download_pathway_table"), "Download Pathway Table")
      ),
      mainPanel(
        shiny::tabsetPanel(
          shiny::tabPanel("Dot Plot",      shiny::plotOutput(ns("dotPlot"))),
          shiny::tabPanel("Emap Plot",     shiny::plotOutput(ns("emapPlot"))),
          shiny::tabPanel("Cnet Plot",     shiny::plotOutput(ns("cnetPlot"))),
          shiny::tabPanel("Circular Plot", shiny::plotOutput(ns("circularPlot"))),
          shiny::tabPanel("Results Table", DT::DTOutput(ns("pathwayTable")))
        )
      )
    )
  )
}
