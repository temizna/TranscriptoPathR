# R/mod_pca_tab_ui.R

#' PCA Cluster UI
#' @param id module id
#' @importFrom shiny NS tagList sidebarLayout sidebarPanel mainPanel
#' @importFrom shiny tabsetPanel tabPanel plotOutput textOutput helpText
#' @importFrom shiny sliderInput numericInput selectInput actionButton downloadButton br
#' @export
mod_pca_tab_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tabPanel(
    "PCA Cluster",
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        shiny::sliderInput(ns("max_pc"), "Max number of PCs", min = 2, max = 20, value = 10),
        shiny::sliderInput(ns("variance_threshold"), "Variance threshold (SD)", min = 0.25, max = 4, value = 1, step = 0.05),
        shiny::sliderInput(ns("percent_of_data"), "Percent of data coverage", min = 0.25, max = 1, value = 0.8, step = 0.01),
        shiny::sliderInput(ns("similarity_threshold"), "Similarity for gene to Component assignment", min = 0, max = 1, value = 0.1),
        shiny::numericInput(ns("max_genes"), "Max genes per PC", value = 500, min = 50, max = 1000),
        shiny::numericInput(ns("min_genes"), "Min genes per PC", value = 50, min = 10, max = 500),
        shiny::selectInput(ns("pca_enrich_method"), "Enrichment Method",
                           choices = c("groupGO", "enrichGO", "enrichKEGG", "enrichDO", "enrichPathway"),
                           selected = "enrichGO"
        ),
        shiny::uiOutput(ns("ann_cols_ui")),
        shiny::actionButton(ns("run_pca"), "Run PCA Clustering"),
        shiny::br(),
        shiny::downloadButton(ns("download_pca_loadings"), "Download Contributing Genes Heatmap"),
        shiny::downloadButton(ns("download_pca_enrichment"), "Download Enrichment Plot"),
        shiny::downloadButton(ns("download_reconstructed_heatmap"), "Download Reconstructed Heatmap"),
        shiny::downloadButton(ns("download_sample_correlation_heatmap"), "Download Sample Correlation"),
        shiny::downloadButton(ns("download_pca_loadings_table"), "Download Contributing Genes Table"),
        shiny::downloadButton(ns("download_pca_enrichment_table"), "Download Enrichment Table")
      ),
      shiny::mainPanel(
        shiny::tabsetPanel(
          shiny::tabPanel(
            "Overview & Selection",
            shiny::helpText(
              "This module performs Principal Component Analysis (PCA) to help you assess ",
              "sample structure and genes contributing to components."
            ),
            shiny::plotOutput(ns("pca_variance_plot")),
            shiny::br(),
            shiny::textOutput(ns("selected_pc_text"))
          ),
          shiny::tabPanel(
            "Contributing Genes Heatmap",
            shiny::plotOutput(ns("pca_loadings_heatmap"), height = "700px")
          ),
          shiny::tabPanel(
            "Sample Correlation Heatmap",
            shiny::plotOutput(ns("pca_sample_heatmap"), height = "600px")
          ),
          shiny::tabPanel(
            "Pathway Enrichment",
            shiny::plotOutput(ns("pca_enrichment_plot"), height = "650px")
          )
        )
      )
    )
  )
}
