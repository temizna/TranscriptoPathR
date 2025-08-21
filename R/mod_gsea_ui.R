# R/mod_gsea_tab_ui.R
#' GSEA UI (tabbed)
#' @param id module id
#' @export
mod_gsea_tab_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    sidebarLayout(
      sidebarPanel(
        shiny::checkboxInput(ns("gsea_split_dotplot"), "Split Dot Plot by Activation State", value = TRUE),
        shiny::selectInput(ns("gsea_color_scale"), "Dot Plot Color By:",
                           choices = c("p.adjust", "pvalue", "qvalue", "NES"), selected = "pvalue"),
        shiny::selectInput(ns("gsea_db"), "Select Database:",
                           choices = c("GO","KEGG","Reactome","Hallmark","Cancer Cell Atlas",
                                       "Cancer Gene Neighbourhoods","Cancer Modules","Txn Factor Targets")),
        shiny::numericInput(ns("gsea_top_n"), "Top N Pathways to Show in GSEA Table:", value = 10, min = 1, max = 50),
        shiny::sliderInput(ns("lfc_threshold"),  "Log2 Fold Change Threshold:", min = 0, max = 8, value = 1, step = 0.25, ticks = TRUE),
        shiny::sliderInput(ns("padj_threshold"), "Adjusted P-Value Threshold:", min = 0, max = 0.5, value = 0.05, step = 0.01, ticks = TRUE),
        shiny::sliderInput(ns("gsea_pvalue"), "GSEA Q-value threshold for enrichment", min = 0, max = 1, value = 0.20, step = 0.01),
        shiny::selectInput(ns("gsea_selected_pathway"), "Select Pathway for Enrichment Plot:", choices = NULL),
        shiny::actionButton(ns("run_gsea"), "Run GSEA"),
        shiny::hr(),
        shiny::downloadButton(ns("download_gsea_table"),            "Download GSEA Table"),
        shiny::downloadButton(ns("download_gsea_dot_plot"),         "Download GSEA Dot Plot"),
        shiny::downloadButton(ns("download_gsea_enrichment_plot"),  "Download GSEA Enrichment Plot"),
        shiny::downloadButton(ns("download_gsea_upset_plot"),       "Download Upset Plot")
      ),
      mainPanel(
        shiny::tabsetPanel(
          shiny::tabPanel("Dot Plot",        shiny::plotOutput(ns("gseaDotPlot"))),
          shiny::tabPanel("Enrichment Plot", shiny::plotOutput(ns("gseaEnrichmentPlot"))),
          shiny::tabPanel("Upset Plot",      shiny::plotOutput(ns("GSEAupsetPlot"))),
          shiny::tabPanel("Results Table",   DT::DTOutput(ns("gseaTable")))
        )
      )
    )
  )
}
