# === Module: mod_cross_tab_ui ===
#' Cross Plot tabPanel (UI module)
#' @param id module id
#' @importFrom shiny NS tabPanel sidebarLayout sidebarPanel mainPanel
#' @importFrom shiny tabsetPanel plotOutput downloadButton
#' @importFrom shiny selectInput sliderInput numericInput actionButton textInput br
#' @export
mod_cross_tab_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tabPanel(
    "Cross Plot",
    sidebarLayout(
      sidebarPanel(
        shiny::textInput(ns("crossplot_gene_label"), "Enter Gene(s) (space-separated):", value = ""),
        shiny::selectInput(ns("metadata_column_x"), "X-axis Variable to test:", choices = NULL),
        shiny::selectInput(ns("reference_condition_x"), "X-axis Reference Condition:", choices = NULL),
        shiny::selectInput(ns("test_condition_x"), "X-axis Test Condition:", choices = NULL),
        
        shiny::selectInput(ns("metadata_column_y"), "Y-axis Variable to test:", choices = NULL),
        shiny::selectInput(ns("reference_condition_y"), "Y-axis Reference Condition:", choices = NULL),
        shiny::selectInput(ns("test_condition_y"), "Y-axis Test Condition:", choices = NULL),
        
        shiny::numericInput(ns("crossplot_gene_count"), "Top N Genes to Plot:",
                            value = 2000, min = 10, max = 5000, step = 10),
        shiny::numericInput(ns("crossplot_topgenes"), "Top N Genes to Label:",
                            value = 10, min = 0, max = 100, step = 1),
        
        shiny::selectInput(ns("cross_enrich_method"), "Enrichment Method",
                           choices = c("groupGO", "enrichGO", "enrichKEGG", "enrichPathway"),
                           selected = "enrichKEGG"),
        
        shiny::actionButton(ns("run_crossplot"), "Run Cross Plot"),
        shiny::br(), shiny::br(),
        
        shiny::textInput(ns("crossplot_filename"), "Cross Plot file name", value = "crossplot.pdf"),
        shiny::textInput(ns("crosspath_filename"),  "Pathway Dotplot file name", value = "crosspath.pdf"),
        
        shiny::downloadButton(ns("download_cross_plot"),              "Download Cross Plot"),
        shiny::downloadButton(ns("download_cross_venn_plot"),         "Download Venn Diagram"),
        shiny::downloadButton(ns("download_overlap_genes"),           "Download Overlapping Genes"),
        shiny::downloadButton(ns("download_cross_category_heatmap"),  "Download Heatmap"),
        shiny::downloadButton(ns("download_cross_pathway_plot"),      "Download Pathway Plot")
      ),
      mainPanel(
        shiny::tabsetPanel(
          shiny::tabPanel("Cross Plot",         shiny::plotOutput(ns("crossPlot"))),
          shiny::tabPanel("Venn",               shiny::plotOutput(ns("crossVennPlot"))),
          shiny::tabPanel("Heatmap",            shiny::plotOutput(ns("crossCategoryHeatmap"))),
          shiny::tabPanel("Pathway Dotplot",    shiny::plotOutput(ns("crosspathplot")))
        )
      )
    )
  )
}
