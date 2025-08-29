#' Differential Expression UI (tab) — supports Comparison Builder
#' @param id module id
#' @export
mod_de_tab_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tabPanel(
    "Differential Expression",
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        # --- NEW: Use comparison builder ---
        shiny::checkboxInput(ns("use_cmp"),
                             "Use Comparison Builder (if available)",
                             value = TRUE),
        shiny::uiOutput(ns("cmp_summary")),   # shows current cmp label + covariates
        
        # Manual controls shown only if NOT using cmp
        shiny::conditionalPanel(
          condition = sprintf("!input['%s']", ns("use_cmp")),
          shiny::selectInput(ns("metadata_column"), "Variable to test:", choices = NULL),
          shiny::selectInput(ns("reference_condition"), "Reference Condition:", choices = NULL),
          shiny::selectInput(ns("test_condition"), "Test Condition:", choices = NULL)
        ),
        
        shiny::hr(),
        shiny::sliderInput(ns("lfc_threshold"), "Log2 Fold Change Threshold:",
                           min = 0, max = 4, value = 1, step = 0.1),
        shiny::sliderInput(ns("padj_threshold"), "Adjusted P-Value Threshold:",
                           min = 0, max = 0.5, value = 0.05, step = 0.005),
        
        shiny::actionButton(ns("run_de"), "Run Differential Expression"),
        shiny::hr(),
        shiny::numericInput(ns("num_genes"), "Number of Top Genes:", value = 100, min = 2, step = 1),
        shiny::checkboxInput(ns("cluster_columns"), "Cluster Columns", value = TRUE),
        shiny::hr(),
        # NEW: Additional column annotations (like ssGSEA)
        shiny::uiOutput(ns("ann_cols_ui")),
        shiny::downloadButton(ns("download_de_table"), "Download DE Table"),
        shiny::downloadButton(ns("download_heatmap"), "Download Heatmap (PDF)")
      ),
      shiny::mainPanel(
        shiny::div(class = "tpr-subtabs",
                   shiny::tabsetPanel(
                     shiny::tabPanel("Results Table", DT::DTOutput(ns("deTable"))),
                     shiny::tabPanel("Heatmap", shiny::plotOutput(ns("heatmapPlot"), height = "820px"))
                   )
        )
      )
    )
  )
}
