#' Differential Expression UI (tab) — supports Comparison Builder
#' @param id module id
#' @export
mod_de_tab_ui <- function(id) {
  ns <- shiny::NS(id)
  
  shiny::tabPanel(
    "Differential Expression",
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        shiny::checkboxInput(
          ns("use_cmp"),
          "Use Comparison Builder (if available)",
          value = TRUE
        ),
        
        shiny::uiOutput(ns("cmp_summary")),
        
        shiny::conditionalPanel(
          condition = sprintf("!input['%s']", ns("use_cmp")),
          ns = ns,
          shiny::selectInput(ns("metadata_column"), "Variable to test:", choices = NULL),
          shiny::selectInput(ns("reference_condition"), "Reference condition:", choices = NULL),
          shiny::selectInput(ns("test_condition"), "Test condition:", choices = NULL)
        ),
        
        shiny::hr(),
        
        shiny::sliderInput(
          ns("lfc_threshold"),
          "Log2 fold change threshold:",
          min = 0, max = 4, value = 1, step = 0.02
        ),
        
        shiny::sliderInput(
          ns("padj_threshold"),
          "Adjusted p-value threshold:",
          min = 0, max = 0.5, value = 0.05, step = 0.001
        ),
        
        shiny::actionButton(ns("run_de"), "Run Differential Expression"),
        
        shiny::hr(),
        
        shiny::numericInput(
          ns("num_genes"),
          "Number of top genes:",
          value = 100, min = 2, step = 1
        ),
        shiny::checkboxInput(
          ns("transpose_heatmap"),
          "Reverse heatmap orientation (samples in rows, genes in columns; only for 100 genes or fewer)",
          value = FALSE
        ),
        shiny::checkboxInput(
          ns("cluster_columns"),
          "Cluster columns",
          value = TRUE
        ),
        
        shiny::hr(),
        
        shiny::uiOutput(ns("ann_cols_ui")),
        
        shiny::hr(),
        
        shiny::downloadButton(ns("download_de_table"), "Download DE Table"),
        shiny::downloadButton(ns("download_heatmap"), "Download Heatmap (PDF)"),
        shiny::downloadButton(ns("download_heatmap_matrix"), "Download Heatmap Matrix (CSV)")
      ),
      
      shiny::mainPanel(
        shiny::div(
          class = "tpr-subtabs",
          shiny::tabsetPanel(
            shiny::tabPanel(
              "Results Table",
              DT::DTOutput(ns("deTable"))
            ),
            shiny::tabPanel(
              "Heatmap",
              shiny::plotOutput(ns("heatmapPlot"), height = "820px")
            )
          )
        )
      )
    )
  )
}