# app/ui.R

bs <- getOption("TranscriptoPathR.bs", default = tpr_theme())

ui <- shiny::fluidPage(
  theme = bs,
  div(class = "tpr-app",
  shiny::titlePanel("TranscriptoPathR: An R Shiny App to Anaylze RNASeq data"),
  shiny::tabsetPanel(
    
    # ---- Load Your Data ----
    shiny::tabPanel(
      "Load Your Data",
      shiny::sidebarLayout(
        shiny::sidebarPanel(
          shiny::radioButtons(
            "data_mode", "Data Mode",
            choices = c("Load Your Own" = "load", "Use Demo Data" = "demo"),
            selected = "load"
          ),
          shiny::conditionalPanel(
            condition = "input.data_mode == 'load'",
            shiny::fileInput("counts_file", "Load RNA-Seq Counts File (CSV or XLSX)",
                             accept = c(".csv", ".xlsx")),
            shiny::fileInput("design_file", "Load Design/Metadata File (CSV or XLSX)",
                             accept = c(".csv", ".xlsx")),
            shiny::textInput(
              "design_formula",
              "Study Design Formula (you can put main column name. If more than one use + between names",
              value = "~ condition"
            ),
            shiny::selectInput(
              "species", "Select Species",
              choices = c("Homo sapiens", "Mus musculus", "Rattus norvegicus"),
              selected = "Homo sapiens"
            ),
            shiny::actionButton("load_data", "Load Data"),
            shiny::downloadButton("download_counts_template", "Download Counts Template"),
            shiny::downloadButton("download_design_template", "Download Sample Template")
          )
        ),
        shiny::mainPanel(
          shiny::helpText(
            "Accepted formats: .csv or .xlsx. First column = gene/sample names. ",
            "Raw count data from TCGA can also be used. Use of TPM, FPKM, RSEM is not recommended. ",
            "Samples should be in columns for the counts data, and in rows as row names for the design sheet. ",
            "Design column names (e.g., '~condition') must be valid in a formula. ",
            "If using more than one column, use + (e.g., ~condition+time). ",
            "If a chosen column does not vary across samples, DESeq2 will error."
          ),
          shiny::tags$hr(),
          shiny::tags$figure(
            shiny::tags$img(src = "sample_sheet_guide.png", class = "guide-img",
                            alt = "Valid sample/metadata sheet"),
            shiny::tags$figcaption("Sample/Design Sheet: rows = samples; columns = metadata.")
          ),
          shiny::tags$figure(
            shiny::tags$img(src = "expression_matrix_guide.png", class = "guide-img",
                            alt = "Valid counts matrix"),
            shiny::tags$figcaption("Counts Matrix: rows = genes; columns = samples. First column = gene IDs.")
          ),
          DT::DTOutput("uploaded_counts"),
          DT::DTOutput("uploaded_design")
        )
      )
    ),
    
    # ---- Sample Select ----
    shiny::tabPanel(
      "Sample Select",
      {
        ns <- shiny::NS("sample_select")
        shiny::sidebarLayout(
          shiny::sidebarPanel(
            shiny::actionButton(ns("select_all"), "Select All Samples"),
            shiny::uiOutput(ns("dynamic_filters")),
            shiny::actionButton(ns("run_filter"), "Apply Filters"),
            shiny::actionButton(ns("deselect_all"), "Deselect All Samples")
          ),
          shiny::mainPanel(DT::DTOutput(ns("sampleTable")))
        )
      }
    ),
    
    # ---- Individual Gene Expression ----
    shiny::tabPanel(
      "Individual Gene Expression",
      shiny::sidebarLayout(
        shiny::sidebarPanel(
          shiny::textInput("gene_select", "Enter Gene(s) (space-separated):", value = ""),
          shiny::selectInput("group_select_geneexpr", "Group by:", choices = NULL),
          shiny::downloadButton("download_gene_plot", "Download Plot"),
          shiny::downloadButton("download_gene_stats", "Download ANOVA Table")
        ),
        shiny::mainPanel(
          shiny::plotOutput("geneExpressionPlot"),
          shiny::br(),
          DT::dataTableOutput("geneExpressionStats")
        )
      )
    ),
    
    # ---- Quality Check ----
    shiny::tabPanel(
      "Quality Check",
      shiny::sidebarLayout(
        shiny::sidebarPanel(
          shiny::selectInput("qc_plot_type", "Select QC Plot:",
                             choices = c("PCA", "Sample Distance", "Mean-Variance", "Variance Histogram")),
          shiny::selectInput("group_select_qc", "Group by:", choices = NULL),
          shiny::checkboxInput("show_labels", "Show Labels", value = TRUE),
          shiny::textInput("qc_plot_filename", "QC Plot Filename:", value = "qc_plot.pdf"),
          shiny::downloadButton("download_qc_plot", "Download Plot")
        ),
        shiny::mainPanel(
          shiny::plotOutput("qcPlot"),
          shiny::br(),
          shiny::textOutput("qcPlotDescription")
        )
      )
    ),
    
    # ---- Genes of Interest Heatmap ----
    shiny::tabPanel(
      "Genes of Interest Heatmap",
      shiny::sidebarLayout(
        shiny::sidebarPanel(
          shiny::fileInput("goi_file", "Upload Gene List (CSV, single column)"),
          shiny::selectInput("goi_group_column", "Group annotation variable:", choices = NULL),
          shiny::checkboxInput("cluster_columns", "Cluster Columns", value = TRUE)
        ),
        shiny::mainPanel(
          shiny::plotOutput("goi_heatmap"),
          shiny::downloadButton("download_goi_heatmap", "Download Heatmap")
        )
      )
    ),
    
    # ---- Comparison Builder (returns tabPanel) ----
    mod_easy_compare_ui("cmp"),
    
    # ---- Differential Expression (returns tabPanel) ----
    mod_de_tab_ui("de"),
    
    # ---- Volcano Plot ----
    shiny::tabPanel(
      "Volcano Plot",
      shiny::sidebarLayout(
        shiny::sidebarPanel(
          shiny::textInput("volcano_select", "Enter Gene(s) (space-separated):", value = ""),
          shiny::numericInput("volcano_gene_label", "Top N Genes to Label:", value = 10, min = 5, max = 50, step = 1),
          shiny::sliderInput("volcano_lfc",  "Log2 Fold Change Cutoff", min = 0, max = 4, value = 1, step = 0.1),
          shiny::sliderInput("volcano_padj", "Adjusted P-value Cutoff",  min = 0, max = 0.1, value = 0.05, step = 0.005),
          shiny::downloadButton("download_volcano_plot", "Download Volcano Plot"),
          shiny::downloadButton("download_ma_plot",      "Download MA Plot")
        ),
        shiny::mainPanel(
          shiny::div(class = "tpr-surface", shiny::plotOutput("volcanoPlot")),
          shiny::br(),
          shiny::div(class = "tpr-surface", shiny::plotOutput("maPlot"))
        )
      )
    ),
    
    # ---- Cross Plot (returns tabPanel) ----
    mod_cross_tab_ui("cross"),
    
    # ---- GSEA (module returns plain UI -> wrap) ----
    shiny::tabPanel("GSEA", mod_gsea_tab_ui("gsea")),
    
    # ---- Pathway Analysis (module returns plain UI -> wrap) ----
    shiny::tabPanel("Pathway Analysis", mod_pathway_ui("pathway")),
    
    # ---- Pathway Plots (module returns plain UI -> wrap) ----
    shiny::tabPanel("Pathway Plots", mod_pathway_plots_ui("pathplots")),
    
    # ---- Non-overlap Genes Pathway Analysis (module returns plain UI -> wrap) ----
    shiny::tabPanel("Non-overlap Genes Pathway Analysis", mod_nonoverlap_tab_ui("nonOL")),
    
    # ---- GSVA / ssGSEA (module returns plain UI -> wrap) ----
    shiny::tabPanel("GSVA / ssGSEA", mod_gsva_ui("gsva")),
    
    # ---- Transcription Factor Enrichment (module returns plain UI -> wrap) ----
    shiny::tabPanel("Transcription Factor Enrichment", mod_tf_tab_ui("tf")),
    
    # ---- Cancer Gene Census ----
    shiny::tabPanel(
      "Cancer Gene Census",
      shiny::sidebarLayout(
        shiny::sidebarPanel(
          shiny::actionButton("run_cancer_gene_census", "Run Cancer Gene Census Analysis"),
          shiny::downloadButton("download_cancer_gene_table", "Download Overlapping Genes Table"),
          shiny::downloadButton("download_cgc_venn_plot",     "Download Overlapping Venn Plot")
        ),
        shiny::mainPanel(
          shiny::plotOutput("vennPlot"),
          shiny::br(),
          DT::DTOutput("overlappingGenesTable")
        )
      )
    ),
    
    # ---- PCA Cluster (returns tabPanel) ----
    mod_pca_tab_ui("pca"),
    
    # ---- Log (module returns plain UI -> wrap) ----
    shiny::tabPanel("Log", mod_logger_ui("logger"))
  )
  )
)
