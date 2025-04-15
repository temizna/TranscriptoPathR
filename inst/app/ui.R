ui <- fluidPage(
  theme = shinythemes::shinytheme("flatly"),
  titlePanel("RNAnalyzeR: An R Shiny App to Anaylze RNASeq data"),
  tabsetPanel(
    tabPanel("Upload Data",
             sidebarLayout(
               sidebarPanel(
                 radioButtons("data_mode", "Data Mode", choices = c("Upload Your Own" = "upload", "Use Demo Data" = "demo"), selected = "upload"),
                 conditionalPanel(
                   condition = "input.data_mode == 'upload'",
                   fileInput("counts_file", "Upload RNA-Seq Counts File (CSV or XLSX)", accept = c(".csv", ".xlsx")),
                   fileInput("design_file", "Upload Design/Metadata File (CSV or XLSX)", accept = c(".csv", ".xlsx")),
                   # Inside your "Upload Data" tabPanel
                   textInput("geo_accession", "Enter GEO Accession (e.g., GSE12345)", value = ""),
                   actionButton("load_geo", "Load Data from GEO"),
                   helpText("Loads raw counts and sample metadata from GEO if available. Only studies with count data are supported."),
                   textInput("design_formula", "Design Formula", value = "~ condition"),
                   selectInput("species_select", "Select Species", 
                               choices = c("Homo sapiens", "Mus musculus"),
                               selected = "Homo sapiens"),
                   downloadButton("download_counts_template", "Download Counts Template"),
                   downloadButton("download_design_template", "Download Design Template"),
                   actionButton("load_data", "Load Data")
                 )
               ),
               mainPanel(
                 helpText("Accepted formats: .csv or .xlsx. First column = gene/sample names. Samples in columns for counts, rows for design. Design column names (e.g., 'condition') must be valid in formula."),
                 DT::DTOutput("uploaded_counts"),
                 DT::DTOutput("uploaded_design")
               )
             )
    ),
    tabPanel("Sample Select",
             sidebarLayout(
               sidebarPanel(
                 actionButton("select_all", "Select All Samples"),
                 actionButton("deselect_all", "Deselect All Samples"),
                 selectInput("sample_select", "Select Samples", choices = NULL, multiple = TRUE)
               ),
               mainPanel(DT::DTOutput("sampleTable"))
             )
    ),
    tabPanel("Individual Gene Expression", 
             sidebarLayout(
               sidebarPanel(
                 textInput("gene_select", "Enter Gene(s) (space-separated):", value = ""),
                 selectInput("group_select_geneexpr", "Group by:", choices = NULL),
                 downloadButton("download_gene_plot", "Download Plot")
               ),
               mainPanel(
                 plotOutput("geneExpressionPlot")
               )
             )
    ),
    tabPanel("Quality Check", 
             sidebarLayout(
               sidebarPanel(
                 selectInput("qc_plot_type", "Select QC Plot:", choices = c("PCA", "Sample Distance", "Mean-Variance", "Variance Histogram")),
                 selectInput("group_select_qc", "Group by:", choices = NULL),
                 checkboxInput("show_labels", "Show Labels", value = TRUE),
                 textInput("qc_plot_filename", "QC Plot Filename:", value = "qc_plot.pdf"),
                 downloadButton("download_qc_plot", "Download Plot")
               ),
               mainPanel(
                 plotOutput("qcPlot")
               )
             )
    ),  
    tabPanel("Differential Expression",
             sidebarLayout(
               sidebarPanel(
                 selectInput("metadata_column", "Metadata Column:", choices = NULL),
                 selectInput("reference_condition", "Reference Condition:", choices = NULL),
                 selectInput("test_condition", "Test Condition:", choices = NULL),
                 sliderInput("lfc_threshold", "Log2 Fold Change Threshold:", min = 0, max = 4, value = 1, step = 0.25, ticks = TRUE),
                 sliderInput("padj_threshold", "Adjusted P-Value Threshold:", min = 0, max = 0.5, value = 0.05, step = 0.01, ticks = TRUE),
                 actionButton("run_de", "Run Differential Expression")
               ),
               mainPanel(
                 DT::DTOutput("deTable"),
                 br(),
                 downloadButton("download_de_table", "Download DE Table")
               )
             )
    ),
    tabPanel("Heatmap", 
             sidebarLayout(
               sidebarPanel(
                 numericInput("num_genes", "Number of Top Genes:", value = 100, min = 10, max = 500, step = 10),
                 checkboxInput("cluster_columns", "Cluster Columns", value = TRUE),
		 actionButton("generate_heatmap", "Generate Heatmap")
               ),
               mainPanel(
                 plotOutput("heatmapPlot"),
                 br(),
                 downloadButton("download_heatmap", "Save as PDF")
               )
             )
    ),
    tabPanel("Volcano Plot",
             sidebarLayout(
               sidebarPanel(
                 textInput("volcano_select", "Enter Gene(s) (space-separated):", value = ""),
                 numericInput("volcano_gene_label", "Top N Genes to Label:", value = 10, min = 5, max = 50, step = 1),
                 sliderInput("volcano_lfc", "Log2 Fold Change Cutoff", min = 0, max = 4, value = 1, step = 0.1),
                 sliderInput("volcano_padj", "Adjusted P-value Cutoff", min = 0, max = 0.1, value = 0.05, step = 0.005),
                 downloadButton("download_volcano_plot", "Download Volcano Plot"),
                 downloadButton("download_ma_plot", "Download MA Plot")
               ),
               mainPanel(
                 plotOutput("volcanoPlot"),
                 br(),
                 plotOutput("maPlot")               
               )
             )
    ),         
    tabPanel("Cross Plot",
             sidebarLayout(
               sidebarPanel(
                 textInput("crossplot_gene_label", "Enter Gene(s) (space-separated):", value = ""),
                 selectInput("metadata_column_x", "X-axis Metadata Column:", choices = NULL),
                 selectInput("reference_condition_x", "X-axis Reference Condition:", choices = NULL),
                 selectInput("test_condition_x", "X-axis Test Condition:", choices = NULL),
                 selectInput("metadata_column_y", "Y-axis Metadata Column:", choices = NULL),
                 selectInput("reference_condition_y", "Y-axis Reference Condition:", choices = NULL),
                 selectInput("test_condition_y", "Y-axis Test Condition:", choices = NULL),
                 numericInput("crossplot_gene_count", "Top N Genes to Plot:", value = 500, min = 10, max = 5000, step = 10),
                 numericInput("crossplot_topgenes", "Top N Genes to Label:", value = 10, min = 1, max = 100, step = 1),
                 actionButton("run_crossplot", "Run Cross Plot"),
                 downloadButton("download_cross_plot", "Download Cross Plot")
               ),
               mainPanel(
                 plotOutput("crossPlot"),
                 br(),
                 plotOutput("crossVennPlot"),
                 br(),
                 plotOutput("crosspathplot")
               )
             )
    ),
    tabPanel("GSEA",
             sidebarLayout(
               sidebarPanel(
                 checkboxInput("gsea_split_dotplot", "Split Dot Plot by Activation State", value = TRUE),
                 #selectInput("gsea_metadata_column", "Metadata Column:", choices = NULL),
                 #selectInput("gsea_reference_condition", "Reference Condition:", choices = NULL),
                 #selectInput("gsea_test_condition", "Test Condition:", choices = NULL),
                 selectInput("gsea_color_scale", "Dot Plot Color By:", choices = c("p.adjust", "pvalue","qvalue"), selected = "pvalue"),
                 selectInput("gsea_db", "Select Database:", choices = c("GO", "KEGG", "Reactome", "Hallmark")),
                 numericInput("gsea_top_n", "Top N Pathways to Show in GSEA Table:", value = 10, min = 1, max = 50),
                 sliderInput("lfc_threshold", "Log2 Fold Change Threshold:", min = 0, max = 8, value = 1, step = 0.25, ticks = TRUE),
                 sliderInput("padj_threshold", "Adjusted P-Value Threshold:", min = 0, max = 0.5, value = 0.05, step = 0.01, ticks = TRUE),
                 sliderInput("gsea_pvalue", "GSEA Q-value Cutoff", min = 0, max = 1, value = 0.20, step = 0.01),
                 downloadButton("download_gsea_table", "Download GSEA Table"),
                 downloadButton("download_gsea_dot_plot", "Download GSEA Dot Plot"),
                 downloadButton("download_gsea_enrichment_plot", "Download GSEA Enrichment Plot"),
                 selectInput("gsea_selected_pathway", "Select Pathway for Enrichment Plot:", choices = NULL),
                 downloadButton("download_gsea_enrichment_plot", "Download Enrichment Plot"),
                 downloadButton("download_gsea_upset_plot", "Download Upset Plot"),
                 actionButton("run_gsea", "Run GSEA")
               ),
               mainPanel(
                 plotOutput("gseaDotPlot"),
                 br(),
                 br(),
                 plotOutput("gseaEnrichmentPlot"),
                 br(),
                 plotOutput("GSEAupsetPlot"),
                 br(),
                 DT::DTOutput("gseaTable")
               )
             )
    ),
    tabPanel("Pathway Analysis", 
             sidebarLayout(
               sidebarPanel(
                 selectInput("pathway_db", "Select Pathway Database:", choices = c("GO", "KEGG", "Reactome")),
                 selectInput("pathway_direction", "Direction:", choices = c("Up", "Down", "Both")),
                 selectInput("circular_layout", "Circular Plot Layout:", choices = c("circle", "kk", "mds"), selected = "circle"),
                 sliderInput("lfc_threshold", "Log2 Fold Change Threshold:", min = 0, max = 4, value = 1, step = 0.25, ticks = TRUE),
                 sliderInput("padj_threshold", "Adjusted P-Value Threshold:", min = 0, max = 0.5, value = 0.05, step = 0.01, ticks = TRUE),
                 sliderInput("pathway.qval", "Pathway Q-value:", min = 0, max = 0.5, value = 0.1, step = 0.01, ticks = TRUE), 
                 numericInput("max_genes", "Max Genes For Pathway Analysis:", value = 1000, min = 100, max = 1500, step = 100),
                 actionButton("run_pathway", "Run Pathway Analysis"),
                 downloadButton("download_dot_plot", "Download Dot Plot"),
                 downloadButton("download_emap_plot", "Download Emap Plot"),
                 downloadButton("download_cnet_plot", "Download Cnet Plot"),
                 downloadButton("download_pathway_table", "Download Pathway Table"),
                 downloadButton("download_circular_plot", "Download Circular Plot")
               ),
               mainPanel(
                 plotOutput("dotPlot"),
                 br(),
                 plotOutput("emapPlot"),
                 br(),
                 plotOutput("cnetPlot"),
                 br(),
                 plotOutput("circularPlot"), 
                 br(),
                 DT::DTOutput("pathwayTable")
               )
             )
    ),
    tabPanel("Pathway Plots",
            # imageOutput("PathwayImage"),
            sidebarLayout(
              sidebarPanel(
                downloadButton("download_heatmap_plot", "Download Heatmap Plot"),
                downloadButton("download_tree_plot", "Download Tree Plot"),
                downloadButton("download_upset_plot", "Download Upset Plot")
              ),
              mainPanel(
             #br(),
             #downloadButton("download_kegg_png", "Download KEGG Pathway PNG"),
             plotOutput("heatmapPlot"),
             br(),
             plotOutput("treePlot"),
             br(),
             plotOutput("upsetPlot")
            )
        )
    ),
    tabPanel("Controls",  # This is the new Controls tab
    	sidebarLayout(
     	   sidebarPanel(
      		actionButton("clear_data", "Clear Data"),  # Clear all data
                actionButton("clear_dds", "Clear DESeq2 Data"),  # Clear DESeq2 data
                actionButton("clear_pathway", "Clear Pathway Data")  # Clear Pathway data
           ),
        mainPanel(
                helpText("Use this tab to clear specific data and reset the app."),
                br(),
                textOutput("clear_status")  # Text to show clear status
        )
      )
    )
  )
)
