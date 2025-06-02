ui <- fluidPage(
  theme = shinythemes::shinytheme("flatly"),
  titlePanel("TranscriptoPathR: An R Shiny App to Anaylze RNASeq data"),
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
                   textInput("design_formula", "Study Design Formula (you can put main column name. If more than one use + between names", value = "~ condition"),
                   selectInput("species", "Select Species", 
                               choices = c("Homo sapiens", "Mus musculus"),
                               selected = "Homo sapiens"),
                   actionButton("load_data", "Load Data"),
                   downloadButton("download_counts_template", "Download Counts Template"),
                   downloadButton("download_design_template", "Download Sample Template"),
                   helpText("Loads raw counts and sample metadata from GEO if available. Only studies with count data are supported. Currently under construction"),
                  # textInput("geo_accession", "Enter GEO Accession (e.g., GSE287953)", value = ""),
                  # actionButton("load_geo", "Load Data from GEO")
                 )
               ),
               mainPanel(
                 helpText("Accepted formats: .csv or .xlsx. First column = gene/sample names. 
                         Raw count data from TCGA can also be used for this analysis. Use of TPM, FPKM, RSEM is not recommended.
                         Samples  should be in columns for the counts data, and in rows as row names for design/sample/experiment sheet. 
                         Design column names (e.g., '~condition') must be valid in formula.
                         If using more than one column name, use + e.g. ~condition+time. 
                         If multiple column names are used and they do not vary accors other conditions and samples DESEQ will return an error!
                         To fix this error remove the column name that does not vary."),
                 DT::DTOutput("uploaded_counts"),
                 DT::DTOutput("uploaded_design")
               )
             )
    ),
    tabPanel("Sample Select",
             { 
               ns <- NS("sample_select")   # <<-- add this inside the panel code block
               
               sidebarLayout(
                 sidebarPanel(
                   actionButton(ns("select_all"), "Select All Samples"),
                   uiOutput(ns("dynamic_filters")),
                   actionButton(ns("run_filter"), "Apply Filters"),
                   actionButton(ns("deselect_all"), "Deselect All Samples")
                 ),
                 mainPanel(
                   DT::DTOutput(ns("sampleTable"))
                 )
               )
             }
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
    tabPanel(
      title = "Genes of Interest Heatmap",
      sidebarLayout(
        sidebarPanel(
          fileInput("goi_file", "Upload Gene List (CSV, single column)"),
          checkboxInput("cluster_columns", "Cluster Columns", value = TRUE)
        ),
        mainPanel(
          plotOutput("goi_heatmap"),
          downloadButton("download_goi_heatmap", "Download Heatmap")
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
                 actionButton("run_de", "Run Differential Expression"),
                 numericInput("num_genes", "Number of Top Genes:", value = 100, min = 10, max = 500, step = 10),
                 checkboxInput("cluster_columns", "Cluster Columns", value = TRUE),
                 downloadButton("download_heatmap", "Save as PDF")
               ),
               mainPanel(
                 plotOutput("heatmapPlot"),
                 br(),
                 DT::DTOutput("deTable"),
                 br(),
                 downloadButton("download_de_table", "Download DE Table")
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
                 numericInput("crossplot_gene_count", "Top N Genes to Plot:", value = 2000, min = 10, max = 5000, step = 10),
                 numericInput("crossplot_topgenes", "Top N Genes to Label:", value = 10, min = 0, max = 100, step = 1),
                 selectInput("cross_enrich_method", "Enrichment Method",
                             choices = c("groupGO", "enrichGO", "enrichKEGG", "enrichPathway"),
                             selected = "enrichKEGG"),
                 actionButton("run_crossplot", "Run Cross Plot"),
                 textInput("crossplot_filename", "Filename (with extension)", value = "crossplot.pdf"),
                 textInput("crosspath_filename", "Filename (with extension)", value = "crosspath.pdf"),
                 downloadButton("download_cross_plot", "Download Cross Plot"),
                 downloadButton("download_cross_venn_plot", "Download Venn Diagram"),
                 downloadButton("download_overlap_genes", "Download Overlapping Genes"),
                 downloadButton("download_cross_category_heatmap", "Download Heatmap"),
                 downloadButton("download_cross_pathway_plot", "Download Pathway Plot")
               ),
               mainPanel(
                 plotOutput("crossPlot"),
                 br(),
                 plotOutput("crossVennPlot"),
                 br(),
                 plotOutput("crossCategoryHeatmap"),
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
                 selectInput("gsea_db", "Select Database:", choices = c("GO", "KEGG", "Reactome", "Hallmark","Cancer Cell Atlas",
                                                                        "Cancer Gene Neighbourhoods", "Cancer Modules","Txn Factor Targets")),
                 numericInput("gsea_top_n", "Top N Pathways to Show in GSEA Table:", value = 10, min = 1, max = 50),
                 sliderInput("lfc_threshold", "Log2 Fold Change Threshold:", min = 0, max = 8, value = 1, step = 0.25, ticks = TRUE),
                 sliderInput("padj_threshold", "Adjusted P-Value Threshold:", min = 0, max = 0.5, value = 0.05, step = 0.01, ticks = TRUE),
                 sliderInput("gsea_pvalue", "GSEA Q-value Cutoff", min = 0, max = 1, value = 0.20, step = 0.01),
                 downloadButton("download_gsea_table", "Download GSEA Table"),
                 downloadButton("download_gsea_dot_plot", "Download GSEA Dot Plot"),
                 selectInput("gsea_selected_pathway", "Select Pathway for Enrichment Plot:", choices = NULL),
                 downloadButton("download_gsea_enrichment_plot", "Download GSEA Enrichment Plot"),
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
                 downloadButton("download_pathheatmap_plot", "Download Heatmap Plot"),
                 downloadButton("download_tree_plot", "Download Tree Plot"),
                 downloadButton("download_upset_plot", "Download Upset Plot")
               ),
               mainPanel(
                 #br(),
                 #downloadButton("download_kegg_png", "Download KEGG Pathway PNG"),
                 plotOutput("pathheatmapPlot"),
                 br(),
                 plotOutput("treePlot"),
                 br(),
                 plotOutput("upsetPlot")
               )
             )
    ),
    tabPanel("Non-overlap Genes Pathway Analysis", 
             sidebarLayout(
               sidebarPanel(
                 actionButton("run_non_overlap_pathway", "Run Non Overlap Pathway Analysis"),
                 downloadButton("download_nonOL_dot_plot", "Download Dot Plot"),
                 downloadButton("download_nonOL_heatmap_plot", "Download Heatmap Plot"),
                 downloadButton("download_nonOL_tree_plot", "Download Tree Plot"),
                 downloadButton("download_nonOL_pathway_table", "Download Pathway Table")
               ),
               mainPanel(
                 plotOutput("dotPlot_nonOL"),
                 br(),
                 plotOutput("heatmapPlot_nonOL"),
                 br(),
                 plotOutput("treePlot_nonOL"),
                 br(),
                 DT::DTOutput("pathwayTable_nonOL")
               )
             )
    ),
    tabPanel("Transcription Factor Enrichment",
             sidebarLayout(
               sidebarPanel(
                 selectInput("tf_data_source", "Select Transcription Factor Dataset:",
                             choices = c(
                               "TRANSFAC and JASPAR PWMs",
                               "ENCODE and ChEA Consensus TFs from ChIP-X",
                               "TRRUST_Transcription_Factors_2019",
                               "TF_Perturbations_Followed_by_Expression",
                               "hTFtarget",
                               "TFLink"), selected = NULL),
                 selectInput("enrichment_method", "Select Enrichment Method:", choices = c("Over_representation","GSEA")),
                 selectInput("gene_direction", "Direction:", choices = c("Up", "Down", "Both")),
                 sliderInput("lfc_threshold", "Log2 Fold Change Threshold:", min = 0, max = 4, value = 1, step = 0.25, ticks = TRUE),
                 sliderInput("padj_threshold", "Adjusted P-Value Threshold:", min = 0, max = 0.5, value = 0.05, step = 0.01, ticks = TRUE),
                 sliderInput("tf.qval", "TF Q-value:", min = 0, max = 0.5, value = 0.1, step = 0.01, ticks = TRUE),
                 actionButton("run_tf_enrichment", "Run TF Enrichment"),
                 downloadButton("download_tf_dotplot", "Download Dot Plot"),
                 downloadButton("download_tf_ridgeplot", "Download Ridge Plot"),
                 downloadButton("download_tf_results_table", "Download Table")
               ),
               mainPanel(
                 plotOutput("tf_dotplot"),
                 br(),
                 plotOutput("tf_ridgeplot"),
                 br(),
                 DT::DTOutput("tf_results_table")
               )
             )
    ),
    tabPanel("Cancer Gene Census",
             sidebarLayout(
               sidebarPanel(
                 actionButton("run_cancer_gene_census", "Run Cancer Gene Census Analysis"),
                 downloadButton("download_cancer_gene_table", "Download Overlapping Genes Table"),
                 downloadButton("download_cgc_venn_plot", "Download Overlapping Venn Plot")
               ),
               mainPanel(
                 plotOutput("vennPlot"),
                 br(),
                 DT::DTOutput("overlappingGenesTable")
               )
             )
    ),
    tabPanel("PCA Cluster",
             sidebarLayout(
               sidebarPanel(
                 sliderInput("max_pc", "Max number of PCs", min = 2, max = 20, value = 10),
                 sliderInput("variance_threshold", "Variance threshold", min = 0.25, max = 4, value = 1, step = 0.05),
                 sliderInput("percent_of_data", "Percent of data coverage", min = 0.25, max = 1, value = 0.8, step = 0.01),
                 sliderInput("similarity_threshold", "Similarity for gene to Component assignment", min = 0, max = 1, value = 0.1),
                 numericInput("max_genes", "Max genes per PC", value = 500, min = 50, max = 1000),
                 numericInput("min_genes", "Min genes per PC", value = 50, min = 10, max = 500),
                 selectInput("pca_enrich_method", "Enrichment Method",
                             choices = c("groupGO", "enrichGO", "enrichKEGG", "enrichDO", "enrichPathway"),
                             selected = "enrichGO"),
                 actionButton("run_pca", "Run PCA Clustering"),
                 downloadButton("download_pca_loadings", "Download Contributing Genes Heatmap"),
                 downloadButton("download_pca_enrichment", "Download Enrichment Plot"),
                 downloadButton("download_reconstructed_heatmap", "Download Reconstructed Heatmap"),
                 downloadButton("download_sample_correlation_heatmap", "Download Sample Correlation"),
                 downloadButton("download_pca_loadings_table", "Download Contributing Genes Table"),
                 downloadButton("download_pca_enrichment_table", "Download Enrichment Table")
               ),
               mainPanel(
                 helpText("This module performs Principal Component Analysis (PCA) to help you determine whether advanced methods like Consensus 
                 Clustering (CC) or Non-negative Matrix Factorization (NMF) may be useful for further exploring sample substructure and 
                 identifying gene contributors to sub-clusters. Since CC and NMF are computationally intensive, they are not included in this app.
                 If needed, we recommend that users consult with an expert to run these analyses externally.

                Inputs and parameters:

                Number of PCs: The initial maximum number of principal components to compute.
                
                Variance threshold: Filters genes by standard deviation; only genes above the threshold are retained.
                
                Percent of data coverage: Sets the cumulative variance cutoff; selects the number of PCs required 
                to explain this percentage of total variance.
                
                Similarity threshold: Determines whether a gene can be associated with multiple PCs based on its relative loading and 
                if the loadings of different PCs are withing this threshold.
                
                Gene limits: You can set minimum and maximum numbers of genes per component for downstream pathway enrichment analysis.

                          "),
                 plotOutput("pca_variance_plot"),
                 br(),
                 textOutput("selected_pc_text"),
                 br(),
                 plotOutput("pca_loadings_heatmap"),
                # br(),
                # plotOutput("pca_reconstructed_heatmap"),
                 br(),
                 plotOutput("pca_sample_heatmap"),
                 br(),
                 plotOutput("pca_enrichment_plot")
               )
             )
    ),
     tabPanel(
      "Log",
      mod_logger_ui("logger")  # 'logger' is the module ID
    )
  )
)
