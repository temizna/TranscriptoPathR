# =========================
# Logger Module - SERVER
# =========================

#' Logger Module Server
#'
#' @param id Shiny module id
#' @param input_all The full Shiny input object (pass the top-level `input`)
#' @return reactiveVal character() with log entries
#' @importFrom shiny moduleServer reactiveVal observeEvent renderText NS downloadHandler
#' @export
mod_logger_server <- function(id, input_all) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    user_log <- reactiveVal(character())
    
    append_log <- function(title, params = NULL) {
      ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      line <- paste0("[", ts, "] ", title)
      if (!is.null(params) && length(params)) {
        kv <- vapply(params, function(p) {
          val <- input_all[[p]]
          val <- if (is.null(val)) NA_character_ else
            if (length(val) > 1) paste(val, collapse = ", ") else as.character(val)
          paste0(p, ": ", val)
        }, character(1))
        line <- paste0(line, "\nParameters: ", paste(kv, collapse = "; "))
      }
      user_log(c(user_log(), line))
    }
    
    # -----------------------
    # Load Your Data (global)
    # -----------------------
    observeEvent(input_all$load_data, {
      append_log("Clicked load_data", c("data_mode","counts_file","design_file","design_formula","species"))
    })
    observeEvent(input_all$load_geo, {
      append_log("Clicked load_geo", c("geo_accession","species"))
    })
    
    # -----------------------
    # Sample Select (module id: sample_select)
    # -----------------------
    observeEvent(input_all[["sample_select-select_all"]], {
      append_log("Clicked sample_select-select_all")
    })
    observeEvent(input_all[["sample_select-run_filter"]], {
      append_log("Clicked sample_select-run_filter")
    })
    observeEvent(input_all[["sample_select-deselect_all"]], {
      append_log("Clicked sample_select-deselect_all")
    })
    
    # -----------------------
    # Comparison Builder (module id: cmp)
    # -----------------------
    observeEvent(input_all[["cmp-apply"]], {
      mode <- input_all[["cmp-mode"]]
      base_params <- c("cmp-group_var","cmp-mode","cmp-covariates")
      extra <- switch(
        mode,
        "simple" = c("cmp-level_test","cmp-level_ref"),
        "pooled" = c("cmp-levels_test","cmp-levels_ref","cmp-keep_levels_only"),
        "paired" = c("cmp-subject_col","cmp-level_test","cmp-level_ref"),
        character(0)
      )
      append_log("Clicked cmp-apply", c(base_params, extra))
    })
    
    # -----------------------
    # Differential Expression (module id: de)
    # -----------------------
    observeEvent(input_all[["de-run_de"]], {
      append_log("Clicked de-run_de", c(
        "de-metadata_column", "de-reference_condition", "de-test_condition",
        "de-lfc_threshold", "de-padj_threshold", "de-num_genes", "de-cluster_columns"
      ))
    })
    # Back-compat:
    observeEvent(input_all$run_de, {
      append_log("Clicked run_de (global)", c(
        "metadata_column","reference_condition","test_condition",
        "lfc_threshold","padj_threshold","num_genes","cluster_columns"
      ))
    })
    
    # -----------------------
    # Cross Plot (module id: cross)
    # -----------------------
    observeEvent(input_all[["cross-run_crossplot"]], {
      append_log("Clicked cross-run_crossplot", c(
        "cross-crossplot_gene_label",
        "cross-metadata_column_x","cross-reference_condition_x","cross-test_condition_x",
        "cross-metadata_column_y","cross-reference_condition_y","cross-test_condition_y",
        "cross-crossplot_gene_count","cross-crossplot_topgenes","cross-cross_enrich_method"
      ))
    })
    # Cross downloads
    observeEvent(input_all[["cross-download_cross_plot"]],                 { append_log("Clicked cross-download_cross_plot") })
    observeEvent(input_all[["cross-download_cross_venn_plot"]],            { append_log("Clicked cross-download_cross_venn_plot") })
    observeEvent(input_all[["cross-download_overlap_genes"]],              { append_log("Clicked cross-download_overlap_genes") })
    observeEvent(input_all[["cross-download_cross_category_heatmap"]],     { append_log("Clicked cross-download_cross_category_heatmap") })
    observeEvent(input_all[["cross-download_cross_pathway_plot"]],         { append_log("Clicked cross-download_cross_pathway_plot") })
    # Back-compat:
    observeEvent(input_all$run_crossplot, {
      append_log("Clicked run_crossplot (global)", c(
        "crossplot_gene_label",
        "metadata_column_x","reference_condition_x","test_condition_x",
        "metadata_column_y","reference_condition_y","test_condition_y",
        "crossplot_gene_count","crossplot_topgenes","cross_enrich_method"
      ))
    })
    
    # -----------------------
    # Volcano (global)
    # -----------------------
    observeEvent(input_all$download_volcano_plot, { append_log("Clicked download_volcano_plot") })
    observeEvent(input_all$download_ma_plot,      { append_log("Clicked download_ma_plot") })
    
    # -----------------------
    # Quality Check (global)
    # -----------------------
    observeEvent(input_all$download_qc_plot, {
      append_log("Clicked download_qc_plot", c("qc_plot_type","group_select_qc","show_labels","qc_plot_filename"))
    })
    
    # -----------------------
    # Genes of Interest Heatmap (global)
    # -----------------------
    observeEvent(input_all$download_goi_heatmap, {
      append_log("Clicked download_goi_heatmap", c("goi_group_column","cluster_columns","goi_file"))
    })
    
    # -----------------------
    # GSEA (module id: gsea)
    # -----------------------
    observeEvent(input_all[["gsea-run_gsea"]], {
      append_log("Clicked gsea-run_gsea", c(
        "gsea-gsea_split_dotplot","gsea-gsea_color_scale","gsea-gsea_db",
        "gsea-gsea_top_n","gsea-lfc_threshold","gsea-padj_threshold","gsea-gsea_pvalue"
      ))
    })
    observeEvent(input_all[["gsea-download_gsea_dot_plot"]],         { append_log("Clicked gsea-download_gsea_dot_plot") })
    observeEvent(input_all[["gsea-download_gsea_enrichment_plot"]],  { append_log("Clicked gsea-download_gsea_enrichment_plot") })
    observeEvent(input_all[["gsea-download_gsea_upset_plot"]],       { append_log("Clicked gsea-download_gsea_upset_plot") })
    observeEvent(input_all[["gsea-download_gsea_table"]],            { append_log("Clicked gsea-download_gsea_table") })
    # Back-compat:
    observeEvent(input_all$run_gsea, {
      append_log("Clicked run_gsea (global)", c(
        "gsea_split_dotplot","gsea_color_scale","gsea_db",
        "gsea_top_n","lfc_threshold","padj_threshold","gsea_pvalue"
      ))
    })
    
    # -----------------------
    # GSVA / ssGSEA (module id: gsva)
    # -----------------------
    observeEvent(input_all[["gsva-run"]], {
      append_log("Clicked gsva-run", c(
        "gsva-gset_source","gsva-gsva_db","gsva-gsva_method","gsva-mx_opt",
        "gsva-min_gs_size","gsva-max_gs_size","gsva-use_de_defaults",
        "gsva-gsva_metadata_column","gsva-gsva_reference_condition","gsva-gsva_test_condition"
      ))
    })
    observeEvent(input_all[["gsva-dl_table"]],    { append_log("Clicked gsva-dl_table") })
    observeEvent(input_all[["gsva-dl_scores"]],   { append_log("Clicked gsva-dl_scores") })
    observeEvent(input_all[["gsva-dl_volcano"]],  { append_log("Clicked gsva-dl_volcano") })
    observeEvent(input_all[["gsva-dl_heatmap"]],  { append_log("Clicked gsva-dl_heatmap") })
    observeEvent(input_all[["gsva-dl_boxplots"]], { append_log("Clicked gsva-dl_boxplots") })
    
    # -----------------------
    # Pathway Analysis (module id: pathway)
    # -----------------------
    observeEvent(input_all[["pathway-run_pathway"]], {
      append_log("Clicked pathway-run_pathway", c(
        "pathway-pathway_db","pathway-pathway_direction","pathway-circular_layout",
        "pathway-lfc_threshold","pathway-padj_threshold","pathway-pathway.qval","pathway-max_genes"
      ))
    })
    observeEvent(input_all[["pathway-download_dot_plot"]],        { append_log("Clicked pathway-download_dot_plot") })
    observeEvent(input_all[["pathway-download_emap_plot"]],       { append_log("Clicked pathway-download_emap_plot") })
    observeEvent(input_all[["pathway-download_cnet_plot"]],       { append_log("Clicked pathway-download_cnet_plot") })
    observeEvent(input_all[["pathway-download_circular_plot"]],   { append_log("Clicked pathway-download_circular_plot") })
    observeEvent(input_all[["pathway-download_pathway_table"]],   { append_log("Clicked pathway-download_pathway_table") })
    # Back-compat:
    observeEvent(input_all$run_pathway, {
      append_log("Clicked run_pathway (global)", c(
        "pathway_db","pathway_direction","circular_layout",
        "lfc_threshold","padj_threshold","pathway.qval","max_genes"
      ))
    })
    
    # -----------------------
    # Pathway Plots (module id: pathplots) downloads
    # -----------------------
    observeEvent(input_all[["pathplots-download_pathheatmap_plot"]], { append_log("Clicked pathplots-download_pathheatmap_plot") })
    observeEvent(input_all[["pathplots-download_tree_plot"]],        { append_log("Clicked pathplots-download_tree_plot") })
    observeEvent(input_all[["pathplots-download_upset_plot"]],       { append_log("Clicked pathplots-download_upset_plot") })
    # Back-compat:
    observeEvent(input_all$download_pathheatmap_plot, { append_log("Clicked download_pathheatmap_plot (global)") })
    observeEvent(input_all$download_tree_plot,        { append_log("Clicked download_tree_plot (global)") })
    observeEvent(input_all$download_upset_plot,       { append_log("Clicked download_upset_plot (global)") })
    
    # -----------------------
    # Non-overlap Pathway (module id: nonOL)
    # -----------------------
    observeEvent(input_all[["nonOL-run_non_overlap_pathway"]], {
      append_log("Clicked nonOL-run_non_overlap_pathway", c(
        "nonOL-pathway_db_nonOL","nonOL-padj_threshold_nonOL",
        "nonOL-pathway_qval_nonOL","nonOL-showCategory_nonOL"
      ))
    }, ignoreInit = TRUE)
    observeEvent(input_all[["nonOL-download_nonOL_dot_plot"]],       { append_log("Clicked nonOL-download_nonOL_dot_plot") })
    observeEvent(input_all[["nonOL-download_nonOL_heatmap_plot"]],   { append_log("Clicked nonOL-download_nonOL_heatmap_plot") })
    observeEvent(input_all[["nonOL-download_nonOL_tree_plot"]],      { append_log("Clicked nonOL-download_nonOL_tree_plot") })
    observeEvent(input_all[["nonOL-download_pathway_nonOL_table"]],  { append_log("Clicked nonOL-download_pathway_nonOL_table") })
    
    # -----------------------
    # TF Enrichment (module id: tf)
    # -----------------------
    observeEvent(input_all[["tf-run_tf_enrichment"]], {
      append_log("Clicked tf-run_tf_enrichment", c(
        "tf-tf_data_source","tf-enrichment_method","tf-gene_direction",
        "tf-lfc_threshold","tf-padj_threshold","tf-tf.qval"
      ))
    })
    observeEvent(input_all[["tf-download_tf_dotplot"]],        { append_log("Clicked tf-download_tf_dotplot") })
    observeEvent(input_all[["tf-download_tf_ridgeplot"]],      { append_log("Clicked tf-download_tf_ridgeplot") })
    observeEvent(input_all[["tf-download_tf_results_table"]],  { append_log("Clicked tf-download_tf_results_table") })
    
    # -----------------------
    # PCA Cluster (module id: pca)
    # -----------------------
    observeEvent(input_all[["pca-run_pca"]], {
      append_log("Clicked pca-run_pca", c(
        "pca-max_pc","pca-variance_threshold","pca-percent_of_data",
        "pca-similarity_threshold","pca-max_genes","pca-min_genes","pca-pca_enrich_method"
      ))
    })
    observeEvent(input_all[["pca-download_pca_loadings"]],               { append_log("Clicked pca-download_pca_loadings") })
    observeEvent(input_all[["pca-download_pca_enrichment"]],             { append_log("Clicked pca-download_pca_enrichment") })
    observeEvent(input_all[["pca-download_reconstructed_heatmap"]],      { append_log("Clicked pca-download_reconstructed_heatmap") })
    observeEvent(input_all[["pca-download_sample_correlation_heatmap"]], { append_log("Clicked pca-download_sample_correlation_heatmap") })
    observeEvent(input_all[["pca-download_pca_loadings_table"]],         { append_log("Clicked pca-download_pca_loadings_table") })
    observeEvent(input_all[["pca-download_pca_enrichment_table"]],       { append_log("Clicked pca-download_pca_enrichment_table") })
    # Legacy global PCA
    observeEvent(input_all$run_pca, {
      append_log("Clicked run_pca (global)", c(
        "max_pc","variance_threshold","percent_of_data",
        "similarity_threshold","max_genes","min_genes","pca_enrich_method"
      ))
    })
    
    # -----------------------
    # Cancer Gene Census (global)
    # -----------------------
    observeEvent(input_all$run_cancer_gene_census, { append_log("Clicked run_cancer_gene_census") })
    observeEvent(input_all$download_cancer_gene_table, { append_log("Clicked download_cancer_gene_table") })
    observeEvent(input_all$download_cgc_venn_plot,     { append_log("Clicked download_cgc_venn_plot") })
    
    # Output and persistence
    output$log_text <- renderText({ paste(user_log(), collapse = "\n\n") })
    output$download_log <- downloadHandler(
      filename = function() paste0("user_session_log_", Sys.Date(), ".txt"),
      content  = function(file) writeLines(user_log(), con = file)
    )
    
    session$onSessionEnded(function() {
      if (nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_"))) return(invisible(NULL))
      log_dir <- file.path(getwd(), "logs")
      if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
      file_path <- file.path(log_dir, paste0("user_session_log_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), ".txt"))
      log_snapshot <- isolate(user_log())
      writeLines(log_snapshot, file_path)
    })
    
    return(user_log)
  })
}

# =========================
# Logger Module - UI
# =========================

#' Logger Module UI
#'
#' @param id Shiny module id
#' @return UI for the log and a download button
#' @importFrom shiny NS tagList h3 div verbatimTextOutput downloadButton
#' @export
mod_logger_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h3("User Session Log"),
    div(
      style = "overflow-y: auto; height: 500px; border: 1px solid #ccc; padding: 10px; white-space: pre-wrap; background-color: #f9f9f9; font-family: monospace; font-size: 12px;",
      verbatimTextOutput(ns("log_text"))
    ),
    downloadButton(ns("download_log"), "Download Log")
  )
}
