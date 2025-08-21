# =========================
# Logger Module - SERVER
# =========================

#' Logger Module Server
#'
#' @param id Shiny module id
#' @param input_all The full Shiny input object (pass the top-level `input`)
#' @return reactiveVal character() with log entries
#' @importFrom shiny moduleServer reactiveVal observeEvent renderText NS
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
    # Differential Expression (module id: de)
    # -----------------------
    observeEvent(input_all[["de-run_de"]], {
      append_log("Clicked de-run_de", c(
        "de-metadata_column", "de-reference_condition", "de-test_condition",
        "de-lfc_threshold", "de-padj_threshold", "de-num_genes", "de-cluster_columns"
      ))
    })
    # Back-compat if still non-modular:
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
    # GSEA (module id: gsea)
    # -----------------------
    observeEvent(input_all[["gsea-run_gsea"]], {
      append_log("Clicked gsea-run_gsea", c(
        "gsea-gsea_split_dotplot","gsea-gsea_color_scale","gsea-gsea_db",
        "gsea-gsea_top_n","gsea-lfc_threshold","gsea-padj_threshold","gsea-gsea_pvalue"
      ))
    })
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
        "gsva-min_gs_size","gsva-max_gs_size",
        "gsva-use_de_defaults",
        # if user overrides DE selections:
        "gsva-gsva_metadata_column","gsva-gsva_reference_condition","gsva-gsva_test_condition"
      ))
    })
    # Optional: log GSVA downloads
    observeEvent(input_all[["gsva-dl_table"]],   { append_log("Clicked gsva-dl_table") })
    observeEvent(input_all[["gsva-dl_scores"]],  { append_log("Clicked gsva-dl_scores") })
    observeEvent(input_all[["gsva-dl_volcano"]], { append_log("Clicked gsva-dl_volcano") })
    observeEvent(input_all[["gsva-dl_heatmap"]], { append_log("Clicked gsva-dl_heatmap") })
    observeEvent(input_all[["gsva-dl_boxplots"]],{ append_log("Clicked gsva-dl_boxplots") })
    
    # -----------------------
    # Pathway Analysis (module id: pathway)
    # -----------------------
    observeEvent(input_all[["pathway-run_pathway"]], {
      append_log("Clicked pathway-run_pathway", c(
        "pathway-pathway_db","pathway-pathway_direction","pathway-circular_layout",
        "pathway-lfc_threshold","pathway-padj_threshold","pathway-pathway.qval","pathway-max_genes"
      ))
    })
    # Back-compat:
    observeEvent(input_all$run_pathway, {
      append_log("Clicked run_pathway (global)", c(
        "pathway_db","pathway_direction","circular_layout",
        "lfc_threshold","padj_threshold","pathway.qval","max_genes"
      ))
    })
    
    # -----------------------
    # Pathway Plots (module id: pathplots) - downloads only
    # -----------------------
    observeEvent(input_all[["pathplots-download_pathheatmap_plot"]], {
      append_log("Clicked pathplots-download_pathheatmap_plot")
    })
    observeEvent(input_all[["pathplots-download_tree_plot"]], {
      append_log("Clicked pathplots-download_tree_plot")
    })
    observeEvent(input_all[["pathplots-download_upset_plot"]], {
      append_log("Clicked pathplots-download_upset_plot")
    })
    # Back-compat:
    observeEvent(input_all$download_pathheatmap_plot, { append_log("Clicked download_pathheatmap_plot (global)") })
    observeEvent(input_all$download_tree_plot,        { append_log("Clicked download_tree_plot (global)") })
    observeEvent(input_all$download_upset_plot,       { append_log("Clicked download_upset_plot (global)") })
    
    # -----------------------
    # 
    # -----------------------
    # Non-overlap Pathway (namespaced)
    observeEvent(input_all[["nonOL-run_non_overlap_pathway"]], {
      timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      params <- c("nonOL-pathway_db_nonOL", "nonOL-padj_threshold_nonOL", "nonOL-pathway_qval_nonOL", "nonOL-showCategory_nonOL")
      entry <- paste0("[", timestamp, "] User clicked nonOL-run_non_overlap_pathway\nParameters: ",
                      paste(sapply(params, function(pn) paste0(pn, ": ", as.character(input_all[[pn]]))), collapse="; "))
      user_log(c(user_log(), entry))
    }, ignoreInit = TRUE)
    
    
    # -----------------------
    # TF Enrichment (namespaced)
    # -----------------------
    observeEvent(input_all[["tf-run_tf_enrichment"]], {
      timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      params <- c("tf-tf_data_source","tf-enrichment_method","tf-gene_direction",
                  "tf-lfc_threshold","tf-padj_threshold","tf-tf.qval")
      entry <- paste0("[", timestamp, "] User clicked tf-run_tf_enrichment\nParameters: ",
                      paste(sapply(params, function(pn) paste0(pn, ": ", as.character(input_all[[pn]]))), collapse="; "))
      user_log(c(user_log(), entry))
    })
    
    
    # -----------------------
    # PCA Cluster (global)
    # -----------------------
    observeEvent(input_all$run_pca, {
      append_log("Clicked run_pca", c(
        "max_pc","variance_threshold","percent_of_data",
        "similarity_threshold","max_genes","min_genes","pca_enrich_method"
      ))
    })
    
    # Outputs
    output$log_text <- renderText({ paste(user_log(), collapse = "\n\n") })
    
    output$download_log <- downloadHandler(
      filename = function() paste0("user_session_log_", Sys.Date(), ".txt"),
      content  = function(file) writeLines(user_log(), con = file)
    )
    
    # Save on session end unless running R CMD check
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
