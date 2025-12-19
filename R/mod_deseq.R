# === Module: mod_de_server (comparison-safe, with extra annotations) ===
#' Differential Expression Server Module (DESeq2)
#'
#' If a comparison builder `cmp` is supplied and has a valid selection,
#' the DE tab runs on that subset with design `~ cmp_group` and a
#' Reference/Test contrast. It NEVER mutates filtered_data_rv$samples.
#'
#' Also supports extra column annotations for the heatmap via input$ann_cols
#' (multi-select UI added in the DE sidebar).
#'
#' @param id module id (must match mod_de_tab_ui)
#' @param filtered_data_rv reactiveValues with $samples, $counts, $norm_counts, $species
#' @param filtered_dds_rv reactiveVal (optional) holding a DESeqDataSet
#' @param res_reactive reactiveVal to store DESeq2 results data.frame
#' @param cmp optional list from mod_easy_compare_server(): included_samples(), cmp_factor(), tag(), label()
#' @return list reactives: group_var, ref_level, test_level, cmp_factor, included_samples
#' @import DESeq2
#' @importFrom DT renderDT datatable
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation draw
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grid gpar
#' @importFrom grDevices hcl.colors
#' @importFrom stats as.formula
#' @importFrom shiny req showNotification downloadHandler renderPlot renderUI observe observeEvent updateSelectInput
#' @export
mod_de_server <- function(id, filtered_data_rv, filtered_dds_rv, res_reactive, cmp = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # ---- helpers -------------------------------------------------------------
    .build_top_annotation <- function(samples, primary_col, extra_cols, col_order) {
      # primary group (fallback if missing)
      if (!primary_col %in% colnames(samples)) {
        primary_col <- colnames(samples)[1]
      }
      ann_df <- data.frame(
        Group = factor(samples[col_order, primary_col, drop = TRUE]),
        check.names = FALSE
      )
      
      # add extra tracks (if UI provides ann_cols)
      extra_cols <- intersect(extra_cols %||% character(0), setdiff(colnames(samples), primary_col))
      for (cn in extra_cols) {
        x <- samples[col_order, cn, drop = TRUE]
        if (is.numeric(x) && length(unique(x)) > 12) {
          q <- cut(x,
                   breaks = unique(stats::quantile(x, probs = seq(0, 1, 0.25), na.rm = TRUE)),
                   include.lowest = TRUE, dig.lab = 6)
          ann_df[[cn]] <- q
        } else {
          ann_df[[cn]] <- factor(as.character(x))
        }
      }
      
      # colors
      col_list <- list()
      glv <- levels(ann_df$Group)
      if (length(glv) <= 8) {
        col_list$Group <- setNames(RColorBrewer::brewer.pal(max(3, length(glv)), "Set2")[seq_along(glv)], glv)
      } else {
        col_list$Group <- setNames(grDevices::hcl.colors(length(glv), "Dark 3"), glv)
      }
      if (length(extra_cols)) {
        for (cn in extra_cols) {
          lv <- levels(ann_df[[cn]])
          pal <- if (length(lv) <= 8) {
            sets <- c("Set3", "Pastel1", "Dark2", "Accent")
            pnm  <- sets[(match(cn, extra_cols) - 1) %% length(sets) + 1]
            base <- try(RColorBrewer::brewer.pal(max(3, length(lv)), pnm), silent = TRUE)
            if (inherits(base, "try-error")) grDevices::hcl.colors(length(lv), "Dark 3") else base
          } else grDevices::hcl.colors(length(lv), "Dark 3")
          col_list[[cn]] <- setNames(pal[seq_along(lv)], lv)
        }
      }
      
      ComplexHeatmap::HeatmapAnnotation(df = ann_df, col = col_list)
    }
    
    `%||%` <- function(a, b) if (!is.null(a)) a else b
    
    # ---- populate variable to test ------------------------------------------
    shiny::observe({
      shiny::req(filtered_data_rv$samples)
      cols <- colnames(filtered_data_rv$samples)
      if (!is.null(cmp) && length(cmp$included_samples()) > 0L) {
        cols <- unique(c("cmp_group", cols))
      }
      shiny::updateSelectInput(
        session, "metadata_column",
        choices  = cols,
        selected = if ("cmp_group" %in% cols) "cmp_group" else cols[1]
      )
    })
    
    # ---- dynamic top/extra-annotation picker (like ssGSEA) -------------------
    output$ann_cols_ui <- shiny::renderUI({
      shiny::req(filtered_data_rv$samples)
      cols <- colnames(filtered_data_rv$samples)
      if (!is.null(cmp) && length(cmp$included_samples()) > 0L) {
        cols <- unique(c("cmp_group", cols))
      }
      default_primary <- if (!is.null(cmp) && length(cmp$included_samples()) > 0) {
        "cmp_group"
      } else if (!is.null(input$metadata_column) && input$metadata_column %in% cols) {
        input$metadata_column
      } else {
        cols[1]
      }
      shiny::tagList(
        shiny::selectInput(ns("ann_primary"), "Top annotation (primary):",
                           choices = cols, selected = default_primary),
        shiny::selectInput(ns("ann_cols"), "Additional column annotations:",
                           choices = setdiff(cols, default_primary), multiple = TRUE)
      )
    })
    shiny::observeEvent(input$ann_primary, {
      shiny::req(filtered_data_rv$samples)
      cols <- colnames(filtered_data_rv$samples)
      if (!is.null(cmp) && length(cmp$included_samples()) > 0L) {
        cols <- unique(c("cmp_group", cols))
      }
      shiny::updateSelectInput(session, "ann_cols",
                               choices = setdiff(cols, input$ann_primary))
    })
    
    # ---- refresh level choices ----------------------------------------------
    shiny::observeEvent(input$metadata_column, {
      shiny::req(filtered_data_rv$samples, input$metadata_column)
      if (identical(input$metadata_column, "cmp_group") && !is.null(cmp)) {
        lv <- c("Reference", "Test")
      } else {
        lv <- unique(as.character(filtered_data_rv$samples[[input$metadata_column]]))
        lv <- lv[!is.na(lv)]
      }
      shiny::updateSelectInput(session, "reference_condition",
                               choices = lv, selected = if (length(lv)) lv[1] else character(0)
      )
      shiny::updateSelectInput(session, "test_condition",
                               choices = lv, selected = if (length(lv) >= 2) lv[2] else if (length(lv)) lv[1] else character(0)
      )
    }, ignoreInit = FALSE)
    
    # ---- run DE --------------------------------------------------------------
    shiny::observeEvent(input$run_de, {
      shiny::req(filtered_data_rv$counts, filtered_data_rv$samples, filtered_data_rv$species)
      
      use_cmp <- !is.null(cmp) &&
        length(cmp$included_samples()) > 0L &&
        identical(input$metadata_column, "cmp_group")
      
      if (!use_cmp) {
        if (identical(input$reference_condition, input$test_condition)) {
          shiny::showNotification("Reference and Test conditions must be different.", type = "error"); return()
        }
        dds <- DESeq2::DESeqDataSetFromMatrix(
          countData = filtered_data_rv$counts,
          colData   = filtered_data_rv$samples,
          design    = stats::as.formula(paste("~", input$metadata_column))
        )
        dds[[input$metadata_column]] <- stats::relevel(
          factor(dds[[input$metadata_column]]), ref = input$reference_condition
        )
        dds <- suppressMessages(DESeq2::DESeq(dds))
        res <- suppressMessages(DESeq2::results(
          dds, contrast = c(input$metadata_column, input$test_condition, input$reference_condition)
        ))
      } else {
        ids <- cmp$included_samples()
        cf  <- cmp$cmp_factor()
        if (length(ids) != length(cf)) {
          shiny::showNotification("Comparison invalid: sample IDs and roles mismatch.", type = "error"); return()
        }
        counts_sub  <- filtered_data_rv$counts[, ids, drop = FALSE]
        coldata_sub <- filtered_data_rv$samples[ids, , drop = FALSE]
        coldata_sub$cmp_group <- factor(as.character(cf), levels = c("Reference", "Test"))
      
        
        dds <- DESeq2::DESeqDataSetFromMatrix(
          countData = counts_sub,
          colData   = coldata_sub,
          design    = ~ cmp_group
        )
        dds$cmp_group <- stats::relevel(dds$cmp_group, ref = "Reference")
        dds <- suppressMessages(DESeq2::DESeq(dds))
        res <- suppressMessages(DESeq2::results(dds, contrast = c("cmp_group", "Test", "Reference")))
      }
      
      # tidy results + SYMBOL column
      rn <- rownames(res)
      res$symbol <- rn
      if (is_ensembl_id(rn)) {
        conv <- convert_ensembl_to_symbol(rn, filtered_data_rv$species)
        res$symbol <- conv[rn]
      }
      res_df <- as.data.frame(res)
      res_df$symbol[is.na(res_df$symbol)] <- rownames(res_df)[is.na(res_df$symbol)]
      res_df <- res_df[!is.na(res_df$padj) & !is.na(res_df$log2FoldChange), , drop = FALSE]
      res_df <- res_df[order(res_df$padj), , drop = FALSE]
      
      res_reactive(res_df)
      filtered_dds_rv(dds)
    })
    
    # ---- results table -------------------------------------------------------
    output$deTable <- DT::renderDT({
      shiny::req(res_reactive())
      DT::datatable(res_reactive(), options = list(scrollX = TRUE, pageLength = 25))
    })
    
    # ---- heatmap (screen) ----------------------------------------------------
    output$heatmapPlot <- shiny::renderPlot({
      shiny::req(res_reactive(), filtered_data_rv$norm_counts, filtered_data_rv$samples)
      
      res_df <- res_reactive()
      keep <- !is.na(res_df$padj) & !is.na(res_df$log2FoldChange) &
        abs(res_df$log2FoldChange) >= input$lfc_threshold &
        res_df$padj <= input$padj_threshold
      res_df <- res_df[keep, , drop = FALSE]
      shiny::validate(shiny::need(nrow(res_df) >= 2, "Not enough DE genes after filtering for heatmap."))
      
      top_n <- input$num_genes
      top   <- head(res_df[order(res_df$padj), ], top_n)
      
      sel_ids <- rownames(top)
      expr <- filtered_data_rv$norm_counts[sel_ids, , drop = FALSE]
      
      use_cmp_now <- !is.null(cmp) && identical(input$metadata_column, "cmp_group") && length(cmp$included_samples()) > 0
      if (use_cmp_now) {
        expr <- expr[, cmp$included_samples(), drop = FALSE]
      }
      expr <- log2(expr + 1)
      
      # SYMBOL row labels (fallback to IDs)
      lab <- top$symbol
      if (is.null(lab) || !length(lab)) lab <- sel_ids
      lab[is.na(lab) | lab == ""] <- sel_ids[is.na(lab) | lab == ""]
      rownames(expr) <- make.unique(lab)
      
      # samples copy for annotations (inject cmp_group if needed)
      samples_for_ann <- filtered_data_rv$samples
      if (use_cmp_now) {
        role <- setNames(as.character(cmp$cmp_factor()), cmp$included_samples())
        samples_for_ann$cmp_group <- factor(role[rownames(samples_for_ann)],
                                            levels = c("Reference","Test"))
      }
      
      # primary column from UI; fallback to cmp_group/metadata_column
      primary_col <- if (!is.null(input$ann_primary) &&
                         input$ann_primary %in% colnames(samples_for_ann)) {
        input$ann_primary
      } else if (use_cmp_now) {
        "cmp_group"
      } else {
        input$metadata_column
      }
      
      # build top annotation (primary + extra tracks)
      ann_cols <- input$ann_cols %||% character(0)
      ha <- .build_top_annotation(
        samples     = samples_for_ann,
        primary_col = primary_col,
        extra_cols  = ann_cols,
        col_order   = colnames(expr)
      )
      
      hp <- ComplexHeatmap::Heatmap(
        expr,
        name = "log2(norm counts)",
        top_annotation   = ha,
        cluster_rows     = TRUE,
        cluster_columns  = isTRUE(input$cluster_columns),
        show_column_names = FALSE,
        show_row_names   = TRUE,
        row_names_gp     = grid::gpar(fontsize = 6, fontface = "bold"),
        column_title     = paste("Top", min(top_n, nrow(expr)), "Diff Genes"),
        column_title_gp  = grid::gpar(fontface = "bold")
      )
      ComplexHeatmap::draw(hp)
    })
    
    # ---- heatmap (PDF) -------------------------------------------------------
    output$download_heatmap <- shiny::downloadHandler(
      filename = function() paste0("heatmap_", Sys.Date(), ".pdf"),
      content = function(file) {
        grDevices::pdf(file, width = 10, height = 8)
        
        res_df <- res_reactive()
        keep <- !is.na(res_df$padj) & !is.na(res_df$log2FoldChange) &
          abs(res_df$log2FoldChange) >= input$lfc_threshold &
          res_df$padj <= input$padj_threshold
        res_df <- res_df[keep, , drop = FALSE]
        if (nrow(res_df) < 2) {
          grid::grid.newpage(); grid::grid.text("Not enough DE genes for heatmap."); grDevices::dev.off(); return()
        }
        
        top_n <- input$num_genes
        top   <- head(res_df[order(res_df$padj), ], top_n)
        sel_ids <- rownames(top)
        expr <- filtered_data_rv$norm_counts[sel_ids, , drop = FALSE]
        
        use_cmp_now <- !is.null(cmp) && identical(input$metadata_column, "cmp_group") && length(cmp$included_samples()) > 0
        if (use_cmp_now) expr <- expr[, cmp$included_samples(), drop = FALSE]
        expr <- log2(expr + 1)
        
        lab <- top$symbol
        if (is.null(lab) || !length(lab)) lab <- sel_ids
        lab[is.na(lab) | lab == ""] <- sel_ids[is.na(lab) | lab == ""]
        rownames(expr) <- make.unique(lab)
        
        samples_for_ann <- filtered_data_rv$samples
        if (use_cmp_now) {
          role <- setNames(as.character(cmp$cmp_factor()), cmp$included_samples())
          samples_for_ann$cmp_group <- factor(role[rownames(samples_for_ann)],
                                              levels = c("Reference","Test"))
        }
        
        primary_col <- if (!is.null(input$ann_primary) &&
                           input$ann_primary %in% colnames(samples_for_ann)) {
          input$ann_primary
        } else if (use_cmp_now) {
          "cmp_group"
        } else {
          input$metadata_column
        }
        
        ann_cols <- input$ann_cols %||% character(0)
        ha <- .build_top_annotation(
          samples     = samples_for_ann,
          primary_col = primary_col,
          extra_cols  = ann_cols,
          col_order   = colnames(expr)
        )
        
        hp <- ComplexHeatmap::Heatmap(
          expr,
          name = "log2(norm counts)",
          top_annotation   = ha,
          cluster_rows     = TRUE,
          cluster_columns  = isTRUE(input$cluster_columns),
          show_column_names = FALSE,
          show_row_names   = TRUE,
          row_names_gp     = grid::gpar(fontsize = 6, fontface = "bold"),
          column_title     = paste("Top", min(top_n, nrow(expr)), "Diff Genes"),
          column_title_gp  = grid::gpar(fontface = "bold")
        )
        ComplexHeatmap::draw(hp)
        grDevices::dev.off()
      }
    )
    
    # ---- results csv ---------------------------------------------------------
    output$download_de_table <- shiny::downloadHandler(
      filename = function() "differential_expression_results.csv",
      content  = function(file) utils::write.csv(res_reactive(), file, row.names = FALSE)
    )
    
    # ---- public API ----------------------------------------------------------
    list(
      group_var        = shiny::reactive(if (identical(input$metadata_column, "cmp_group")) "cmp_group" else input$metadata_column),
      ref_level        = shiny::reactive(if (identical(input$metadata_column, "cmp_group")) "Reference" else input$reference_condition),
      test_level       = shiny::reactive(if (identical(input$metadata_column, "cmp_group")) "Test" else input$test_condition),
      cmp_factor       = shiny::reactive(if (!is.null(cmp) && identical(input$metadata_column, "cmp_group")) cmp$cmp_factor() else NULL),
      included_samples = shiny::reactive(if (!is.null(cmp) && identical(input$metadata_column, "cmp_group")) cmp$included_samples() else NULL)
    )
  })
}
