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
    .current_ref_label <- function() {
      if (!is.null(cmp) &&
          identical(input$metadata_column, "cmp_group") &&
          length(cmp$included_samples()) > 0L) {
        return(tryCatch(cmp$ref_level(), error = function(e) NULL))
      }
      input$reference_condition
    }
    
    .current_test_label <- function() {
      if (!is.null(cmp) &&
          identical(input$metadata_column, "cmp_group") &&
          length(cmp$included_samples()) > 0L) {
        return(tryCatch(cmp$test_level(), error = function(e) NULL))
      }
      input$test_condition
    }
    
    .current_cmp_tag <- function() {
      test <- .current_test_label()
      ref  <- .current_ref_label()
      if (is.null(test) || !nzchar(as.character(test)) ||
          is.null(ref)  || !nzchar(as.character(ref))) {
        return("contrast")
      }
      paste0(.safe_tag(test), "_vs_", .safe_tag(ref))
    }
    .use_cmp_now <- function() {
      isTRUE(input$use_cmp) &&
        !is.null(cmp) &&
        length(cmp$included_samples()) > 0L
    }
    .heatmap_payload <- function() {
      shiny::req(res_reactive(), filtered_data_rv$norm_counts, filtered_data_rv$samples)
      
      res_df <- res_reactive()
      keep <- !is.na(res_df$padj) & !is.na(res_df$log2FoldChange) &
        abs(res_df$log2FoldChange) >= input$lfc_threshold &
        res_df$padj <= input$padj_threshold
      res_df <- res_df[keep, , drop = FALSE]
      
      shiny::validate(shiny::need(nrow(res_df) >= 2, "Not enough DE genes after filtering for heatmap."))
      
      top_n <- input$num_genes
      top   <- head(res_df[order(res_df$padj), , drop = FALSE], top_n)
      
      sel_ids <- rownames(top)
      expr <- filtered_data_rv$norm_counts[sel_ids, , drop = FALSE]
      
      use_cmp_now <- .use_cmp_now() &&
        identical(input$metadata_column, "cmp_group")
      
      if (use_cmp_now) {
        expr <- expr[, cmp$included_samples(), drop = FALSE]
      }
      
      expr <- log2(expr + 1)
      
      lab <- top$symbol
      if (is.null(lab) || !length(lab)) lab <- sel_ids
      lab[is.na(lab) | lab == ""] <- sel_ids[is.na(lab) | lab == ""]
      rownames(expr) <- make.unique(lab)
      
      samples_for_ann <- filtered_data_rv$samples
      if (use_cmp_now) {
        role <- setNames(as.character(cmp$cmp_factor()), cmp$included_samples())
        samples_for_ann$cmp_group <- factor(
          role[rownames(samples_for_ann)],
          levels = c("Reference", "Test")
        )
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
      
      transpose <- isTRUE(input$transpose_heatmap) && nrow(expr) <= 100
      
      if (isTRUE(input$transpose_heatmap) && nrow(expr) > 100) {
        shiny::showNotification(
          "Heatmap reversal is only available when 100 genes or fewer are displayed.",
          type = "warning"
        )
      }
      
      list(
        expr = expr,
        samples_for_ann = samples_for_ann,
        primary_col = primary_col,
        ann_cols = ann_cols,
        top_n = top_n,
        transpose = transpose
      )
    }
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
    .build_row_annotation <- function(samples, primary_col, extra_cols, row_order) {
      if (!primary_col %in% colnames(samples)) {
        primary_col <- colnames(samples)[1]
      }
      
      ann_df <- data.frame(
        Group = factor(samples[row_order, primary_col, drop = TRUE]),
        check.names = FALSE
      )
      
      extra_cols <- intersect(extra_cols %||% character(0), setdiff(colnames(samples), primary_col))
      for (cn in extra_cols) {
        x <- samples[row_order, cn, drop = TRUE]
        if (is.numeric(x) && length(unique(x)) > 12) {
          q <- cut(
            x,
            breaks = unique(stats::quantile(x, probs = seq(0, 1, 0.25), na.rm = TRUE)),
            include.lowest = TRUE, dig.lab = 6
          )
          ann_df[[cn]] <- q
        } else {
          ann_df[[cn]] <- factor(as.character(x))
        }
      }
      
      col_list <- list()
      glv <- levels(ann_df$Group)
      if (length(glv) <= 8) {
        col_list$Group <- setNames(
          RColorBrewer::brewer.pal(max(3, length(glv)), "Set2")[seq_along(glv)],
          glv
        )
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
          } else {
            grDevices::hcl.colors(length(lv), "Dark 3")
          }
          col_list[[cn]] <- setNames(pal[seq_along(lv)], lv)
        }
      }
      
      ComplexHeatmap::rowAnnotation(df = ann_df, col = col_list)
    }
    
    `%||%` <- function(a, b) if (!is.null(a)) a else b
    
    # ---- populate variable to test ------------------------------------------
    shiny::observe({
      shiny::req(filtered_data_rv$samples)
      cols <- colnames(filtered_data_rv$samples)
      
      if (.use_cmp_now()) {
        cols <- unique(c("cmp_group", cols))
        shiny::updateSelectInput(
          session, "metadata_column",
          choices  = cols,
          selected = "cmp_group"
        )
      } else {
        current <- isolate(input$metadata_column)
        selected_col <- if (!is.null(current) && current %in% cols) current else cols[1]
        
        shiny::updateSelectInput(
          session, "metadata_column",
          choices  = cols,
          selected = selected_col
        )
      }
    })
    
    
    # ---- dynamic top/extra-annotation picker (like ssGSEA) -------------------
    output$ann_cols_ui <- shiny::renderUI({
      shiny::req(filtered_data_rv$samples)
      cols <- colnames(filtered_data_rv$samples)
      if (.use_cmp_now()) {
        cols <- unique(c("cmp_group", cols))
      }
     default_primary <- if (.use_cmp_now()) {
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
      if (.use_cmp_now()) {
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
      
      use_cmp <- .use_cmp_now()
      
      if (!use_cmp) {
        shiny::req(input$metadata_column, input$reference_condition, input$test_condition)
        
        if (!input$metadata_column %in% colnames(filtered_data_rv$samples)) {
          shiny::showNotification("Selected metadata column is not present in sample metadata.", type = "error")
          return()
        }
        
        grp <- filtered_data_rv$samples[[input$metadata_column]]
        if (all(is.na(grp))) {
          shiny::showNotification("Selected metadata column contains only missing values.", type = "error")
          return()
        }
        
        grp_chr <- as.character(grp)
        if (!input$reference_condition %in% grp_chr) {
          shiny::showNotification("Reference condition is not present in the selected metadata column.", type = "error")
          return()
        }
        if (!input$test_condition %in% grp_chr) {
          shiny::showNotification("Test condition is not present in the selected metadata column.", type = "error")
          return()
        }
        if (identical(input$reference_condition, input$test_condition)) {
          shiny::showNotification("Reference and Test conditions must be different.", type = "error"); return()
        }
        keep <- !is.na(filtered_data_rv$samples[[input$metadata_column]])
        keep <- keep & as.character(filtered_data_rv$samples[[input$metadata_column]]) %in%
          c(input$reference_condition, input$test_condition)
        
        if (sum(keep) < 2) {
          shiny::showNotification("Not enough samples remain after filtering to the selected reference/test conditions.", type = "error")
          return()
        }
        
        counts_sub  <- filtered_data_rv$counts[, keep, drop = FALSE]
        coldata_sub <- filtered_data_rv$samples[keep, , drop = FALSE]
        
        dds <- DESeq2::DESeqDataSetFromMatrix(
          countData = counts_sub,
          colData   = coldata_sub,
          design    = stats::as.formula(paste("~", input$metadata_column))
        )
        dds[[input$metadata_column]] <- stats::relevel(
          droplevels(factor(dds[[input$metadata_column]])),
          ref = input$reference_condition
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
        role_vec <- as.character(cf)
        names(role_vec) <- ids  # just in case not already named
        coldata_sub$cmp_group <- droplevels(
          factor(role_vec[rownames(coldata_sub)], levels = c("Reference", "Test"))
        ) 
        message("=== DEBUG: cmp_group column ===")
        print(table(coldata_sub$cmp_group))
        print(head(coldata_sub[, c("cmp_group")], 10))
        message("=== DEBUG: cmp_factor() ===")
        print(head(cmp$cmp_factor()))
        print(table(cmp$cmp_factor()))
        ids <- cmp$included_samples()
        cf <- cmp$cmp_factor()
        stopifnot(all(ids %in% colnames(filtered_data_rv$counts)))  # critical
        stopifnot(length(ids) == length(cf))  # required
        
        # Check matching
        true_group <- NULL
        if (!is.null(input$metadata_column) &&
            input$metadata_column %in% colnames(filtered_data_rv$samples)) {
          true_group <- as.character(filtered_data_rv$samples[ids, input$metadata_column, drop = TRUE])
        }
        
        role_df <- data.frame(
          sample = ids,
          role = as.character(cf),
          true_group = true_group,
          stringsAsFactors = FALSE
        )
        print(head(role_df))
        
        
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
      hpdat <- .heatmap_payload()
      
      plot_mat <- if (isTRUE(hpdat$transpose)) t(hpdat$expr) else hpdat$expr
      
      if (isTRUE(hpdat$transpose)) {
        row_ha <- .build_row_annotation(
          samples     = hpdat$samples_for_ann,
          primary_col = hpdat$primary_col,
          extra_cols  = hpdat$ann_cols,
          row_order   = rownames(plot_mat)
        )
        
        hp <- ComplexHeatmap::Heatmap(
          plot_mat,
          name = "log2(norm counts)",
          left_annotation   = row_ha,
          cluster_rows      = isTRUE(input$cluster_columns),
          cluster_columns   = TRUE,
          show_row_names    = FALSE,
          show_column_names = TRUE,
          column_names_gp   = grid::gpar(fontsize = 6, fontface = "bold"),
          row_title         = "Samples",
          column_title      = paste("Top", min(hpdat$top_n, ncol(plot_mat)), "Diff Genes"),
          column_title_gp   = grid::gpar(fontface = "bold")
        )
      } else {
        ha <- .build_top_annotation(
          samples     = hpdat$samples_for_ann,
          primary_col = hpdat$primary_col,
          extra_cols  = hpdat$ann_cols,
          col_order   = colnames(plot_mat)
        )
        
        hp <- ComplexHeatmap::Heatmap(
          plot_mat,
          name = "log2(norm counts)",
          top_annotation    = ha,
          cluster_rows      = TRUE,
          cluster_columns   = isTRUE(input$cluster_columns),
          show_column_names = FALSE,
          show_row_names    = TRUE,
          row_names_gp      = grid::gpar(fontsize = 6, fontface = "bold"),
          column_title      = paste("Top", min(hpdat$top_n, nrow(plot_mat)), "Diff Genes"),
          column_title_gp   = grid::gpar(fontface = "bold")
        )
      }
      
      ComplexHeatmap::draw(hp)
    })
    # ---- heatmap (PDF) -------------------------------------------------------
    output$download_heatmap <- shiny::downloadHandler(
      filename = function() paste0("heatmap_", .current_cmp_tag(), ".pdf"),
      content = function(file) {
        grDevices::pdf(file, width = 10, height = 8)
        
        hpdat <- tryCatch(.heatmap_payload(), error = function(e) NULL)
        if (is.null(hpdat)) {
          grid::grid.newpage()
          grid::grid.text("Not enough DE genes for heatmap.")
          grDevices::dev.off()
          return()
        }
        
        plot_mat <- if (isTRUE(hpdat$transpose)) t(hpdat$expr) else hpdat$expr
        
        if (isTRUE(hpdat$transpose)) {
          row_ha <- .build_row_annotation(
            samples     = hpdat$samples_for_ann,
            primary_col = hpdat$primary_col,
            extra_cols  = hpdat$ann_cols,
            row_order   = rownames(plot_mat)
          )
          
          hp <- ComplexHeatmap::Heatmap(
            plot_mat,
            name = "log2(norm counts)",
            left_annotation   = row_ha,
            cluster_rows      = isTRUE(input$cluster_columns),
            cluster_columns   = TRUE,
            show_row_names    = FALSE,
            show_column_names = TRUE,
            column_names_gp   = grid::gpar(fontsize = 6, fontface = "bold"),
            row_title         = "Samples",
            column_title      = paste("Top", min(hpdat$top_n, ncol(plot_mat)), "Diff Genes"),
            column_title_gp   = grid::gpar(fontface = "bold")
          )
        } else {
          ha <- .build_top_annotation(
            samples     = hpdat$samples_for_ann,
            primary_col = hpdat$primary_col,
            extra_cols  = hpdat$ann_cols,
            col_order   = colnames(plot_mat)
          )
          
          hp <- ComplexHeatmap::Heatmap(
            plot_mat,
            name = "log2(norm counts)",
            top_annotation    = ha,
            cluster_rows      = TRUE,
            cluster_columns   = isTRUE(input$cluster_columns),
            show_column_names = FALSE,
            show_row_names    = TRUE,
            row_names_gp      = grid::gpar(fontsize = 6, fontface = "bold"),
            column_title      = paste("Top", min(hpdat$top_n, nrow(plot_mat)), "Diff Genes"),
            column_title_gp   = grid::gpar(fontface = "bold")
          )
        }
        
        ComplexHeatmap::draw(hp)
        grDevices::dev.off()
      }
    )
    output$download_heatmap_matrix <- shiny::downloadHandler(
      filename = function() paste0("heatmap_matrix_", .current_cmp_tag(), ".csv"),
      content = function(file) {
        hpdat <- .heatmap_payload()
        out_mat <- if (isTRUE(hpdat$transpose)) t(hpdat$expr) else hpdat$expr
        utils::write.csv(
          as.data.frame(out_mat),
          file = file,
          row.names = TRUE
        )
      }
    )
    # ---- results csv ---------------------------------------------------------
    output$download_de_table <- shiny::downloadHandler(
      filename = function() paste0("differential_expression_", .current_cmp_tag(), "_results.csv"),
      content  = function(file) utils::write.csv(res_reactive(), file, row.names = FALSE)
    )
    
    # ---- public API ----------------------------------------------------------
    list(
      group_var        = shiny::reactive(if (identical(input$metadata_column, "cmp_group")) "cmp_group" else input$metadata_column),
      ref_level        = shiny::reactive(.current_ref_label()),
      test_level       = shiny::reactive(.current_test_label()),
      cmp_factor       = shiny::reactive(if (.use_cmp_now() && identical(input$metadata_column, "cmp_group")) cmp$cmp_factor() else NULL),
      included_samples = shiny::reactive(if (.use_cmp_now() && identical(input$metadata_column, "cmp_group")) cmp$included_samples() else NULL)
    )
  })
}
