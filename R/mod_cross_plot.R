# === Module: mod_cross_server ===
#' Cross Plot and Venn Diagram Server
#'
#' Tabs: Cross Plot, Venn, Heatmap (with Ensembl->SYMBOL labeling), Pathway Dotplot.
#'
#' @param id module id (must match mod_cross_tab_ui)
#' @param filtered_data_rv reactiveValues with $counts, $samples, $norm_counts, $species
#' @param filtered_dds_rv reactiveVal with a DESeqDataSet (subset already applied)
#' @importFrom shiny moduleServer req observe observeEvent updateSelectInput
#' @importFrom shiny showNotification reactive reactiveVal renderPlot downloadHandler
#' @importFrom ggplot2 ggplot aes geom_point geom_vline geom_hline theme_minimal labs annotate scale_color_manual ggsave
#' @importFrom ggrepel geom_text_repel
#' @importFrom gridExtra grid.arrange
#' @importFrom grid grid.draw grid.grabExpr gpar
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation draw
#' @importFrom DESeq2 DESeq DESeqDataSetFromMatrix results resultsNames
#' @importFrom stats as.formula cor complete.cases
#' @importFrom utils write.csv
#' @export
mod_cross_server <- function(id, filtered_data_rv, filtered_dds_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    crossplot_data <- shiny::reactiveVal(NULL)
    
    # 1) Populate metadata selects ----
    shiny::observe({
      shiny::req(filtered_data_rv$samples)
      cols <- colnames(filtered_data_rv$samples)
      shiny::updateSelectInput(session, "metadata_column_x",
                               choices = cols, selected = if (length(cols)) cols[1] else character(0))
      shiny::updateSelectInput(session, "metadata_column_y",
                               choices = cols, selected = if (length(cols) >= 2) cols[2] else if (length(cols)) cols[1] else character(0))
    })
    
    # X levels
    shiny::observeEvent(input$metadata_column_x, {
      shiny::req(filtered_data_rv$samples, input$metadata_column_x)
      lv <- unique(as.character(filtered_data_rv$samples[[input$metadata_column_x]]))
      lv <- lv[!is.na(lv)]
      shiny::updateSelectInput(session, "reference_condition_x",
                               choices = lv, selected = if (length(lv)) lv[1] else character(0))
      shiny::updateSelectInput(session, "test_condition_x",
                               choices = lv, selected = if (length(lv) >= 2) lv[2] else if (length(lv)) lv[1] else character(0))
    }, ignoreInit = FALSE)
    
    # Y levels
    shiny::observeEvent(input$metadata_column_y, {
      shiny::req(filtered_data_rv$samples, input$metadata_column_y)
      lv <- unique(as.character(filtered_data_rv$samples[[input$metadata_column_y]]))
      lv <- lv[!is.na(lv)]
      shiny::updateSelectInput(session, "reference_condition_y",
                               choices = lv, selected = if (length(lv)) lv[1] else character(0))
      shiny::updateSelectInput(session, "test_condition_y",
                               choices = lv, selected = if (length(lv) >= 2) lv[2] else if (length(lv)) lv[1] else character(0))
    }, ignoreInit = FALSE)
    
    # One-time hint
    shiny::observeEvent(input$run_crossplot, {
      shiny::showNotification(
        "Tip: add highlight genes (space-separated) before running if you want labels.",
        type = "message", duration = 5
      )
    }, once = TRUE)
    
    # ---- 2) Run contrasts and merge ----
    shiny::observeEvent(input$run_crossplot, {
      shiny::req(filtered_dds_rv(), filtered_data_rv,
                 input$metadata_column_x, input$reference_condition_x, input$test_condition_x,
                 input$metadata_column_y, input$reference_condition_y, input$test_condition_y)
      
      if (identical(input$reference_condition_x, input$test_condition_x)) {
        shiny::showNotification("X-axis: Reference and Test must be different.", type = "error"); return()
      }
      if (identical(input$reference_condition_y, input$test_condition_y)) {
        shiny::showNotification("Y-axis: Reference and Test must be different.", type = "error"); return()
      }
      
      shiny::showNotification("Running Cross Plot analysis...", type = "message")
      
      samples <- filtered_data_rv$samples
      counts  <- filtered_data_rv$counts
      
      # Compose contrasts
      x_comp <- c(input$metadata_column_x, input$test_condition_x, input$reference_condition_x)
      y_comp <- c(input$metadata_column_y, input$test_condition_y, input$reference_condition_y)
      
      # Build DE objects as needed
      if (!paste(x_comp, collapse = "_") %in% DESeq2::resultsNames(filtered_dds_rv())) {
        dds_x <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = samples,
                                                design = stats::as.formula(paste("~", input$metadata_column_x)))
        dds_x <- DESeq2::DESeq(dds_x)
      } else dds_x <- filtered_dds_rv()
      
      if (!paste(y_comp, collapse = "_") %in% DESeq2::resultsNames(filtered_dds_rv())) {
        dds_y <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = samples,
                                                design = stats::as.formula(paste("~", input$metadata_column_y)))
        dds_y <- DESeq2::DESeq(dds_y)
      } else dds_y <- filtered_dds_rv()
      
      # Results
      res_x <- DESeq2::results(dds_x, contrast = x_comp)
      res_y <- DESeq2::results(dds_y, contrast = y_comp)
      res_x <- res_x[!is.na(res_x$padj) & !is.na(res_x$log2FoldChange), ]
      res_y <- res_y[!is.na(res_y$padj) & !is.na(res_y$log2FoldChange), ]
      
      df_x <- data.frame(gene = rownames(res_x),
                         log2FoldChange_x = res_x$log2FoldChange,
                         padj_x           = res_x$padj,
                         stringsAsFactors = FALSE)
      df_y <- data.frame(gene = rownames(res_y),
                         log2FoldChange_y = res_y$log2FoldChange,
                         padj_y           = res_y$padj,
                         stringsAsFactors = FALSE)
      
      merged_df <- merge(df_x, df_y, by = "gene", all = TRUE)
      
      # Labels
      nlab <- tryCatch({ n <- as.numeric(input$crossplot_topgenes); if (is.na(n) || n < 0) 0 else n }, error = function(e) 0)
      merged_df$combined_padj <- pmax(merged_df$padj_x, merged_df$padj_y, na.rm = TRUE)
      
      top_genes <- character(0)
      if (nlab > 0 && nrow(merged_df)) {
        ord <- merged_df[order(merged_df$combined_padj), "gene"]
        top_genes <- head(ord, nlab)
        if (is_ensembl_id(top_genes)) {
          conv <- convert_ensembl_to_symbol(top_genes, filtered_data_rv$species)
          conv[is.na(conv)] <- top_genes[is.na(conv)]
          top_genes <- conv
        }
      }
      
      highlight <- input$crossplot_gene_label
      highlight <- if (nchar(trimws(highlight)) > 0) unlist(stringr::str_split(highlight, "[\\s,]+")) else character(0)
      
      merged_df$label <- ifelse(merged_df$gene %in% union(top_genes, highlight), merged_df$gene, NA)
      merged_df <- merged_df[order(merged_df$padj_x), ]
      crossplot_data(merged_df)
    })
    
    # ---- 3) Cross plot tab ----
    generate_cross_plot <- function(df, n_keep, tx, rx, ty, ry) {
      df$combined_padj <- pmax(df$padj_x, df$padj_y, na.rm = TRUE)
      df$category <- "Other"
      df$category[df$log2FoldChange_x >=  1 & df$log2FoldChange_y >=  1] <- "Up-Up"
      df$category[df$log2FoldChange_x <= -1 & df$log2FoldChange_y <= -1] <- "Down-Down"
      df$category[df$log2FoldChange_x >=  1 & df$log2FoldChange_y <= -1] <- "Up-Down"
      df$category[df$log2FoldChange_x <= -1 & df$log2FoldChange_y >=  1] <- "Down-Up"
      df$category[abs(df$log2FoldChange_x) >  1 & abs(df$log2FoldChange_y) <  1] <- "Comp1-only"
      df$category[abs(df$log2FoldChange_x) <  1 & abs(df$log2FoldChange_y) >  1] <- "Comp2-only"
      
      df <- head(df, n_keep)
      
      pearson  <- round(stats::cor(df$log2FoldChange_x, df$log2FoldChange_y, method = "pearson",  use = "complete.obs"), 3)
      spearman <- round(stats::cor(df$log2FoldChange_x, df$log2FoldChange_y, method = "spearman", use = "complete.obs"), 3)
      
      p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$log2FoldChange_x, y = .data$log2FoldChange_y, color = .data$category)) +
        ggplot2::annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -1.1,
                          label = paste0("Pearson r = ", pearson, "   Spearman r = ", spearman),
                          size = 4) +
        ggplot2::geom_point(alpha = 0.7) +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "black") +
        ggplot2::geom_vline(xintercept = c(-1, 1), linetype = "solid", color = "blue") +
        ggplot2::geom_hline(yintercept = c(-1, 1), linetype = "solid", color = "blue") +
        ggplot2::scale_color_manual(values = c(
          "Up-Up" = "firebrick",
          "Down-Down" = "royalblue",
          "Up-Down" = "goldenrod",
          "Down-Up" = "purple",
          "Comp1-only" = "darkorange",
          "Comp2-only" = "darkgreen",
          "Other" = "gray"
        )) +
        ggplot2::theme_minimal() +
        ggplot2::labs(
          title = "Cross Plot: log2FC vs log2FC",
          x = paste("log2FC:", tx, "vs", rx),
          y = paste("log2FC:", ty, "vs", ry),
          color = "Regulation"
        )
      
      if ("label" %in% names(df) && any(!is.na(df$label) & df$label != "")) {
        lab_df <- df[!is.na(df$label) & df$label != "", ]
        if (nrow(lab_df)) {
          p <- p + ggrepel::geom_text_repel(
            data = lab_df, ggplot2::aes(label = .data$label),
            max.overlaps = Inf, size = 5, box.padding = 0.1
          )
        }
      }
      p
    }
    
    output$crossPlot <- shiny::renderPlot({
      shiny::req(crossplot_data())
      p <- generate_cross_plot(
        crossplot_data(),
        input$crossplot_gene_count,
        input$test_condition_x,  input$reference_condition_x,
        input$test_condition_y,  input$reference_condition_y
      )
      print(p)
    })
    
    output$download_cross_plot <- shiny::downloadHandler(
      filename = function() input$crossplot_filename,
      content  = function(file) {
        shiny::req(crossplot_data())
        p <- generate_cross_plot(
          crossplot_data(),
          input$crossplot_gene_count,
          input$test_condition_x,  input$reference_condition_x,
          input$test_condition_y,  input$reference_condition_y
        )
        ggplot2::ggsave(file, p, width = 10, height = 8)
      }
    )
    
    # ---- 4) Pathway Dotplot tab ----
    generate_cross_pathway_plot <- function(df, species) {
      if (is.null(species) || is.na(species) || identical(species, "")) {
        species <- "Homo sapiens"
        shiny::showNotification("Species not specified. Defaulting to Homo sapiens.", type = "warning")
      }
      orgdb <- get_orgdb(species)
      
      df$combined_padj <- pmax(df$padj_x, df$padj_y, na.rm = TRUE)
      df$category <- "Other"
      df$category[df$log2FoldChange_x >=  1 & df$log2FoldChange_y >=  1] <- "Up-Up"
      df$category[df$log2FoldChange_x <= -1 & df$log2FoldChange_y <= -1] <- "Down-Down"
      df$category[df$log2FoldChange_x >=  1 & df$log2FoldChange_y <= -1] <- "Up-Down"
      df$category[df$log2FoldChange_x <= -1 & df$log2FoldChange_y >=  1] <- "Down-Up"
      df <- df[df$category != "Other", , drop = FALSE]
      if (!nrow(df)) return(NULL)
      
      if (is_symbol(df$gene)) {
        df_ids <- clusterProfiler::bitr(df$gene, fromType = "SYMBOL",  toType = "ENTREZID", OrgDb = orgdb)
      } else {
        df_ids <- clusterProfiler::bitr(df$gene, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = orgdb)
      }
      df_merged <- merge(df, df_ids, by.x = "gene", by.y = 1)
      df_merged <- dplyr::distinct(df_merged, ENTREZID, .keep_all = TRUE)
      df_merged2 <- df_merged[, c("ENTREZID", "log2FoldChange_x", "category")]
      df_merged2 <- df_merged2[stats::complete.cases(df_merged2), , drop = FALSE]
      if (!nrow(df_merged2)) return(NULL)
      
      m <- input$cross_enrich_method
      if (m == "enrichGO") {
        formula_res <- clusterProfiler::compareCluster(
          ENTREZID ~ category, data = df_merged2, fun = "enrichGO",
          OrgDb = orgdb, keyType = "ENTREZID", ont = "BP",
          pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE
        )
      } else if (m == "groupGO") {
        formula_res <- clusterProfiler::compareCluster(
          ENTREZID ~ category, data = df_merged2, fun = "groupGO",
          OrgDb = orgdb, keyType = "ENTREZID", ont = "BP", readable = TRUE
        )
      } else if (m == "enrichKEGG") {
        kegg_sp <- get_kegg_code(species); if (is.null(kegg_sp)) return(NULL)
        formula_res <- clusterProfiler::compareCluster(
          ENTREZID ~ category, data = df_merged2, fun = "enrichKEGG",
          organism = kegg_sp, pvalueCutoff = 0.05, qvalueCutoff = 0.2
        )
      } else if (m == "enrichPathway") {
        rc <- get_reactome_code(species); if (is.null(rc)) return(NULL)
        formula_res <- clusterProfiler::compareCluster(
          ENTREZID ~ category, data = df_merged2, fun = ReactomePA::enrichPathway,
          organism = rc, pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE
        )
      } else {
        shiny::showNotification("Unsupported enrichment method selected.", type = "error")
        return(NULL)
      }
      
      if (is.null(formula_res) || !nrow(as.data.frame(formula_res))) return(NULL)
      enrichplot::dotplot(formula_res, x = "category")
    }
    
    output$crosspathplot <- shiny::renderPlot({
      shiny::req(crossplot_data(), filtered_data_rv$species)
      p <- generate_cross_pathway_plot(crossplot_data(), filtered_data_rv$species)
      if (!is.null(p)) print(p)
    })
    
    output$download_cross_pathway_plot <- shiny::downloadHandler(
      filename = function() {
        ext <- tools::file_ext(input$crosspath_filename)
        if (identical(ext, "")) paste0("cross_pathway_", Sys.Date(), ".pdf") else input$crosspath_filename
      },
      content = function(file) {
        p <- generate_cross_pathway_plot(crossplot_data(), filtered_data_rv$species)
        if (!is.null(p)) ggplot2::ggsave(file, plot = p, width = 10, height = 8, units = "in")
      }
    )
    
    # ---- 5) Venn tab ----
    generate_cross_venn_plot <- function(df) {
      lfc_cutoff <- 1; padj_cutoff <- 0.05
      up_x   <- df$gene[df$log2FoldChange_x >  lfc_cutoff & df$padj_x < padj_cutoff]
      up_y   <- df$gene[df$log2FoldChange_y >  lfc_cutoff & df$padj_y < padj_cutoff]
      down_x <- df$gene[df$log2FoldChange_x < -lfc_cutoff & df$padj_x < padj_cutoff]
      down_y <- df$gene[df$log2FoldChange_y < -lfc_cutoff & df$padj_y < padj_cutoff]
      
      up_venn <- VennDiagram::venn.diagram(
        x = list(X_Up = up_x, Y_Up = up_y),
        filename = NULL, fill = c("darkorange", "darkgreen"), alpha = 0.5,
        main = "Upregulated Genes Venn", disable.logging = TRUE
      )
      down_venn <- VennDiagram::venn.diagram(
        x = list(X_Down = down_x, Y_Down = down_y),
        filename = NULL, fill = c("royalblue", "purple"), alpha = 0.5,
        main = "Downregulated Genes Venn", disable.logging = TRUE
      )
      list(up_venn = up_venn, down_venn = down_venn)
    }
    
    output$crossVennPlot <- shiny::renderPlot({
      shiny::req(crossplot_data())
      vp <- generate_cross_venn_plot(crossplot_data())
      gridExtra::grid.arrange(
        grid::grid.grabExpr(grid::grid.draw(vp$up_venn)),
        grid::grid.grabExpr(grid::grid.draw(vp$down_venn)),
        ncol = 2
      )
    })
    
    output$download_cross_venn_plot <- shiny::downloadHandler(
      filename = function() paste0("cross_venn_plot_", Sys.Date(), ".pdf"),
      content  = function(file) {
        grDevices::pdf(file, width = 10, height = 6)
        vp <- generate_cross_venn_plot(crossplot_data())
        gridExtra::grid.arrange(
          grid::grid.grabExpr(grid::grid.draw(vp$up_venn)),
          grid::grid.grabExpr(grid::grid.draw(vp$down_venn)),
          ncol = 2
        )
        grDevices::dev.off()
      }
    )
    
    output$download_overlap_genes <- shiny::downloadHandler(
      filename = function() {
        paste0("cross_plot_overlap_genes_", Sys.Date(), "_",
               input$test_condition_x, "_", input$test_condition_y, ".csv")
      },
      content = function(file) {
        shiny::req(crossplot_data())
        df <- crossplot_data()
        lfc_cutoff <- 1; padj_cutoff <- 0.05
        up_x   <- df$gene[df$log2FoldChange_x >  lfc_cutoff & df$padj_x < padj_cutoff]
        up_y   <- df$gene[df$log2FoldChange_y >  lfc_cutoff & df$padj_y < padj_cutoff]
        down_x <- df$gene[df$log2FoldChange_x < -lfc_cutoff & df$padj_x < padj_cutoff]
        down_y <- df$gene[df$log2FoldChange_y < -lfc_cutoff & df$padj_y < padj_cutoff]
        
        out <- rbind(
          data.frame(Gene = intersect(up_x, up_y),   Direction = "Upregulated"),
          data.frame(Gene = intersect(down_x, down_y), Direction = "Downregulated")
        )
        utils::write.csv(out, file, row.names = FALSE)
      }
    )
    
    # ---- 6) Heatmap tab (with Ensembl->Symbol row labels) ----
    generate_cross_heat <- function(df, norm_counts, samples, meta_col_x, species) {
      df$category <- "Other"
      df$category[df$log2FoldChange_x >=  1 & df$log2FoldChange_y >=  1] <- "Up-Up"
      df$category[df$log2FoldChange_x <= -1 & df$log2FoldChange_y <= -1] <- "Down-Down"
      df$category[df$log2FoldChange_x >=  1 & df$log2FoldChange_y <= -1] <- "Up-Down"
      df$category[df$log2FoldChange_x <= -1 & df$log2FoldChange_y >=  1] <- "Down-Up"
      
      sel <- unique(df$gene[!(df$category %in% c("Other", "Comp1-only", "Comp2-only"))])
      if (!length(sel)) return(NULL)
      
      # Subset matrix and log transform
      expr_mat <- log2(norm_counts[sel, , drop = FALSE] + 1)
      
      # --- Row label conversion: Ensembl -> Symbol (fallback to Ensembl) ---
      if (is_ensembl_id(sel)) {
        sym <- convert_ensembl_to_symbol(sel, species)
        row_labels <- ifelse(is.na(sym) | sym == "", sel, sym)
        rownames(expr_mat) <- make.unique(row_labels)
      } else {
        # already symbols; just ensure uniqueness
        rownames(expr_mat) <- make.unique(rownames(expr_mat))
      }
      
      # Annotation by selected X metadata
      group_values <- as.character(samples[[meta_col_x]])
      group_levels <- unique(group_values)
      pal <- if (length(group_levels) <= 8) {
        setNames(RColorBrewer::brewer.pal(max(length(group_levels), 3), "Set2")[seq_along(group_levels)], group_levels)
      } else {
        setNames(colorspace::rainbow_hcl(length(group_levels)), group_levels)
      }
      
      ha <- ComplexHeatmap::HeatmapAnnotation(
        df = data.frame(Group = factor(group_values, levels = group_levels)),
        col = list(Group = pal)
      )
      
      ComplexHeatmap::Heatmap(
        expr_mat,
        name = "log2(norm counts)",
        top_annotation = ha,
        cluster_columns = FALSE,
        cluster_rows = TRUE,
        show_column_names = FALSE,
        row_names_gp = grid::gpar(fontsize = 4, fontface = "bold"),
        column_title = "Crossplot Gene Set",
        column_title_gp = grid::gpar(fontsize = 10, fontface = "bold")
      )
    }
    
    output$crossCategoryHeatmap <- shiny::renderPlot({
      shiny::req(crossplot_data(), filtered_data_rv$norm_counts, filtered_data_rv$samples)
      hp <- generate_cross_heat(crossplot_data(), filtered_data_rv$norm_counts,
                                filtered_data_rv$samples, input$metadata_column_x, filtered_data_rv$species)
      if (is.null(hp)) {
        grid::grid.newpage(); grid::grid.text("No genes for heatmap.")
      } else {
        ComplexHeatmap::draw(hp)
      }
    })
    
    output$download_cross_category_heatmap <- shiny::downloadHandler(
      filename = function() paste0("cross_category_heatmap_", Sys.Date(), ".pdf"),
      content  = function(file) {
        shiny::req(crossplot_data(), filtered_data_rv$norm_counts, filtered_data_rv$samples)
        hp <- generate_cross_heat(crossplot_data(), filtered_data_rv$norm_counts,
                                  filtered_data_rv$samples, input$metadata_column_x, filtered_data_rv$species)
        grDevices::pdf(file, width = 8.5, height = 11)
        if (is.null(hp)) {
          grid::grid.newpage(); grid::grid.text("No genes for heatmap.")
        } else {
          ComplexHeatmap::draw(hp)
        }
        grDevices::dev.off()
      }
    )
  })
}
