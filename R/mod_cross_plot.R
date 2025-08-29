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
    
    # ---------- helpers ----------
    `%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a
    sanitize_vec <- function(x) {
      x <- as.character(x)
      x <- x[!is.na(x) & nzchar(x)]
      unique(x)
    }
    parse_gene_ratio <- function(gr) {
      if (!length(gr)) return(NA_real_)
      if (is.numeric(gr)) return(as.numeric(gr))
      num <- suppressWarnings(as.numeric(sub("/.*", "", gr)))
      den <- suppressWarnings(as.numeric(sub(".*/", "", gr)))
      ifelse(!is.na(num) & !is.na(den) & den > 0, num/den, NA_real_)
    }
    # prefer app helpers if present
    is_symbol <- if (exists("is_symbol", mode = "function")) get("is_symbol") else function(ids) !grepl("^ENS[A-Z]*\\d+", ids)
    is_ensembl_id <- if (exists("is_ensembl_id", mode = "function")) get("is_ensembl_id") else function(ids) grepl("^ENS[A-Z]*\\d+", ids)
    get_orgdb <- if (exists("get_orgdb", mode = "function")) get("get_orgdb") else function(species) {
      sp <- as.character(species)
      if (sp %in% c("Homo sapiens","human")) { if (requireNamespace("org.Hs.eg.db", quietly=TRUE)) return(org.Hs.eg.db::org.Hs.eg.db) }
      if (sp %in% c("Mus musculus","mouse")) { if (requireNamespace("org.Mm.eg.db", quietly=TRUE)) return(org.Mm.eg.db::org.Mm.eg.db) }
      if (sp %in% c("Rattus norvegicus","rat")) { if (requireNamespace("org.Rn.eg.db", quietly=TRUE)) return(org.Rn.eg.db::org.Rn.eg.db) }
      NULL
    }
    get_kegg_code <- if (exists("get_kegg_code", mode = "function")) get("get_kegg_code") else function(species) {
      sp <- as.character(species)
      if (sp %in% c("Homo sapiens","human")) return("hsa")
      if (sp %in% c("Mus musculus","mouse")) return("mmu")
      if (sp %in% c("Rattus norvegicus","rat")) return("rno")
      if (sp %in% c("Danio rerio")) return("dre")
      if (sp %in% c("Bos taurus")) return("bta")
      NULL
    }
    get_reactome_code <- if (exists("get_reactome_code", mode = "function")) get("get_reactome_code") else function(species) {
      sp <- as.character(species)
      if (sp %in% c("Homo sapiens","human")) return("human")
      if (sp %in% c("Mus musculus","mouse")) return("mouse")
      if (sp %in% c("Rattus norvegicus","rat")) return("rat")
      NULL
    }
    
    # Safe dotplot for compareCluster (no crash if p.adjust missing)
    safe_comparecluster_dotplot <- function(cc, showCategory = 10, x = "Cluster") {
      if (is.null(cc)) return(NULL)
      df <- as.data.frame(cc)
      if (!nrow(df)) return(NULL)
      
      # Ensure columns needed for plotting exist
      if (!"p.adjust" %in% names(df)) {
        if ("pvalue" %in% names(df)) df$p.adjust <- df$pvalue else df$p.adjust <- NA_real_
      }
      # Count & GeneRatio
      if (!"Count" %in% names(df) && "core_enrichment" %in% names(df)) {
        df$Count <- vapply(strsplit(df$core_enrichment, "/"), length, integer(1))
      }
      if (!"GeneRatio" %in% names(df) || all(is.na(df$GeneRatio))) {
        if ("setSize" %in% names(df) && "Count" %in% names(df)) {
          df$GeneRatio <- paste0(df$Count, "/", df$setSize)
        } else {
          df$GeneRatio <- NA_character_
        }
      }
      df$GeneRatioNum <- parse_gene_ratio(df$GeneRatio)
      
      # Order by significance within each cluster
      if ("Cluster" %in% names(df)) {
        df <- df[order(df$Cluster, df$p.adjust, df$pvalue %||% Inf), ]
        # keep up to showCategory per cluster
        df <- do.call(rbind, lapply(split(df, df$Cluster), function(s) head(s, showCategory)))
      } else {
        df <- head(df[order(df$p.adjust, df$pvalue %||% Inf), ], showCategory)
      }
      
      # Manual ggplot (avoids enrichplot aesthetics issues)
      color_var <- if (all(is.na(df$p.adjust))) NULL else "p.adjust"
      p <- ggplot2::ggplot(
        df,
        ggplot2::aes(
          x = .data[[x]],
          y = reorder(.data$Description, .data$p.adjust, function(v) -rank(v, na.last = "keep")),
          size = .data$Count
        )
      ) +
        ggplot2::geom_point(ggplot2::aes(color = !!(if (!is.null(color_var)) rlang::sym(color_var) else rlang::sym("Count")))) +
        ggplot2::labs(x = x, y = NULL, size = "Gene Count", color = if (is.null(color_var)) "Count" else "Adj. P") +
        ggplot2::theme_minimal(base_size = 11)
      if (!is.null(color_var)) {
        p <- p + ggplot2::scale_color_continuous(low = "#276EF1", high = "#F04438", trans = "reverse")
      }
      p
    }
    
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
      lv <- sanitize_vec(filtered_data_rv$samples[[input$metadata_column_x]])
      shiny::updateSelectInput(session, "reference_condition_x",
                               choices = lv, selected = if (length(lv)) lv[1] else character(0))
      shiny::updateSelectInput(session, "test_condition_x",
                               choices = lv, selected = if (length(lv) >= 2) lv[2] else if (length(lv)) lv[1] else character(0))
    }, ignoreInit = FALSE)
    
    # Y levels
    shiny::observeEvent(input$metadata_column_y, {
      shiny::req(filtered_data_rv$samples, input$metadata_column_y)
      lv <- sanitize_vec(filtered_data_rv$samples[[input$metadata_column_y]])
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
      
      x_comp <- c(input$metadata_column_x, input$test_condition_x, input$reference_condition_x)
      y_comp <- c(input$metadata_column_y, input$test_condition_y, input$reference_condition_y)
      
      # Build DE objects (don’t rely on resultsNames format)
      dds_x <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = samples,
                                              design = stats::as.formula(paste("~", input$metadata_column_x)))
      dds_x[[input$metadata_column_x]] <- stats::relevel(factor(dds_x[[input$metadata_column_x]]),
                                                         ref = input$reference_condition_x)
      dds_x <- DESeq2::DESeq(dds_x)
      
      dds_y <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = samples,
                                              design = stats::as.formula(paste("~", input$metadata_column_y)))
      dds_y[[input$metadata_column_y]] <- stats::relevel(factor(dds_y[[input$metadata_column_y]]),
                                                         ref = input$reference_condition_y)
      dds_y <- DESeq2::DESeq(dds_y)
      
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
      
      # Display label column (SYMBOL if possible)
      merged_df$display_gene <- merged_df$gene
      if (is_ensembl_id(merged_df$gene)) {
        conv <- convert_ensembl_to_symbol(merged_df$gene, filtered_data_rv$species)
        merged_df$display_gene <- ifelse(is.na(conv) | !nzchar(conv), merged_df$gene, conv)
      }
      
      # Labels: top by worst padj & user highlights (match against both SYMBOL and raw ids)
      nlab <- tryCatch({ n <- as.numeric(input$crossplot_topgenes); if (is.na(n) || n < 0) 0 else n }, error = function(e) 0)
      merged_df$combined_padj <- pmax(merged_df$padj_x, merged_df$padj_y, na.rm = TRUE)
      
      top_genes <- character(0)
      if (nlab > 0 && nrow(merged_df)) {
        top_ids <- head(merged_df[order(merged_df$combined_padj), "gene"], nlab)
        # convert to SYMBOL for label if possible
        if (is_ensembl_id(top_ids)) {
          conv <- convert_ensembl_to_symbol(top_ids, filtered_data_rv$species)
          conv[is.na(conv)] <- top_ids[is.na(conv)]
          top_genes <- conv
        } else {
          top_genes <- top_ids
        }
      }
      highlight <- input$crossplot_gene_label
      highlight <- if (nchar(trimws(highlight)) > 0) unlist(strsplit(highlight, "[\\s,]+")) else character(0)
      
      merged_df$label <- ifelse(merged_df$display_gene %in% union(top_genes, highlight), merged_df$display_gene, NA)
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
            max.overlaps = Inf, size = 4, box.padding = 0.15
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
      species <- species %||% "Homo sapiens"
      orgdb <- get_orgdb(species)
      if (is.null(orgdb)) return(NULL)
      
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
      if (is.null(df_ids) || !nrow(df_ids)) return(NULL)
      df_merged <- merge(df, df_ids, by.x = "gene", by.y = 1)
      df_merged <- dplyr::distinct(df_merged, ENTREZID, .keep_all = TRUE)
      df_merged2 <- df_merged[, c("ENTREZID", "log2FoldChange_x", "category")]
      df_merged2 <- df_merged2[stats::complete.cases(df_merged2), , drop = FALSE]
      if (!nrow(df_merged2)) return(NULL)
      
      m <- input$cross_enrich_method
      cc <- NULL
      if (m == "enrichGO") {
        cc <- clusterProfiler::compareCluster(
          ENTREZID ~ category, data = df_merged2, fun = "enrichGO",
          OrgDb = orgdb, keyType = "ENTREZID", ont = "BP",
          pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE
        )
      } else if (m == "groupGO") {
        cc <- clusterProfiler::compareCluster(
          ENTREZID ~ category, data = df_merged2, fun = "groupGO",
          OrgDb = orgdb, keyType = "ENTREZID", ont = "BP", readable = TRUE
        )
      } else if (m == "enrichKEGG") {
        kegg_sp <- get_kegg_code(species); if (is.null(kegg_sp)) return(NULL)
        cc <- clusterProfiler::compareCluster(
          ENTREZID ~ category, data = df_merged2, fun = "enrichKEGG",
          organism = kegg_sp, pvalueCutoff = 0.05, qvalueCutoff = 0.2
        )
      } else if (m == "enrichPathway") {
        rc <- get_reactome_code(species); if (is.null(rc)) return(NULL)
        cc <- clusterProfiler::compareCluster(
          ENTREZID ~ category, data = df_merged2, fun = ReactomePA::enrichPathway,
          organism = rc, pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE
        )
      } else {
        shiny::showNotification("Unsupported enrichment method selected.", type = "error")
        return(NULL)
      }
      
      if (is.null(cc) || !nrow(as.data.frame(cc))) return(NULL)
      # robust dot plot (no fill = p.adjust crash)
      safe_comparecluster_dotplot(cc, showCategory = input$crossplot_topgenes %||% 10, x = "Cluster")
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
      up_x   <- sanitize_vec(df$gene[df$log2FoldChange_x >  lfc_cutoff & df$padj_x < padj_cutoff])
      up_y   <- sanitize_vec(df$gene[df$log2FoldChange_y >  lfc_cutoff & df$padj_y < padj_cutoff])
      down_x <- sanitize_vec(df$gene[df$log2FoldChange_x < -lfc_cutoff & df$padj_x < padj_cutoff])
      down_y <- sanitize_vec(df$gene[df$log2FoldChange_y < -lfc_cutoff & df$padj_y < padj_cutoff])
      
      # If both empty, draw a friendly placeholder instead of throwing
      if (!length(up_x) && !length(up_y)) up_venn  <- grid::grid.text("No upregulated genes")
      else up_venn <- VennDiagram::venn.diagram(
        x = list(X_Up = up_x, Y_Up = up_y),
        filename = NULL, fill = c("darkorange", "darkgreen"), alpha = 0.5,
        main = "Upregulated Genes Venn", disable.logging = TRUE
      )
      if (!length(down_x) && !length(down_y)) down_venn <- grid::grid.text("No downregulated genes")
      else down_venn <- VennDiagram::venn.diagram(
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
        up_x   <- sanitize_vec(df$gene[df$log2FoldChange_x >  lfc_cutoff & df$padj_x < padj_cutoff])
        up_y   <- sanitize_vec(df$gene[df$log2FoldChange_y >  lfc_cutoff & df$padj_y < padj_cutoff])
        down_x <- sanitize_vec(df$gene[df$log2FoldChange_x < -lfc_cutoff & df$padj_x < padj_cutoff])
        down_y <- sanitize_vec(df$gene[df$log2FoldChange_y < -lfc_cutoff & df$padj_y < padj_cutoff])
        
        out <- rbind(
          data.frame(Gene = intersect(up_x, up_y),   Direction = "Upregulated"),
          data.frame(Gene = intersect(down_x, down_y), Direction = "Downregulated")
        )
        utils::write.csv(out, file, row.names = FALSE)
      }
    )
    
    # ---- 6) Heatmap tab (SYMBOL row labels if possible) ----
    generate_cross_heat <- function(df, norm_counts, samples, meta_col_x, species) {
      df$category <- "Other"
      df$category[df$log2FoldChange_x >=  1 & df$log2FoldChange_y >=  1] <- "Up-Up"
      df$category[df$log2FoldChange_x <= -1 & df$log2FoldChange_y <= -1] <- "Down-Down"
      df$category[df$log2FoldChange_x >=  1 & df$log2FoldChange_y <= -1] <- "Up-Down"
      df$category[df$log2FoldChange_x <= -1 & df$log2FoldChange_y >=  1] <- "Down-Up"
      
      sel <- unique(df$gene[!(df$category %in% c("Other", "Comp1-only", "Comp2-only"))])
      if (!length(sel)) return(NULL)
      
      expr_mat <- log2(norm_counts[sel, , drop = FALSE] + 1)
      
      # Rename rows to SYMBOL where possible
      if (is_ensembl_id(sel)) {
        sym <- convert_ensembl_to_symbol(sel, species)
        row_labels <- ifelse(is.na(sym) | !nzchar(sym), sel, sym)
        rownames(expr_mat) <- make.unique(row_labels)
      } else {
        rownames(expr_mat) <- make.unique(rownames(expr_mat))
      }
      
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
