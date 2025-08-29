# === Module: mod_pca_server (auto PC selection + compareCluster across PCs) ===
#' PCA Cluster Server
#' @param id module id
#' @param filtered_data_rv reactiveValues with $norm_counts, $samples, $species
#' @importFrom shiny moduleServer req observeEvent renderPlot renderText downloadHandler showNotification
#' @importFrom stats prcomp sd cor
#' @importFrom ggplot2 ggplot aes geom_col geom_line scale_y_continuous labs theme_minimal
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation draw
#' @importFrom circlize colorRamp2
#' @importFrom grid grid.newpage grid.text gpar
#' @importFrom grDevices pdf dev.off
#' @export
mod_pca_server <- function(id, filtered_data_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    
    ## ---- tiny local helpers (safe fallbacks) --------------------------------
    is_ensembl_like <- function(x) {
      # vector-safe
      x <- as.character(x)
      if (!length(x)) return(FALSE)
      grepl("^ENS[A-Z]*\\d+", x[1])
    }
    to_symbol_safe <- function(ids, species) {
      ids <- as.character(ids)
      out <- setNames(rep(NA_character_, length(ids)), ids)
      if (exists("convert_ensembl_to_symbol", mode = "function")) {
        conv <- convert_ensembl_to_symbol(ids, species)
        out[names(conv)] <- conv
      }
      out
    }
    get_orgdb_safe <- function(species) {
      if (exists("get_orgdb", mode = "function")) return(get("get_orgdb")(species))
      NULL
    }
    get_kegg_code_safe <- function(species) {
      if (exists("get_kegg_code", mode = "function")) return(get("get_kegg_code")(species))
      NULL
    }
    get_reactome_code_safe <- function(species) {
      if (exists("get_reactome_code", mode = "function")) return(get("get_reactome_code")(species))
      NULL
    }
    
    ## ---- main PCA reactive ---------------------------------------------------
    pca_res <- shiny::eventReactive(input$run_pca, {
      shiny::req(filtered_data_rv$norm_counts, filtered_data_rv$samples)
      
      mat <- as.matrix(filtered_data_rv$norm_counts)
      if (!is.numeric(mat[1, 1])) {
        shiny::showNotification("Counts matrix is not numeric.", type = "error")
        return(NULL)
      }
      
      # Filter by gene-wise SD
      sds  <- apply(mat, 1, stats::sd, na.rm = TRUE)
      keep <- which(sds >= input$variance_threshold)
      if (!length(keep)) {
        shiny::showNotification("No genes pass the variance threshold.", type = "error")
        return(NULL)
      }
      mat_f <- log2(mat[keep, , drop = FALSE] + 1)
      
      # PCA on samples (transpose)
      pr   <- stats::prcomp(t(mat_f), center = TRUE, scale. = TRUE)
      s2   <- pr$sdev^2
      prop <- s2 / sum(s2)
      cprop <- cumsum(prop)
      
      # choose number of PCs to reach coverage
      target <- input$percent_of_data
      max_k  <- min(input$max_pc, length(prop))
      k      <- which(cprop >= target)[1]
      if (is.na(k)) k <- max_k else k <- min(k, max_k)
      
      list(
        prcomp = pr,
        prop = prop,
        cprop = cprop,
        k = k,
        mat_filtered = mat_f
      )
    })
    
    ## ---- variance plot & text ------------------------------------------------
    output$pca_variance_plot <- shiny::renderPlot({
      shiny::req(pca_res())
      pr <- pca_res()
      df <- data.frame(pc = seq_along(pr$prop), var = pr$prop, cum = pr$cprop)
      ggplot2::ggplot(df, ggplot2::aes(x = pc, y = var)) +
        ggplot2::geom_col(fill = "#4C78A8") +
        ggplot2::geom_line(ggplot2::aes(y = cum), color = "#F58518", linewidth = 1) +
        ggplot2::scale_y_continuous(name = "Proportion of variance", limits = c(0, max(pr$cprop))) +
        ggplot2::labs(
          title = "Explained variance per PC (bars) and cumulative (line)",
          x = "Principal component",
          caption = paste0("Selected PCs: 1..", pr$k, " (coverage >= ", round(pr$cprop[pr$k]*100, 1), "%)")
        ) +
        ggplot2::theme_minimal(base_size = 12)
    })
    output$selected_pc_text <- shiny::renderText({
      shiny::req(pca_res())
      pr <- pca_res()
      paste0("Using PCs 1..", pr$k, " covering ", round(pr$cprop[pr$k]*100, 1), "% variance.")
    })
    
    ## ---- helper: pick contributing genes across PCs 1..k --------------------
    top_genes_by_pc <- function(loadings, k, min_genes, max_genes) {
      k <- min(k, ncol(loadings))
      sel <- character(0)
      for (i in seq_len(k)) {
        w <- abs(loadings[, i])
        ord <- order(w, decreasing = TRUE)
        n_take <- max(min_genes, min(max_genes, sum(w > 0)))
        gi <- rownames(loadings)[ord][seq_len(n_take)]
        sel <- union(sel, gi)
      }
      sel
    }
    
    ## ---- helper: build contributing-genes heatmap object ---------------------
    build_contrib_heatmap <- function(pr, species, samples, min_genes, max_genes) {
      load <- pr$prcomp$rotation
      genes <- top_genes_by_pc(load, pr$k, min_genes, max_genes)
      if (!length(genes)) return(NULL)
      
      # expression for selected genes
      expr <- pr$mat_filtered[genes, , drop = FALSE]
      
      # row label conversion: ENSEMBL -> SYMBOL (fallback to ENSEMBL)
      labs <- rownames(expr)
      if (is_ensembl_like(labs)) {
        conv <- to_symbol_safe(labs, species)
        labs[!is.na(conv) & nzchar(conv)] <- conv[!is.na(conv) & nzchar(conv)]
      }
      rownames(expr) <- make.unique(labs)
      
      # z-score by gene
      z <- t(scale(t(expr)))
      z[is.na(z)] <- 0
      
      # top column annotation (cmp_group if present, else first metadata column)
      group_col <- if ("cmp_group" %in% colnames(samples)) "cmp_group" else colnames(samples)[1]
      group_values <- factor(samples[[group_col]])
      ann_df <- data.frame(Group = group_values, check.names = FALSE)
      rownames(ann_df) <- colnames(z)
      lev <- levels(group_values)
      pal <- setNames(scales::hue_pal()(length(lev)), lev)
      top_ha <- ComplexHeatmap::HeatmapAnnotation(df = ann_df, col = list(Group = pal))
      
      hm_col <- circlize::colorRamp2(c(-2, 0, 2), c("#2b8cbe", "white", "#e34a33"))
      ComplexHeatmap::Heatmap(
        z, name = "Z", col = hm_col,
        top_annotation = top_ha,
        show_column_names = FALSE,
        cluster_rows = TRUE, cluster_columns = TRUE,
        row_names_gp = grid::gpar(fontsize = 4, fontface = "bold")
      )
    }
    
    ## ---- contributing genes heatmap (screen) ---------------------------------
    output$pca_loadings_heatmap <- shiny::renderPlot({
      shiny::req(pca_res(), filtered_data_rv$samples, filtered_data_rv$species)
      hp <- build_contrib_heatmap(
        pr       = pca_res(),
        species  = filtered_data_rv$species,
        samples  = filtered_data_rv$samples,
        min_genes = input$min_genes,
        max_genes = input$max_genes
      )
      if (is.null(hp)) { grid::grid.newpage(); grid::grid.text("No contributing genes selected."); return() }
      ComplexHeatmap::draw(hp, heatmap_legend_side = "right", annotation_legend_side = "right")
    })
    
    ## ---- sample correlation heatmap (based on selected PCs) ------------------
    output$pca_sample_heatmap <- shiny::renderPlot({
      shiny::req(pca_res())
      pr <- pca_res()
      scores  <- pr$prcomp$x[, seq_len(pr$k), drop = FALSE]
      cor_mat <- stats::cor(t(scores))
      hm_col  <- circlize::colorRamp2(c(-1, 0, 1), c("#2166ac", "white", "#b2182b"))
      ComplexHeatmap::draw(
        ComplexHeatmap::Heatmap(
          cor_mat, name = "r", col = hm_col, show_row_names = FALSE, show_column_names = FALSE
        ),
        heatmap_legend_side = "right"
      )
    })
    
    ## ---- compareCluster across PCs 1..k (no UI) ------------------------------
    compute_pca_compare_enrichment <- function(pr, species, min_genes, max_genes, method) {
      load <- pr$prcomp$rotation
      k <- min(pr$k, ncol(load))
      pcs <- paste0("PC", seq_len(k))
      
      # per-PC gene lists (top loadings)
      by_pc <- lapply(seq_len(k), function(i) {
        w <- abs(load[, i]); ord <- order(w, decreasing = TRUE)
        n_take <- max(min_genes, min(max_genes, sum(w > 0)))
        rownames(load)[ord][seq_len(n_take)]
      })
      names(by_pc) <- pcs
      
      orgdb <- get_orgdb_safe(species)
      if (is.null(orgdb)) return(NULL)
      
      # map to ENTREZ per PC
      dfe_list <- list()
      for (lab in names(by_pc)) {
        g <- by_pc[[lab]]
        # detect SYMBOL vs ENSEMBL for this PC
        is_sym <- if (exists("is_symbol", mode = "function")) is_symbol(g) else !is_ensembl_like(g)
        map <- tryCatch(
          {
            if (is_sym) {
              clusterProfiler::bitr(g, fromType = "SYMBOL",  toType = "ENTREZID", OrgDb = orgdb)
            } else {
              clusterProfiler::bitr(g, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = orgdb)
            }
          },
          error = function(e) NULL
        )
        if (!is.null(map) && nrow(map)) {
          dfe_list[[lab]] <- data.frame(ENTREZID = unique(map$ENTREZID), Cluster = lab, row.names = NULL)
        }
      }
      if (!length(dfe_list)) return(NULL)
      dfe <- do.call(rbind, dfe_list)
      dfe <- dfe[stats::complete.cases(dfe), , drop = FALSE]
      if (!nrow(dfe)) return(NULL)
      
      # run compareCluster
      if (method == "groupGO") {
        clusterProfiler::compareCluster(
          ENTREZID ~ Cluster, data = dfe, fun = "groupGO",
          OrgDb = orgdb, keyType = "ENTREZID", ont = "BP", readable = TRUE
        )
      } else if (method == "enrichGO") {
        clusterProfiler::compareCluster(
          ENTREZID ~ Cluster, data = dfe, fun = "enrichGO",
          OrgDb = orgdb, keyType = "ENTREZID", ont = "BP",
          pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE
        )
      } else if (method == "enrichKEGG") {
        kc <- get_kegg_code_safe(species); if (is.null(kc)) return(NULL)
        clusterProfiler::compareCluster(
          ENTREZID ~ Cluster, data = dfe, fun = "enrichKEGG",
          organism = kc, pvalueCutoff = 0.05, qvalueCutoff = 0.2
        )
      } else if (method == "enrichPathway") {
        rc <- get_reactome_code_safe(species); if (is.null(rc)) return(NULL)
        clusterProfiler::compareCluster(
          ENTREZID ~ Cluster, data = dfe, fun = ReactomePA::enrichPathway,
          organism = rc, pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE
        )
      } else {
        NULL
      }
    }
    
    ## ---- enrichment dotplot across PCs (screen) ------------------------------
    output$pca_enrichment_plot <- shiny::renderPlot({
      shiny::req(pca_res(), filtered_data_rv$species)
      pr <- pca_res()
      cmp_res <- compute_pca_compare_enrichment(
        pr, filtered_data_rv$species,
        min_genes = input$min_genes,
        max_genes = input$max_genes,
        method    = input$pca_enrich_method
      )
      if (is.null(cmp_res) || !nrow(as.data.frame(cmp_res))) {
        shiny::showNotification("No enriched terms for selected PCs.", type = "warning"); return(NULL)
      }
      # compareCluster groups by "Cluster"
      enrichplot::dotplot(cmp_res, x = "Cluster") +
        ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6, face = "bold"))
    })
    
    ## ---------------------------- Downloads -----------------------------------
    output$download_pca_loadings <- shiny::downloadHandler(
      filename = function() "pca_contributing_genes_heatmap.pdf",
      content = function(file) {
        shiny::req(pca_res())
        grDevices::pdf(file, width = 8.5, height = 11)
        hp <- build_contrib_heatmap(
          pr       = pca_res(),
          species  = filtered_data_rv$species,
          samples  = filtered_data_rv$samples,
          min_genes = input$min_genes,
          max_genes = input$max_genes
        )
        if (is.null(hp)) { grid::grid.newpage(); grid::grid.text("No contributing genes selected.") } else {
          ComplexHeatmap::draw(hp, heatmap_legend_side = "right", annotation_legend_side = "right")
        }
        grDevices::dev.off()
      }
    )
    
    # Reuse contributing genes heatmap as "reconstructed" proxy
    output$download_reconstructed_heatmap <- shiny::downloadHandler(
      filename = function() "pca_reconstructed_heatmap.pdf",
      content = function(file) {
        shiny::req(pca_res())
        grDevices::pdf(file, width = 8.5, height = 11)
        hp <- build_contrib_heatmap(
          pr       = pca_res(),
          species  = filtered_data_rv$species,
          samples  = filtered_data_rv$samples,
          min_genes = input$min_genes,
          max_genes = input$max_genes
        )
        if (is.null(hp)) { grid::grid.newpage(); grid::grid.text("No contributing genes selected.") } else {
          ComplexHeatmap::draw(hp, heatmap_legend_side = "right", annotation_legend_side = "right")
        }
        grDevices::dev.off()
      }
    )
    
    output$download_sample_correlation_heatmap <- shiny::downloadHandler(
      filename = function() "pca_sample_correlation_heatmap.pdf",
      content = function(file) {
        shiny::req(pca_res())
        grDevices::pdf(file, width = 8, height = 8)
        pr <- pca_res()
        scores  <- pr$prcomp$x[, seq_len(pr$k), drop = FALSE]
        cor_mat <- stats::cor(t(scores))
        hm_col  <- circlize::colorRamp2(c(-1, 0, 1), c("#2166ac", "white", "#b2182b"))
        ComplexHeatmap::draw(
          ComplexHeatmap::Heatmap(
            cor_mat, name = "r", col = hm_col, show_row_names = FALSE, show_column_names = FALSE
          ),
          heatmap_legend_side = "right"
        )
        grDevices::dev.off()
      }
    )
    
    output$download_pca_enrichment <- shiny::downloadHandler(
      filename = function() "pca_compareCluster_dotplot.pdf",
      content = function(file) {
        shiny::req(pca_res(), filtered_data_rv$species)
        grDevices::pdf(file, width = 10, height = 8)
        pr <- pca_res()
        cmp_res <- compute_pca_compare_enrichment(
          pr, filtered_data_rv$species,
          min_genes = input$min_genes,
          max_genes = input$max_genes,
          method    = input$pca_enrich_method
        )
        if (!is.null(cmp_res) && nrow(as.data.frame(cmp_res))) {
          print(enrichplot::dotplot(cmp_res, x = "Cluster") +
                  ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6, face = "bold")))
        } else {
          grid::grid.newpage(); grid::grid.text("No enriched terms for selected PCs.")
        }
        grDevices::dev.off()
      }
    )
    
    output$download_pca_loadings_table <- shiny::downloadHandler(
      filename = function() "pca_contributing_genes.csv",
      content = function(file) {
        shiny::req(pca_res())
        pr <- pca_res()
        load <- pr$prcomp$rotation
        df <- as.data.frame(load[, seq_len(pr$k), drop = FALSE])
        df <- df[order(rowSums(abs(df)), decreasing = TRUE), , drop = FALSE]
        utils::write.csv(df, file, row.names = TRUE)
      }
    )
    
    output$download_pca_enrichment_table <- shiny::downloadHandler(
      filename = function() "pca_compare_cluster_table.csv",
      content = function(file) {
        shiny::req(pca_res(), filtered_data_rv$species)
        pr <- pca_res()
        cmp_res <- compute_pca_compare_enrichment(
          pr, filtered_data_rv$species,
          min_genes = input$min_genes,
          max_genes = input$max_genes,
          method    = input$pca_enrich_method
        )
        if (is.null(cmp_res)) {
          utils::write.csv(data.frame(), file, row.names = FALSE)
        } else {
          utils::write.csv(as.data.frame(cmp_res), file, row.names = FALSE)
        }
      }
    )
  })
}
