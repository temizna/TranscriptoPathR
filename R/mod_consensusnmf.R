#' Consensus Clustering and Classification Module Server with NMF (Enhanced Styling + coefmap + marker_heatmap fixed)
#'
#' @param input Shiny input object
#' @param output Shiny output object
#' @param session Shiny session object
#' @param filtered_data_rv A reactiveValues object containing:
#'   \itemize{
#'     \item \code{norm_counts}: normalized gene expression matrix (genes x samples)
#'     \item \code{species}: species name (e.g., "Homo sapiens")
#'     \item \code{metadata}: data frame with sample metadata (columns annotated on heatmap)
#'   }
#'
#' @return None (called for its side effects to render plots and outputs in Shiny)
#' @importFrom ComplexHeatmap Heatmap draw HeatmapAnnotation
#' @importFrom clusterProfiler compareCluster enrichGO enrichKEGG enrichDO groupGO bitr
#' @importFrom ReactomePA enrichPathway
#' @importFrom enrichplot dotplot
#' @importFrom DT renderDT datatable
#' @export
mod_pca_cluster_server <- function(input, output, session, filtered_data_rv) {
  gene_pc_map <- reactiveVal(NULL)
  pca_res <- reactiveVal(NULL)
  enrichment_res <- reactiveVal(NULL)
  
  observeEvent(input$run_pca, {
    req(filtered_data_rv$norm_counts)
    mat <- as.matrix(filtered_data_rv$norm_counts)
    
    # Filter by variance threshold
    gene_vars <- apply(mat, 1, var, na.rm = TRUE)
    mat <- mat[gene_vars >= input$variance_threshold, , drop = FALSE]
    
    # Median-center and log2 transform
    mat_centered <- t(scale(t(mat), center = apply(mat, 1, median, na.rm = TRUE), scale = FALSE))
    mat_log2 <- log2(mat_centered + 1)
    mat_log2 <- mat_log2[apply(mat_log2, 1, function(x) all(is.finite(x))), , drop = FALSE]
    
    # PCA
    pca_obj <- prcomp(t(mat_log2), center = FALSE, scale. = FALSE)
    pca_res(pca_obj)
    
    # CDF of variance explained
    var_explained <- summary(pca_obj)$importance["Proportion of Variance", 1:input$max_pc]
    cdf_var <- cumsum(var_explained)
    
    output$pca_variance_plot <- renderPlot({
      plot(cdf_var, type = "l", lwd = 2, col = "darkblue", xlab = "Principal Components", 
           ylab = "Cumulative Variance Explained", ylim=c(0,1))
      abline(h = input$variance_threshold, col = "red", lty = 2)
    })
    
    output$download_pca_variance_plot <- downloadHandler(
      filename = function() paste0("PCA_variance_CDF_", Sys.Date(), ".pdf"),
      content = function(file) {
        pdf(file)
        plot(cdf_var, type = "l", lwd = 2, col = "darkblue", xlab = "Principal Components", 
             ylab = "Cumulative Variance Explained", ylim=c(0,1))
        abline(h = input$variance_threshold, col = "red", lty = 2)
        dev.off()
      }
    )
    
    # Select number of PCs
    num_pcs <- max(which(cdf_var <= input$percent_of_data))  
    if (is.na(num_pcs) || num_pcs < 1) num_pcs <- 1
    if (num_pcs > input$max_pc) num_pcs <- input$max_pc
    
    output$selected_pc_text <- renderText({ paste("Selected number of PCs:", num_pcs) })
    
    # Gene selection per PC
    loadings <- pca_obj$rotation[, 1:num_pcs, drop = FALSE]
    abs_loadings <- abs(loadings)
    
    top_genes_per_pc <- lapply(1:num_pcs, function(pc_idx) {
      sorted <- sort(abs_loadings[, pc_idx], decreasing = TRUE)
      threshold <- sorted[min(input$max_genes, length(sorted))]
      selected <- names(sorted[sorted >= threshold])
      if (length(selected) < input$min_genes) warning(paste("PC", pc_idx, "has < min genes"))
      selected
    })
    names(top_genes_per_pc) <- paste0("PC", 1:num_pcs)
    
    # Include overlapping genes
    gene_pc_table <- lapply(rownames(abs_loadings), function(gene) {
      max_loading <- max(abs_loadings[gene, ])
      pcs <- which(abs_loadings[gene, 1:num_pcs] >= (1 - input$similarity_threshold) * max_loading)
      data.frame(gene = gene, pc = paste0("PC", pcs))
    })
    gene_pc_table <- do.call(rbind, gene_pc_table)
    gene_pc_table <- gene_pc_table[gene_pc_table$gene %in% unlist(top_genes_per_pc), ]
    gene_pc_map(gene_pc_table)
    
    # Loadings heatmap
    output$pca_loadings_heatmap <- renderPlot({
      selected_genes <- unique(gene_pc_table$gene)
      ComplexHeatmap::Heatmap(loadings[selected_genes, , drop = FALSE], name = "PC Loadings",
                              cluster_columns = TRUE, cluster_rows = TRUE,
                              row_names_gp = grid::gpar(fontsize = 4),
                              column_names_gp = grid::gpar(fontsize = 4))
    })
    output$download_pca_loadings <- downloadHandler(
      filename = function() paste0("PCA_Loadings_Heatmap_", Sys.Date(), ".pdf"),
      content = function(file) {
        pdf(file)
        selected_genes <- unique(gene_pc_map()$gene)
        loadings <- pca_res()$rotation[, 1:num_pcs, drop = FALSE]
        mat_selected <- loadings[selected_genes, , drop = FALSE]
        ComplexHeatmap::draw(
          ComplexHeatmap::Heatmap(mat_selected,
                                  name = "PC Loadings",
                                  cluster_rows = TRUE,
                                  cluster_columns = TRUE,
                                  column_names_gp = grid::gpar(fontsize = 4),
                                  row_names_gp = grid::gpar(fontsize = 4))
        )
        dev.off()
      }
    )
    output$download_pca_loadings_table <- downloadHandler(
      filename = function() paste0("PCA_Loadings_GMT_", Sys.Date(), ".gmt"),
      content = function(file) {
        req(gene_pc_map())
        gene_pc_df <- gene_pc_map()
        # Safeguard
        if (is.null(gene_pc_df) || nrow(gene_pc_df) == 0 || !all(c("gene", "pc") %in% colnames(gene_pc_df))) {
          writeLines("No data available for PCA loadings export.", con = file)
          return()
        }
        
        # Ensure 'pc' is a character for GMT format
        gene_pc_df$pc <- as.character(gene_pc_df$pc)
        pcs <- unique(gene_pc_df$pc)
        
        lines <- lapply(pcs, function(pc) {
          genes <- gene_pc_df$gene[gene_pc_df$pc == pc]
          paste(c(pc, "PCA gene set", genes), collapse = "\t")
        })
        
        writeLines(lines, con = file)
      }
    )
    
    
    # Sample correlation heatmap using selected PCs
    scores <- pca_obj$x[, 1:num_pcs, drop = FALSE]
    sample_corr <- cor(t(scores))
    output$pca_sample_heatmap <- renderPlot({
      ht <- ComplexHeatmap::Heatmap(
        sample_corr,
        name = "Sample Correlation",
        cluster_rows = TRUE, cluster_columns = TRUE,
        column_names_gp = grid::gpar(fontsize = 5),
        row_names_gp = grid::gpar(fontsize = 5)
      )
      ComplexHeatmap::draw(ht)
    }, height = 800, width = 800)
    
    
    # Reconstructed matrix using selected PCs
    reconstruction <- pca_obj$x[, 1:num_pcs, drop = FALSE] %*% t(pca_obj$rotation[, 1:num_pcs, drop = FALSE])
    reconstruction <- t(reconstruction) 
    rownames(reconstruction) <- rownames(mat_log2)
    colnames(reconstruction) <- colnames(mat_log2)
    
    output$pca_reconstructed_heatmap <- renderPlot({
      ComplexHeatmap::Heatmap(reconstruction,
                              name = "Reconstructed",
                              col = circlize::colorRamp2(c(min(reconstruction), mean(reconstruction), max(reconstruction)), c("blue", "white", "yellow")),
                              show_row_names = FALSE,
                              column_names_gp = grid::gpar(fontsize = 4),
                              cluster_columns = TRUE,
                              cluster_rows = TRUE)
    })
    output$download_sample_correlation_heatmap <- downloadHandler(
      filename = function() paste0("PCA_sample_correlation_heatmap_", Sys.Date(), ".pdf"),
      content = function(file) {
        pdf(file, width = 11, height = 14)  # larger canvas
        ht <- ComplexHeatmap::Heatmap(
          sample_corr,
          name = "Sample Corr",
          cluster_rows = TRUE,
          cluster_columns = TRUE,
          column_names_gp = grid::gpar(fontsize = 5),
          row_names_gp = grid::gpar(fontsize = 5)
        )
        ComplexHeatmap::draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
        dev.off()
      }
    )
    
    output$download_reconstructed_heatmap <- downloadHandler(
      filename = function() paste0("PCA_reconstructed_heatmap_", Sys.Date(), ".pdf"),
      content = function(file) {
        pdf(file)
        reconstructed <- pca_obj$x[, 1:num_pcs, drop = FALSE] %*% 
          t(pca_obj$rotation[, 1:num_pcs, drop = FALSE])
        reconstructed <- t(reconstructed)  # genes x samples
        rownames(reconstructed) <- rownames(mat_log2)
        colnames(reconstructed) <- colnames(mat_log2)
        ComplexHeatmap::draw(
          ComplexHeatmap::Heatmap(
            reconstructed,
            name = "Reconstructed",
            cluster_rows = TRUE,
            cluster_columns = TRUE,
            col = circlize::colorRamp2(c(min(reconstructed), 0, max(reconstructed)), c("blue", "white", "yellow")),
            column_names_gp = grid::gpar(fontsize = 4),
            show_row_names = FALSE
          )
        )
        dev.off()
      }
    )
    
    
    # Enrichment
    species <- filtered_data_rv$species
    orgdb <- get_orgdb(species)
    genes <- unique(gene_pc_table$gene)
    
    if (!is_symbol(genes)) {
      converted <- convert_ensembl_to_symbol(genes, species = species)
      genes <- ifelse(!is.na(converted), converted, genes)
    }
    
    entrez <- suppressMessages({
      clusterProfiler::bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = orgdb)
    })
    
    term2gene <- dplyr::left_join(gene_pc_table, entrez, by = c("gene" = "SYMBOL")) %>%
      dplyr::select(pc, ENTREZID) %>%
      dplyr::filter(!is.na(ENTREZID))
    
    enrich_fun <- switch(input$pca_enrich_method,
                         "groupGO" = if (species %in% c("Homo sapiens", "Mus musculus", "Rattus norvegicus", "Canis familiaris", "Saccharomyces cerevisiae")) {
                           function(...) enrichGO(OrgDb = orgdb, ...)
                         } else {
                           showNotification("groupGO is not supported for this species", type = "error")
                           return(NULL)
                         },
                         "enrichGO" = if (species %in% c("Homo sapiens", "Mus musculus", "Rattus norvegicus", "Canis familiaris", "Saccharomyces cerevisiae")) {
                           function(...) enrichGO(OrgDb = orgdb, ...)
                         } else {
                           showNotification("enrichGO is not supported for this species", type = "error")
                           return(NULL)
                         },
                         "enrichKEGG" = {
                           kegg_code <- tryCatch(get_kegg_code(species), error = function(e) NULL)
                           if (is.null(kegg_code)) {
                             showNotification("KEGG enrichment is not supported for this species.", type = "error")
                             return(NULL)
                           }
                           function(...) enrichKEGG(organism = kegg_code, ...)
                         },
                         "enrichDO" = if (species == "Homo sapiens") {
                           enrichDO
                         } else {
                           showNotification("DO enrichment is only available for Homo sapiens.", type = "error")
                           return(NULL)
                         },
                         "enrichPathway" = if (species %in% c("Homo sapiens", "Mus musculus")) {
                           function(...) enrichPathway(organism = get_reactome_code(species), ...)
                         } else {
                           showNotification("Reactome enrichment is only supported for human and mouse.", type = "error")
                           return(NULL)
                         }
    )
    
    
    res <- tryCatch({
      compareCluster(ENTREZID ~ pc, data = term2gene, fun = enrich_fun)
    }, error = function(e) {
      showNotification("Pathway enrichment failed.", type = "error")
      NULL
    })
    enrichment_res(res)
    
    
    output$pca_enrichment_plot <- renderPlot({
      req(enrichment_res())
      if (inherits(enrichment_res(), "compareClusterResult") && nrow(as.data.frame(enrichment_res())) > 0) {
        enrichplot::dotplot(enrichment_res(), x = "pc") + ggplot2::theme_minimal()
      } else {
        showNotification("No enrichment results to display.", type = "warning")
        NULL
      }
    })
    
    output$download_pca_enrichment <- downloadHandler(
      filename = function() paste0("PCA_enrichment_", input$pca_enrich_method, "_", Sys.Date(), ".pdf"),
      content = function(file) {
        pdf(file)
        enrichplot::dotplot(enrichment_res(), x = "pc") + ggplot2::theme_minimal()
        dev.off()
      }
    )
    
    output$download_pca_enrichment_table <- downloadHandler(
      filename = function() paste0("PCA_enrichment_table_", input$pca_enrich_method, "_", Sys.Date(), ".csv"),
      content = function(file) {
        write.csv(as.data.frame(enrichment_res()), file, row.names = FALSE)
      }
    )
  })
}
