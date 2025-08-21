#' GSEA Module (server)
#'
#' Runs Gene Set Enrichment Analysis (GSEA) using clusterProfiler and msigdbr,
#' supports multiple databases (GO, KEGG, Reactome, Hallmark) and visualizes
#' the top enriched pathways with dot plots and enrichment scores.
#'
#' @param id module id (must match mod_gsea_tab_ui)
#' @param filtered_data_rv reactiveValues list containing at least $species
#' @param res_reactive reactiveVal containing DE results (data.frame)
#' @param de_sel (optional) list returned by mod_de_server(), providing ref/test levels for filenames
#' @import clusterProfiler msigdbr enrichplot
#' @importFrom utils head write.csv
#' @importFrom stats as.formula dist model.matrix prcomp quantile relevel var
#' @importFrom grDevices dev.off pdf
#' @importFrom grid gpar
#' @importFrom shiny req renderPlot showNotification downloadHandler
#' @export
mod_gsea_server <- function(id, filtered_data_rv, res_reactive, de_sel = NULL) {
  moduleServer(id, function(input, output, session) {
    
    observeEvent(input$run_gsea, {
      shiny::req(
        res_reactive(),
        input$gsea_db,
        input$lfc_threshold,
        input$padj_threshold,
        input$gsea_pvalue,
        filtered_data_rv$species
      )
      shiny::showNotification("Starting GSEA Analysis", type = "message")
      
      # Prepare DE results
      res <- res_reactive()
      res <- res[!is.na(res$log2FoldChange) & !is.na(res$padj), , drop = FALSE]
      res$gene <- rownames(res)
      
      # ID handling (SYMBOL vs ENSEMBL); use your app helpers
      if (!is_symbol(rownames(res))) {
        conv <- convert_ensembl_to_symbol(rownames(res), filtered_data_rv$species)
        valid_idx <- !is.na(conv)
        res <- res[valid_idx, , drop = FALSE]
        res$gene <- conv[valid_idx]
      } else {
        res$gene <- rownames(res)
      }
      
      # Filter and rank vector
      res <- res[abs(res$log2FoldChange) >= input$lfc_threshold & res$padj <= input$padj_threshold, , drop = FALSE]
      res <- res[!is.na(res$gene), , drop = FALSE]
      res$gene <- make.unique(res$gene)
      
      geneList <- res$log2FoldChange
      names(geneList) <- res$gene
      geneList <- sort(geneList, decreasing = TRUE)
      geneList <- geneList[!duplicated(names(geneList))]
      
      if (length(geneList) < 10) {
        shiny::showNotification("Too few genes left after filtering or ID conversion for GSEA.", type = "error")
        return()
      }
      
      # Species mapping for msigdbr
      species_map <- list(
        "Homo sapiens"             = "Homo sapiens",
        "Mus musculus"             = "Mus musculus",
        "Rattus norvegicus"        = "Rattus norvegicus",
        "Canis familiaris"         = "dog",
        "Saccharomyces cerevisiae" = "Saccharomyces cerevisiae"
      )
      species_name <- species_map[[filtered_data_rv$species]]
      if (is.null(species_name)) {
        shiny::showNotification("Selected species not supported by msigdbr.", type = "error")
        return()
      }
      
      # DB -> collection/subcollection
      db_collection <- switch(
        input$gsea_db,
        "GO"                         = "C5",
        "KEGG"                       = "C2",
        "Reactome"                   = "C2",
        "Hallmark"                   = "H",
        "Cancer Cell Atlas"          = "C4",
        "Cancer Gene Neighbourhoods" = "C4",
        "Cancer Modules"             = "C4",
        "Txn Factor Targets"         = "C3"
      )
      db_subcollection <- switch(
        input$gsea_db,
        "GO"                         = "GO:BP",
        "KEGG"                       = "CP:KEGG_MEDICUS",
        "Reactome"                   = "CP:REACTOME",
        "Hallmark"                   = NULL,
        "Cancer Cell Atlas"          = "3CA",
        "Cancer Gene Neighbourhoods" = "CGN",
        "Cancer Modules"             = "CM",
        "Txn Factor Targets"         = "TFT:GTRD"
      )
      
      m_df <- msigdbr::msigdbr(
        species = species_name,
        collection = db_collection,
        subcollection = db_subcollection
      )
      term2gene <- m_df[, c("gs_name", "gene_symbol")]
      colnames(term2gene) <- c("ID", "gene")
      
      # Run GSEA
      gsea_result <- tryCatch({
        clusterProfiler::GSEA(
          geneList = geneList,
          TERM2GENE = term2gene,
          pvalueCutoff = input$gsea_pvalue,
          verbose = FALSE
        )
      }, error = function(e) {
        shiny::showNotification(paste("GSEA failed:", e$message), type = "error")
        NULL
      })
      
      if (is.null(gsea_result) || nrow(as.data.frame(gsea_result)) < 1) {
        shiny::showNotification("No enriched terms found in GSEA or analysis failed.", type = "warning")
        return()
      }
      
      updateSelectInput(session, "gsea_selected_pathway", choices = gsea_result@result$ID)
      
      # Dot plot (screen)
      output$gseaDotPlot <- shiny::renderPlot({
        shiny::req(gsea_result)
        df <- gsea_result@result
        df$.sign <- ifelse(df$NES > 0, "Activated", "Inhibited")
        color_by <- input$gsea_color_scale
        top_n <- input$gsea_top_n
        df <- head(df[order(df$p.adjust), ], top_n)
        
        if (identical(color_by, "NES")) {
          df$coreSize  <- sapply(strsplit(df$core_enrichment, "/"), length)
          df$geneRatio <- df$coreSize / df$setSize
          p <- ggplot2::ggplot(df, ggplot2::aes(x = geneRatio, y = reorder(Description, geneRatio),
                                                size = coreSize, color = NES)) +
            ggplot2::geom_point() +
            ggplot2::scale_color_gradient2(low = "blue", high = "red", midpoint = 0) +
            ggplot2::labs(x = "Gene Ratio (coreSize / setSize)", y = NULL, color = "NES", size = "Gene Count") +
            ggplot2::theme_minimal() +
            ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6, face = "bold"))
          if (isTRUE(input$gsea_split_dotplot) && any(df$NES < 0, na.rm = TRUE)) {
            p <- p + ggplot2::facet_wrap(~.sign)
          }
          p
        } else {
          if (isTRUE(input$gsea_split_dotplot) && any(df$NES < 0, na.rm = TRUE)) {
            enrichplot::dotplot(gsea_result, showCategory = top_n, split = ".sign", color = color_by) +
              ggplot2::facet_grid(. ~ .sign) +
              ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6, face = "bold"))
          } else {
            enrichplot::dotplot(gsea_result, showCategory = top_n, color = color_by) +
              ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6, face = "bold"))
          }
        }
      })
      tag <- .contrast_tag_from(de_sel)
      # Dot plot (PDF)
      output$download_gsea_dot_plot <- shiny::downloadHandler(
        filename = function() paste0("GSEA_", input$gsea_db,  "_", tag, "_dot_plot.pdf"),
        content  = function(file) {
          grDevices::pdf(file)
          print(enrichplot::dotplot(gsea_result, showCategory = input$gsea_top_n,
                                    color = input$gsea_color_scale))
          grDevices::dev.off()
        }
      )
      
      # Enrichment plot (screen + PDF)
      output$gseaEnrichmentPlot <- shiny::renderPlot({
        shiny::req(input$gsea_selected_pathway, gsea_result)
        if (!(input$gsea_selected_pathway %in% gsea_result@result$ID)) {
          shiny::showNotification("Selected pathway not found in GSEA results.", type = "error")
          return(NULL)
        }
        enrichplot::gseaplot2(gsea_result, geneSetID = input$gsea_selected_pathway)
      })
      output$download_gsea_enrichment_plot <- shiny::downloadHandler(
        filename = function() paste0("gsea_enrichment_plot_",
                                     gsub("\\W+", "_", input$gsea_selected_pathway), "_", tag,  ".pdf"),
        content  = function(file) {
          grDevices::pdf(file, width = 8, height = 6)
          print(enrichplot::gseaplot2(gsea_result, geneSetID = input$gsea_selected_pathway))
          grDevices::dev.off()
        }
      )
      
      # Upset plot (screen + PDF)
      output$GSEAupsetPlot <- shiny::renderPlot({
        enrichplot::upsetplot(gsea_result) +
          ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6, face = "bold"))
      })
      output$download_gsea_upset_plot <- shiny::downloadHandler(
        filename = function() paste0("GSEA_", input$gsea_db, "_", tag,  "_upset_plot.pdf"),
        content  = function(file) {
          grDevices::pdf(file)
          print(enrichplot::upsetplot(gsea_result) +
                  ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6, face = "bold")))
          grDevices::dev.off()
        }
      )
      
      # Table (screen + CSV)
      output$gseaTable <- DT::renderDT({
        DT::datatable(as.data.frame(gsea_result), options = list(scrollX = TRUE))
      })
      output$download_gsea_table <- shiny::downloadHandler(
        filename = function() paste0("GSEA_", input$gsea_db, "_", tag,  "_result.csv"),
        content  = function(file) utils::write.csv(as.data.frame(gsea_result), file, row.names = FALSE)
      )
    })
  })
}
