# === Module: mod_gsea_analysis ===

#' GSEA Module
#'
#' Runs Gene Set Enrichment Analysis (GSEA) using clusterProfiler and msigdbr,
#' supports multiple databases (GO, KEGG, Reactome, Hallmark) and visualizes
#' the top enriched pathways with dot plots and enrichment scores.
#'
#' @param input Shiny input object
#' @param output Shiny output object
#' @param session Shiny session object
#' @param filtered_data_rv a list containing species and normalized expression
#' @param res_reactive ReactiveVal containing DE results
#' @import clusterProfiler msigdbr enrichplot
#' @importFrom utils head read.csv write.csv str
#' @importFrom stats as.formula dist model.matrix prcomp quantile relevel var
#' @importFrom grDevices dev.off pdf colorRampPalette
#' @importFrom grid gpar
#' @importFrom shiny isolate req renderPlot renderImage showNotification
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom shinythemes shinytheme
#' @export
mod_gsea_analysis <- function(input, output, session, filtered_data_rv, res_reactive) {

  observeEvent(input$run_gsea, {
    req(res_reactive(), input$gsea_db, input$lfc_threshold, input$padj_threshold, input$gsea_pvalue,filtered_data_rv$species)
    showNotification("Starting GSEA Analysis", type = "message")
    #print(filtered_data_rv$species)
    res <- res_reactive()
    res <- res[!is.na(res$log2FoldChange) & !is.na(res$padj), ]
    res$gene <- rownames(res)
    
    # Convert gene identifiers to symbols
    if (!is_symbol(rownames(res))) {
      conv <- convert_ensembl_to_symbol(rownames(res), filtered_data_rv$species)
      valid_idx <- !is.na(conv)
      res <- res[valid_idx, , drop = FALSE]
      res$gene <- conv[valid_idx]
    } else {
      res$gene <- rownames(res)
    }
    
    # Keep only entries passing thresholds
    res <- res[abs(res$log2FoldChange) >= input$lfc_threshold & res$padj <= input$padj_threshold, ]
    res <- res[!is.na(res$gene), , drop = FALSE]
    res$gene <- make.unique(res$gene)
    
    # Now create the gene list
    geneList <- res$log2FoldChange
    names(geneList) <- res$gene
    geneList <- sort(geneList, decreasing = TRUE)
    geneList <- geneList[!duplicated(names(geneList))]
    
    
    if (length(geneList) < 10) {
      showNotification("Too few genes left after filtering or ID conversion for GSEA.", type = "error")
      return()
    }
    
    species_name <- if (filtered_data_rv$species == "Mus musculus") "Mus musculus" else "Homo sapiens"
    db_collection <- switch(input$gsea_db,
      "GO" = "C5",
      "KEGG" = "C2",
      "Reactome" = "C2",
      "Hallmark" = "H",
      "Cancer Cell Atlas" = "C4",
      "Cancer Gene Neighbourhoods" ="C4",
      "Cancer Modules" = "C4",
      "Txn Factor Targets" = "C3"
    )
    db_subcollection <- switch(input$gsea_db,
      "GO" = "GO:BP",
      "KEGG" = "CP:KEGG_MEDICUS",
      "Reactome" = "CP:REACTOME",
      "Hallmark" = NULL,
      "Cancer Cell Atlas" = "3CA",
      "Cancer Gene Neighbourhoods" ="CGN",
      "Cancer Modules" = "CM",
      "Txn Factor Targets" = "TFT:GTRD"
    )

    m_df <- msigdbr(species = species_name, collection = db_collection, subcollection = db_subcollection)
    term2gene <- m_df[, c("gs_name", "gene_symbol")]
    colnames(term2gene) <- c("ID", "gene")
    #print(head(geneList))
    gsea_result <- tryCatch({
      clusterProfiler::GSEA(
        geneList = geneList,
        TERM2GENE = term2gene,
        pvalueCutoff = input$gsea_pvalue,
        verbose = FALSE
      )
    }, error = function(e) {
      showNotification(paste("GSEA failed:", e$message), type = "error")
      return(NULL)
    })
    
    if (is.null(gsea_result) || nrow(as.data.frame(gsea_result)) < 1) {
      showNotification("No enriched terms found in GSEA or analysis failed.", type = "warning")
      return()
    }
    
    updateSelectInput(session, "gsea_selected_pathway", choices = gsea_result@result$ID)
    output$gseaDotPlot <- renderPlot({
      gsea_result@result$.sign <- ifelse(gsea_result@result$NES > 0, "Activated", "Inhibited")
      color_by <- input$gsea_color_scale

      if (input$gsea_split_dotplot && min(gsea_result@result$NES, na.rm = TRUE) < 0) {
        enrichplot::dotplot(gsea_result, showCategory = input$gsea_top_n, split = ".sign", color = color_by) +
          facet_grid(. ~ .sign) +
          theme(axis.text.y = element_text(size = 6, face = "bold"))
      } else {
        enrichplot::dotplot(gsea_result, showCategory = input$gsea_top_n, color = color_by) +
          theme(axis.text.y = element_text(size = 6, face = "bold"))
      }
    })

    output$download_gsea_dot_plot <- downloadHandler(
      filename = function() { paste0("GSEA_",input$gsea_db,"_dot_plot.pdf", sep="")  },
      content = function(file) {
        pdf(file)
        gsea_result@result$.sign <- ifelse(gsea_result@result$NES > 0, "Activated", "Inhibited")
        color_by <- input$gsea_color_scale

        p <- if (input$gsea_split_dotplot && min(gsea_result@result$NES, na.rm = TRUE) < 0) {
          enrichplot::dotplot(gsea_result, showCategory = input$gsea_top_n, split = ".sign", color = color_by) +
            facet_grid(. ~ .sign) +
            theme(axis.text.y = element_text(size = 6, face = "bold"))
        } else {
          enrichplot::dotplot(gsea_result, showCategory = input$gsea_top_n, color = color_by) +
            theme(axis.text.y = element_text(size = 6, face = "bold"))
        }
        print(p)
        dev.off()
      }
    )

    output$gseaEnrichmentPlot <- renderPlot({
      req(input$gsea_selected_pathway, gsea_result)
      if (!(input$gsea_selected_pathway %in% gsea_result@result$ID)) {
        showNotification("Selected pathway not found in GSEA results.", type = "error")
        return(NULL)
      }
      enrichplot::gseaplot2(gsea_result, geneSetID = input$gsea_selected_pathway)
    })
    
    output$download_gsea_enrichment_plot <- downloadHandler(
      filename = function() {
        paste0("gsea_enrichment_plot_", gsub("\\W+", "_", input$gsea_selected_pathway), ".pdf")
      },
      content = function(file) {
        pdf(file, width = 8, height = 6)
        print(enrichplot::gseaplot2(gsea_result, geneSetID = input$gsea_selected_pathway))
        dev.off()
      }
    )
    
    output$GSEAupsetPlot <- renderPlot({
      #req(pathway_result)
      enrichplot::upsetplot(gsea_result) + theme(axis.text.y = element_text(size = 6, face = "bold"))
    })
    
    output$download_gsea_upset_plot <- downloadHandler(
      filename = function() { paste0("GSEA_",input$gsea_db,"_upset_plot.pdf", sep="") },
      content = function(file) {
        pdf(file)
        print(enrichplot::upsetplot(gsea_result) + theme(axis.text.y = element_text(size = 6, face = "bold")))
        dev.off()
      }
    )
    output$gseaTable <- renderDT({
      datatable(as.data.frame(gsea_result), options = list(scrollX = TRUE))
    })
    output$download_gsea_table <- downloadHandler(
      filename = function() { paste0("GSEA_",input$gsea_db,"_result.csv", sep="")  },
      content = function(file) {
        write.csv(as.data.frame(gsea_result), file, row.names = FALSE)
      }
    )
  })
}



