#' GSEA Module (server)
#'
#' Runs Gene Set Enrichment Analysis (GSEA) using clusterProfiler and msigdbr,
#' supports multiple databases (GO, KEGG, Reactome, Hallmark) and visualizes
#' the top enriched pathways with dot plots and enrichment scores.
#'
#' @param id module id (must match mod_gsea_tab_ui)
#' @param filtered_data_rv reactiveValues list containing at least $species
#' @param res_reactive reactiveVal containing DE results (data.frame)
#' @param de_sel optional list returned by mod_de_server(), used for filenames if cmp is NULL
#' @param cmp optional comparison bridge from make_cmp_bridge(), used for filenames if present
#'
#' @import clusterProfiler msigdbr enrichplot
#' @importFrom DT renderDT datatable
#' @importFrom utils write.csv
#' @importFrom grDevices dev.off pdf
#' @importFrom shiny req renderPlot showNotification downloadHandler updateSelectInput moduleServer
#' @export
mod_gsea_server <- function(id, filtered_data_rv, res_reactive, de_sel = NULL, cmp = NULL) {
  moduleServer(id, function(input, output, session) {
      
      observeEvent(input$run_gsea, {
        shiny::req(
          input$gsea_db,
          input$lfc_threshold,
          input$padj_threshold,
          input$gsea_pvalue,
          filtered_data_rv$species
        )
        shiny::showNotification("Starting GSEA Analysis", type = "message")
        
        use_custom <- isTRUE(input$gsea_use_custom_file)
        res <- NULL
        
        if (use_custom) {
          file <- input$gsea_custom_file
          if (is.null(file)) {
            shiny::showNotification("No file uploaded.", type = "error")
            return()
          }
          
          ext <- tools::file_ext(file$datapath)
          res <- tryCatch({
            if (ext %in% c("csv")) {
              read.csv(file$datapath, stringsAsFactors = FALSE)
            } else {
              read.delim(file$datapath, stringsAsFactors = FALSE)
            }
          }, error = function(e) {
            shiny::showNotification(paste("File read failed:", e$message), type = "error")
            return(NULL)
          })
          
          if (is.null(res)) return()
          
          required_cols <- c("gene", "log2FoldChange", "padj")
          if (!all(required_cols %in% colnames(res))) {
            shiny::showNotification("Uploaded file must contain 'gene', 'log2FoldChange', and 'padj' columns.", type = "error")
            return()
          }
          
          res <- res[!is.na(res$log2FoldChange) & !is.na(res$padj), ]
          res$gene <- make.unique(as.character(res$gene))
          rownames(res) <- res$gene
          
          shiny::showNotification("Using custom DE results", type = "message")
          
        } else {
          res <- res_reactive()
          shiny::req(res)
          shiny::showNotification("Using default DESeq2 results", type = "message")
        }
        
        req(res)
        res <- res[!is.na(res$log2FoldChange) & !is.na(res$padj), , drop = FALSE]
        
        if (!"gene" %in% colnames(res)) {
          res$gene <- rownames(res)
        }
        
        # ID handling
        if (!is_symbol(res$gene)) {
          conv <- convert_ensembl_to_symbol(res$gene, filtered_data_rv$species)
          keep <- !is.na(conv)
          res <- res[keep, , drop = FALSE]
          res$gene <- conv[keep]
        }
        
        res$gene <- make.unique(res$gene)
        res <- res[abs(res$log2FoldChange) >= input$lfc_threshold & res$padj <= input$padj_threshold, ]
        res <- res[!duplicated(res$gene), ]
        
        geneList <- res$log2FoldChange
        names(geneList) <- res$gene
        geneList <- sort(geneList, decreasing = TRUE)
        
        if (length(geneList) < 10) {
          shiny::showNotification("Too few genes for GSEA after filtering.", type = "error")
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
        species     = species_name,
        collection  = db_collection,
        subcollection = db_subcollection
      )
      if (is.null(m_df) || nrow(m_df) == 0) {
        shiny::showNotification("No gene sets returned from msigdbr for this selection.", type = "warning")
        return()
      }
      term2gene <- m_df[, c("gs_name", "gene_symbol")]
      colnames(term2gene) <- c("ID", "gene")
      
      # Run GSEA
      gsea_result <- tryCatch({
        clusterProfiler::GSEA(
          geneList     = geneList,
          TERM2GENE    = term2gene,
          pvalueCutoff = input$gsea_pvalue,
          verbose      = FALSE
        )
      }, error = function(e) {
        shiny::showNotification(paste("GSEA failed:", e$message), type = "error")
        NULL
      })
      
      if (is.null(gsea_result) || nrow(as.data.frame(gsea_result)) < 1) {
        shiny::showNotification("No enriched terms found in GSEA or analysis failed.", type = "warning")
        return()
      }
      
      shiny::updateSelectInput(session, "gsea_selected_pathway", choices = gsea_result@result$ID)
      
      # Determine filename tag from cmp (preferred) or de_sel fallback
      tag <- tryCatch({
        if (!is.null(cmp)) cmp$tag() else .contrast_tag_from(de_sel)
      }, error = function(e) .contrast_tag_from(de_sel))
      
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
      
      # Dot plot (PDF)
      output$download_gsea_dot_plot <- shiny::downloadHandler(
        filename = function() paste0("GSEA_", input$gsea_db, "_", tag, "_dot_plot.pdf"),
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
                                     gsub("\\W+", "_", input$gsea_selected_pathway), "_", tag, ".pdf"),
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
        filename = function() paste0("GSEA_", input$gsea_db, "_", tag, "_upset_plot.pdf"),
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
        filename = function() paste0("GSEA_", input$gsea_db, "_", tag, "_result.csv"),
        content  = function(file) utils::write.csv(as.data.frame(gsea_result), file, row.names = FALSE)
      )
    })
  })
}
