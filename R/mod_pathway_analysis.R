#' Pathway Analysis Server (namespaced, full logic)
#'
#' @param id module id (must match mod_pathway_tab_ui)
#' @param filtered_data_rv reactiveValues with $counts, $norm_counts, $samples, $species
#' @param res_reactive reactiveVal holding DESeq2 results (data.frame)
#' @param geneList_rv reactiveVal to store named log2FC vector (ENTREZID -> log2FC)
#' @param kegg_pathway_results reactiveVal (kept for compatibility; not required here)
#' @param d1_merged_rv reactiveVal to store merged df (gene, log2FC, padj, ENTREZID)
#' @param pathway_result_rv reactiveVal to store the enrichment result object
##' @param de_sel (optional) list returned by mod_de_server(), providing ref/test levels for filenames
#' @export
mod_pathway_server <- function(
    id,
    filtered_data_rv,
    res_reactive,
    geneList_rv,
    kegg_pathway_results,
    d1_merged_rv,
    de_sel = NULL,
    pathway_result_rv
) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # ---- Run enrichment ------------------------------------------------------
    observeEvent(input$run_pathway, {
      shiny::req(filtered_data_rv$species, res_reactive())
      shiny::showNotification("Starting Pathway Analysis", type = "message")
      
      species   <- filtered_data_rv$species
      orgdb     <- .get_orgdb(species)
      if (is.null(orgdb)) {
        shiny::showNotification("Could not load OrgDb for selected species.", type = "error")
        return()
      }
      
      res <- res_reactive()
      res <- res[!is.na(res$log2FoldChange) & !is.na(res$padj), ]
      d1  <- res[, c("log2FoldChange", "padj"), drop = FALSE]
      d1$gene <- rownames(res)
      
      # Gene ID mapping -> ENTREZID
      from_type <- if (.is_symbol(d1$gene)) "SYMBOL" else "ENSEMBL"
      d1_ids <- suppressMessages({
        clusterProfiler::bitr(d1$gene, fromType = from_type, toType = "ENTREZID", OrgDb = orgdb)
      })
      if (is.null(d1_ids) || nrow(d1_ids) == 0) {
        shiny::showNotification("No genes could be mapped to ENTREZ IDs.", type = "error"); return()
      }
      
      d1_merged <- merge(d1, d1_ids, by.x = "gene", by.y = 1)
      d1_merged <- d1_merged[!duplicated(d1_merged$ENTREZID), ]
      d1_merged_rv(d1_merged)
      
      # Direction subset
      dir <- input$pathway_direction
      keep <- switch(
        dir,
        "Up"   = d1_merged$log2FoldChange >=  input$lfc_threshold & d1_merged$padj <= input$padj_threshold,
        "Down" = d1_merged$log2FoldChange <= -input$lfc_threshold & d1_merged$padj <= input$padj_threshold,
        # Both
        abs(d1_merged$log2FoldChange) >= input$lfc_threshold & d1_merged$padj <= input$padj_threshold
      )
      gene_vector <- d1_merged[keep & !is.na(d1_merged$ENTREZID), , drop = FALSE]
      
      geneList <- gene_vector$log2FoldChange
      names(geneList) <- gene_vector$ENTREZID
      geneList <- geneList[!duplicated(names(geneList))]
      geneList <- sort(geneList, decreasing = TRUE)
      
      max_genes <- if (!is.null(input$max_genes)) input$max_genes else 1000
      if (length(geneList) > max_genes) geneList <- head(geneList, max_genes)
      
      if (length(geneList) < 10) {
        shiny::showNotification("Too few mapped genes for pathway analysis.", type = "error")
        return()
      }
      geneList_rv(geneList)
      selected_genes <- names(geneList)
      
      # ---- Enrichment per DB ----
      pathway_result <- NULL
      if (identical(input$pathway_db, "GO")) {
        pathway_result <- tryCatch({
          clusterProfiler::enrichGO(
            gene = selected_genes, OrgDb = orgdb, keyType = "ENTREZID",
            ont = "BP", pAdjustMethod = "BH",
            pvalueCutoff = input$padj_threshold, qvalueCutoff = input$pathway.qval,
            readable = TRUE
          )
        }, error = function(e) {
          shiny::showNotification(paste("GO enrichment failed:", e$message), type = "error"); NULL
        })
        
      } else if (identical(input$pathway_db, "KEGG")) {
        kegg_code <- .get_kegg_code(species)
        if (is.null(kegg_code)) {
          shiny::showNotification("KEGG not supported for this species.", type = "error"); return()
        }
        kegg_raw <- tryCatch({
          clusterProfiler::enrichKEGG(
            gene = selected_genes, organism = kegg_code,
            pvalueCutoff = input$padj_threshold, qvalueCutoff = input$pathway.qval
          )
        }, error = function(e) {
          shiny::showNotification(paste("KEGG enrichment failed:", e$message), type = "error"); NULL
        })
        if (!is.null(kegg_raw) && nrow(as.data.frame(kegg_raw)) > 0) {
          pathway_result <- tryCatch({
            clusterProfiler::setReadable(kegg_raw, OrgDb = orgdb, keyType = "ENTREZID")
          }, error = function(e) {
            shiny::showNotification(paste("setReadable failed:", e$message), type = "error"); kegg_raw
          })
        }
        
      } else if (identical(input$pathway_db, "Reactome")) {
        reactome_code <- .get_reactome_code(species)
        if (is.null(reactome_code)) {
          shiny::showNotification("Reactome enrichment supported only for human/mouse/rat.", type = "error"); return()
        }
        pathway_result <- tryCatch({
          ReactomePA::enrichPathway(
            gene = selected_genes, organism = reactome_code,
            pvalueCutoff = input$padj_threshold, qvalueCutoff = input$pathway.qval,
            readable = TRUE
          )
        }, error = function(e) {
          shiny::showNotification(paste("Reactome enrichment failed:", e$message), type = "error"); NULL
        })
        
      } else {
        shiny::showNotification("Unknown pathway database selected.", type = "error"); return()
      }
      
      # Term-similarity (for emapplot)
      if (!is.null(pathway_result) && nrow(as.data.frame(pathway_result)) > 0) {
        # make gene symbols readable if not yet
        if (!"geneSymbol" %in% colnames(as.data.frame(pathway_result))) {
          pathway_result <- tryCatch({
            clusterProfiler::setReadable(pathway_result, OrgDb = orgdb, keyType = "ENTREZID")
          }, error = function(e) pathway_result)
        }
        # pairwise termsim
        pathway_result <- tryCatch({
          enrichplot::pairwise_termsim(pathway_result)
        }, error = function(e) pathway_result)
      } else {
        shiny::showNotification("No enriched pathways found.", type = "warning")
      }
      
      pathway_result_rv(pathway_result)
    })
    
    # ---- Plots & Downloads ---------------------------------------------------
    tag <- .contrast_tag_from(de_sel)
    output$dotPlot <- shiny::renderPlot({
      shiny::req(pathway_result_rv())
      pr <- pathway_result_rv()
      if (is.null(pr) || nrow(as.data.frame(pr)) == 0) {
        shiny::showNotification("No enrichment terms available for dotplot.", type = "warning"); return(NULL)
      }
      enrichplot::dotplot(pr) + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6, face = "bold"))
    })
    
    output$download_dot_plot <- shiny::downloadHandler(
      filename = function() paste0("Pathway_", input$pathway_db, "_", input$pathway_direction,  "_", tag, "_dot_plot.pdf"),
      content  = function(file) {
        shiny::req(pathway_result_rv())
        pr <- pathway_result_rv()
        grDevices::pdf(file)
        print(enrichplot::dotplot(pr) + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6, face = "bold")))
        grDevices::dev.off()
      }
    )
    
    output$emapPlot <- shiny::renderPlot({
      shiny::req(pathway_result_rv())
      pr <- pathway_result_rv()
      df <- as.data.frame(pr)
      if (is.null(df) || nrow(df) < 2 || !"ID" %in% colnames(df) || !"geneID" %in% colnames(df)) {
        shiny::showNotification("Insufficient enrichment terms or missing columns for emapplot.", type = "warning"); return(NULL)
      }
      if (is.null(pr@termsim) || all(is.na(pr@termsim)) || all(pr@termsim == 0)) {
        shiny::showNotification("No valid term similarity matrix for emapplot.", type = "warning"); return(NULL)
      }
      enrichplot::emapplot(pr, showCategory = 10)
    })
    
    output$download_emap_plot <- shiny::downloadHandler(
      filename = function() paste0("Pathway_", input$pathway_db, "_", input$pathway_direction, "_", tag,  "_emap_plot.pdf"),
      content  = function(file) {
        shiny::req(pathway_result_rv())
        pr <- pathway_result_rv()
        df <- as.data.frame(pr)
        if (is.null(df) || nrow(df) < 2 || is.null(pr@termsim) || all(pr@termsim == 0)) {
          shiny::showNotification("No valid term similarity for export.", type = "warning"); return()
        }
        grDevices::pdf(file); print(enrichplot::emapplot(pr, showCategory = 10)); grDevices::dev.off()
      }
    )
    
    output$cnetPlot <- shiny::renderPlot({
      shiny::req(pathway_result_rv())
      pr <- pathway_result_rv()
      if (is.null(pr) || nrow(as.data.frame(pr)) == 0) {
        shiny::showNotification("No enrichment terms available for cnetplot.", type = "warning"); return(NULL)
      }
      enrichplot::cnetplot(pr, showCategory = 10)
    })
    
    output$download_cnet_plot <- shiny::downloadHandler(
      filename = function() paste0("Pathway_", input$pathway_db, "_", input$pathway_direction, "_", tag,  "_cnet_plot.pdf"),
      content  = function(file) {
        shiny::req(pathway_result_rv())
        pr <- pathway_result_rv()
        grDevices::pdf(file); print(enrichplot::cnetplot(pr, showCategory = 10)); grDevices::dev.off()
      }
    )
    
    output$circularPlot <- shiny::renderPlot({
      shiny::req(pathway_result_rv(), geneList_rv())
      pr <- pathway_result_rv()
      df <- as.data.frame(pr)
      if (is.null(df) || nrow(df) == 0) {
        shiny::showNotification("No enriched terms to show circular plot.", type = "warning"); return(NULL)
      }
      enrichplot::cnetplot(pr,
                           layout = input$circular_layout, foldChange = geneList_rv(),
                           showCategory = 5, circular = TRUE, colorEdge = TRUE
      )
    })
    
    output$download_circular_plot <- shiny::downloadHandler(
      filename = function() paste0("Pathway_", input$pathway_db, "_", input$pathway_direction, "_", input$circular_layout, "_", tag,  "_circular_plot.pdf"),
      content  = function(file) {
        shiny::req(pathway_result_rv(), geneList_rv())
        pr <- pathway_result_rv()
        grDevices::pdf(file)
        print(enrichplot::cnetplot(pr,
                                   layout = input$circular_layout, foldChange = geneList_rv(),
                                   showCategory = 5, circular = TRUE, colorEdge = TRUE
        ))
        grDevices::dev.off()
      }
    )
    
    output$pathwayTable <- DT::renderDT({
      shiny::req(pathway_result_rv())
      DT::datatable(as.data.frame(pathway_result_rv()), options = list(scrollX = TRUE, pageLength = 25))
    })
    
    output$download_pathway_table <- shiny::downloadHandler(
      filename = function() paste0("Pathway_", input$pathway_db, "_", input$pathway_direction, "_", tag,  "_results.csv"),
      content  = function(file) {
        shiny::req(pathway_result_rv())
        utils::write.csv(as.data.frame(pathway_result_rv()), file, row.names = FALSE)
      }
    )
  })
}

# ---- helpers (internal) ------------------------------------------------------

.is_symbol <- function(ids) {
  # Use app's is_symbol() if present; else heuristic (not ENSEMBL-like)
  if (exists("is_symbol", mode = "function")) return(get("is_symbol")(ids))
  !grepl("^ENS[A-Z]*\\d+", ids)
}

.get_orgdb <- function(species) {
  # Try app helper first
  if (exists("get_orgdb", mode = "function")) return(get("get_orgdb")(species))
  # Fallback: support Hs/Mm/Rn (load if installed)
  sp <- as.character(species)
  if (sp %in% c("Homo sapiens", "human")) {
    if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) return(org.Hs.eg.db::org.Hs.eg.db)
  } else if (sp %in% c("Mus musculus", "mouse")) {
    if (requireNamespace("org.Mm.eg.db", quietly = TRUE)) return(org.Mm.eg.db::org.Mm.eg.db)
  } else if (sp %in% c("Rattus norvegicus", "rat")) {
    if (requireNamespace("org.Rn.eg.db", quietly = TRUE)) return(org.Rn.eg.db::org.Rn.eg.db)
  }
  NULL
}

.get_kegg_code <- function(species) {
  if (exists("get_kegg_code", mode = "function")) return(get("get_kegg_code")(species))
  sp <- as.character(species)
  if (sp %in% c("Homo sapiens", "human")) return("hsa")
  if (sp %in% c("Mus musculus", "mouse")) return("mmu")
  if (sp %in% c("Rattus norvegicus", "rat")) return("rno")
  if (sp %in% c("Danio rerio")) return("dre")
  if (sp %in% c("Bos taurus")) return("bta")
  NULL
}

.get_reactome_code <- function(species) {
  if (exists("get_reactome_code", mode = "function")) return(get("get_reactome_code")(species))
  sp <- as.character(species)
  if (sp %in% c("Homo sapiens", "human")) return("human")
  if (sp %in% c("Mus musculus", "mouse")) return("mouse")
  if (sp %in% c("Rattus norvegicus", "rat")) return("rat")
  NULL
}
