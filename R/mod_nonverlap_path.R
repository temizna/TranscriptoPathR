# --- Module Server: Non-overlap Genes Pathway Analysis -----------------------
#' Non-overlap Pathway Analysis Server
#'
#' Removes genes already present in current enrichment results and re-runs
#' enrichment on the remaining (non-overlap) set.
#'
#' @param id module id
#' @param filtered_data_rv reactiveValues with $species
#' @param pathway_result_rv reactiveVal holding the *initial* enrichment result (enrichResult)
#' @param geneList_rv reactiveVal named numeric vector (log2FC), names = ENTREZID
#' @param geneListU_rv reactiveVal to store the updated non-overlap gene list (named by ENTREZID)
#'
#' @importFrom shiny moduleServer req observeEvent showNotification renderPlot downloadHandler
#' @importFrom DT renderDT datatable
#' @importFrom utils write.csv
#' @importFrom grDevices pdf dev.off
#' @importFrom clusterProfiler enrichGO enrichKEGG setReadable bitr
#' @importFrom ReactomePA enrichPathway
#' @importFrom enrichplot dotplot heatplot treeplot pairwise_termsim
#' @export
mod_nonoverlap_server <- function(
    id,
    filtered_data_rv,
    pathway_result_rv,
    geneList_rv,
    geneListU_rv
) {
  shiny::moduleServer(id, function(input, output, session) {
    
    # -- Run analysis
    shiny::observeEvent(input$run_non_overlap_pathway, {
      shiny::req(geneList_rv(), pathway_result_rv(), filtered_data_rv$species)
      
      # Pull current enrichment results as a data.frame
      pr <- pathway_result_rv()
      pr_df <- as.data.frame(pr)
      if (!nrow(pr_df) || !"geneID" %in% colnames(pr_df)) {
        shiny::showNotification("No prior pathway result with geneID column.", type = "error")
        return()
      }
      
      # Genes included in the current pathways (ENTREZ IDs as 'a/b/c')
      pathway_genes <- unique(unlist(strsplit(pr_df$geneID, "/", fixed = TRUE)))
      
      # Original gene list (named by ENTREZID)
      gl0 <- geneList_rv()
      if (!length(gl0)) {
        shiny::showNotification("Initial gene list is empty.", type = "error"); return()
      }
      
      # Non-overlap: remove genes already in enriched pathways
      gl_non <- gl0[setdiff(names(gl0), pathway_genes)]
      if (!length(gl_non)) {
        shiny::showNotification("No genes left after removing overlap.", type = "warning"); return()
      }
      geneListU_rv(gl_non)
      
      species <- filtered_data_rv$species
      orgdb   <- get_orgdb(species)               # your existing helper
      showN   <- max(3, min(input$showCategory_nonOL, 50))
      padjC   <- input$padj_threshold_nonOL
      qvalC   <- input$pathway_qval_nonOL
      
      selected_genes <- names(gl_non)
      if (length(selected_genes) < 10) {
        shiny::showNotification("Too few genes left for enrichment (<10).", type = "warning")
        return()
      }
      
      # -- Enrichment by chosen DB
      pathway_result <- NULL
      db <- input$pathway_db_nonOL
      
      if (identical(db, "GO")) {
        pathway_result <- tryCatch({
          clusterProfiler::enrichGO(
            gene = selected_genes,
            OrgDb = orgdb, keyType = "ENTREZID",
            ont = "BP", pAdjustMethod = "BH",
            pvalueCutoff = padjC, qvalueCutoff = qvalC,
            readable = TRUE
          )
        }, error = function(e) { shiny::showNotification(e$message, type = "error"); NULL })
      } else if (identical(db, "KEGG")) {
        kegg_code <- get_kegg_code(species)
        if (is.null(kegg_code)) {
          shiny::showNotification("KEGG not supported for this species.", type = "error"); return()
        }
        x <- tryCatch({
          clusterProfiler::enrichKEGG(
            gene = selected_genes, organism = kegg_code,
            pvalueCutoff = padjC, qvalueCutoff = qvalC
          )
        }, error = function(e) { shiny::showNotification(e$message, type = "error"); NULL })
        if (!is.null(x)) {
          pathway_result <- tryCatch({
            clusterProfiler::setReadable(x, OrgDb = orgdb, keyType = "ENTREZID")
          }, error = function(e) { shiny::showNotification(e$message, type = "error"); NULL })
        }
      } else if (identical(db, "Reactome")) {
        rc <- get_reactome_code(species)
        if (is.null(rc)) {
          shiny::showNotification("Reactome not supported for this species.", type = "error"); return()
        }
        pathway_result <- tryCatch({
          ReactomePA::enrichPathway(
            gene = selected_genes, organism = rc,
            pvalueCutoff = padjC, qvalueCutoff = qvalC, readable = TRUE
          )
        }, error = function(e) { shiny::showNotification(e$message, type = "error"); NULL })
      }
      
      if (is.null(pathway_result) || !nrow(as.data.frame(pathway_result))) {
        shiny::showNotification("No enriched pathways found for non-overlap set.", type = "warning")
        # Clear outputs
        output$dotPlot_nonOL        <- shiny::renderPlot({})
        output$heatmapPlot_nonOL    <- shiny::renderPlot({})
        output$treePlot_nonOL       <- shiny::renderPlot({})
        output$pathwayTable_nonOL   <- DT::renderDT({})
        return()
      }
      
      # Try term similarity (guarded)
      pr_termsim <- try(enrichplot::pairwise_termsim(pathway_result), silent = TRUE)
      if (!inherits(pr_termsim, "try-error")) pathway_result <- pr_termsim
      
      # ---- Plots/Tables
      output$dotPlot_nonOL <- shiny::renderPlot({
        enrichplot::dotplot(pathway_result, showCategory = showN)
      })
      
      output$heatmapPlot_nonOL <- shiny::renderPlot({
        enrichplot::heatplot(pathway_result, foldChange = gl_non, showCategory = min(5, showN))
      })
      
      output$treePlot_nonOL <- shiny::renderPlot({
        # treeplot can be picky; protect with tryCatch
        tryCatch({
          enrichplot::treeplot(pathway_result, showCategory = showN)
        }, error = function(e) {
          shiny::showNotification(
            paste("Treeplot unavailable:", e$message),
            type = "warning"
          )
          NULL
        })
      })
      
      output$pathwayTable_nonOL <- DT::renderDT({
        DT::datatable(as.data.frame(pathway_result), options = list(scrollX = TRUE))
      })
      
      # ---- Downloads
      output$download_nonOL_dot_plot <- shiny::downloadHandler(
        filename = function() paste0("Pathway_nonoverlap_", db, "_dotplot.pdf"),
        content  = function(file) {
          grDevices::pdf(file, width = 8, height = 6)
          print(enrichplot::dotplot(pathway_result, showCategory = showN))
          grDevices::dev.off()
        }
      )
      
      output$download_nonOL_heatmap_plot <- shiny::downloadHandler(
        filename = function() paste0("Pathway_nonoverlap_", db, "_heatmap.pdf"),
        content  = function(file) {
          grDevices::pdf(file, width = 8, height = 6)
          print(enrichplot::heatplot(pathway_result, foldChange = gl_non, showCategory = min(5, showN)))
          grDevices::dev.off()
        }
      )
      
      output$download_nonOL_tree_plot <- shiny::downloadHandler(
        filename = function() paste0("Pathway_nonoverlap_", db, "_treeplot.pdf"),
        content  = function(file) {
          grDevices::pdf(file, width = 8, height = 6)
          p <- tryCatch({
            enrichplot::treeplot(pathway_result, showCategory = showN)
          }, error = function(e) NULL)
          if (!is.null(p)) print(p)
          grDevices::dev.off()
        }
      )
      
      output$download_nonOL_pathway_table <- shiny::downloadHandler(
        filename = function() paste0("Pathway_nonoverlap_", db, "_results.csv"),
        content  = function(file) utils::write.csv(as.data.frame(pathway_result), file, row.names = FALSE)
      )
    })
  })
}
