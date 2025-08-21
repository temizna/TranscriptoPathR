# R/mod_pathway_plots_server.R
#' Pathway Plots Server (heatmap / tree / upset)
#'
#' Consumes an enrichment result (from Pathway Analysis) and ranked gene vector,
#' draws enrichplot::heatplot / treeplot / upsetplot, and provides PDF downloads.
#'
#' @param id module id (must match mod_pathway_plots_ui)
#' @param pathway_result_rv reactiveVal holding a clusterProfiler/ReactomePA result
#' @param geneList_rv reactiveVal with a named numeric vector of log2FC (ENTREZID names)
#' @param de_sel list returned by mod_de_server() with ref_level()/test_level() (optional; for filenames)
#'
#' @importFrom shiny moduleServer req renderPlot downloadHandler showNotification
#' @importFrom grDevices pdf dev.off
#' @importFrom utils write.csv
#' @importFrom enrichplot heatplot treeplot upsetplot pairwise_termsim
#' @export
mod_pathway_plots_server <- function(id, pathway_result_rv, geneList_rv, de_sel = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    # ---- Heatmap (screen)
    output$pathheatmapPlot <- shiny::renderPlot({
      shiny::req(pathway_result_rv(), geneList_rv())
      pr <- pathway_result_rv()
      gl <- geneList_rv()
      enrichplot::heatplot(pr, foldChange = gl, showCategory = 5)
    })
    
    # ---- Tree plot (screen)  (needs termsim; compute if missing)
    output$treePlot <- shiny::renderPlot({
      shiny::req(pathway_result_rv())
      pr <- pathway_result_rv()
      
      df <- as.data.frame(pr)
      if (is.null(df) || nrow(df) < 2) {
        shiny::showNotification("No enrichment terms available for treeplot.", type = "warning")
        return(NULL)
      }
      
      # ensure term similarity is present
      if (is.null(pr@termsim) || all(is.na(pr@termsim)) || all(pr@termsim == 0)) {
        pr_try <- try(enrichplot::pairwise_termsim(pr), silent = TRUE)
        if (inherits(pr_try, "try-error")) {
          shiny::showNotification("No valid term similarity for treeplot.", type = "warning"); return(NULL)
        } else {
          pr <- pr_try
        }
      }
      
      # choose safe parameters
      n_terms  <- nrow(df)
      show_cat <- min(30, max(5, n_terms))    # plot up to 30
      nCluster <- min(4, max(1, ceiling(show_cat/10)))  # 1..4 clusters
      nWords   <- 4                            # clamp to [1,4]
      
      p <- try(
        enrichplot::treeplot(
          pr,
          showCategory = show_cat,
          nCluster     = nCluster,
          nWords       = nWords
        ),
        silent = TRUE
      )
      
      if (inherits(p, "try-error")) {
        # graceful fallback
        shiny::showNotification("treeplot failed on this result; showing emapplot instead.", type = "message")
        p <- enrichplot::emapplot(pr, showCategory = min(20, n_terms))
      }
      print(p)
    })
    
    # ---- Upset plot (screen)
    output$upsetPlot <- shiny::renderPlot({
      shiny::req(pathway_result_rv())
      pr <- pathway_result_rv()
      df <- as.data.frame(pr)
      if (is.null(df) || nrow(df) == 0) {
        shiny::showNotification("No enrichment terms available for upsetplot.", type = "warning"); return(NULL)
      }
      enrichplot::upsetplot(pr)
    })
    
    # ---- Downloads (PDF)
    output$download_pathheatmap_plot <- shiny::downloadHandler(
      filename = function() {
        tag <- .contrast_tag_from(de_sel)
        paste0("Pathway_", tag, "_heatmap_plot.pdf")
      },
      content = function(file) {
        shiny::req(pathway_result_rv(), geneList_rv())
        pr <- pathway_result_rv(); gl <- geneList_rv()
        grDevices::pdf(file, width = 8.5, height = 11)
        print(enrichplot::heatplot(pr, foldChange = gl, showCategory = 5))
        grDevices::dev.off()
      }
    )
    
    output$download_tree_plot <- shiny::downloadHandler(
      filename = function() {
        tag <- .contrast_tag_from(de_sel)
        paste0("Pathway_", tag, "_tree_plot.pdf")
      },
      content = function(file) {
        shiny::req(pathway_result_rv())
        pr <- pathway_result_rv()
        df <- as.data.frame(pr)
        grDevices::pdf(file, width = 8.5, height = 11)
        if (is.null(df) || nrow(df) < 2) {
          grid::grid.newpage(); grid::grid.text("No enrichment terms for treeplot."); grDevices::dev.off(); return()
        }
        if (is.null(pr@termsim) || all(is.na(pr@termsim)) || all(pr@termsim == 0)) {
          pr_try <- try(enrichplot::pairwise_termsim(pr), silent = TRUE)
          if (!inherits(pr_try, "try-error")) pr <- pr_try else {
            grid::grid.newpage(); grid::grid.text("No valid term similarity for treeplot."); grDevices::dev.off(); return()
          }
        }
        print(enrichplot::treeplot(pr))
        grDevices::dev.off()
      }
    )
    
    output$download_upset_plot <- shiny::downloadHandler(
      filename = function() {
        tag <- .contrast_tag_from(de_sel)
        paste0("Pathway_", tag, "_upset_plot.pdf")
      },
      content = function(file) {
        shiny::req(pathway_result_rv())
        pr <- pathway_result_rv()
        df <- as.data.frame(pr)
        grDevices::pdf(file, width = 8.5, height = 11)
        if (is.null(df) || nrow(df) == 0) {
          grid::grid.newpage(); grid::grid.text("No enrichment terms for upsetplot."); grDevices::dev.off(); return()
        }
        print(enrichplot::upsetplot(pr))
        grDevices::dev.off()
      }
    )
  })
}
