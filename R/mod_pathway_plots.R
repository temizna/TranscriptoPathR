# R/mod_pathway_plots_server.R

#' Pathway Plots Server (heatmap / tree / upset)
#'
#' Consumes an enrichment result (from Pathway Analysis) and ranked gene vector,
#' draws enrichplot::heatplot / treeplot / upsetplot, and provides PDF downloads.
#'
#' @param id module id (must match mod_pathway_plots_ui)
#' @param pathway_result_rv reactiveVal holding a clusterProfiler/ReactomePA result
#' @param geneList_rv reactiveVal with a named numeric vector of log2FC (ENTREZID names)
#' @param de_sel optional list returned by mod_de_server() with ref_level()/test_level() for filenames
#' @param cmp optional comparison bridge from make_cmp_bridge(); if provided, used for filenames
#'
#' @importFrom shiny moduleServer req renderPlot downloadHandler showNotification
#' @importFrom grDevices pdf dev.off
#' @importFrom utils write.csv
#' @importFrom enrichplot heatplot treeplot upsetplot pairwise_termsim emapplot
#' @importFrom grid grid.newpage grid.text
#' @export
mod_pathway_plots_server <- function(id, pathway_result_rv, geneList_rv, de_sel = NULL, cmp = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # -------- tag helper (cmp preferred, then de_sel, else "contrast")
    .safe_tag <- function(x) gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
    .get_tag <- function() {
      # cmp$tag() if available
      t <- try(if (!is.null(cmp) && is.function(cmp$tag)) cmp$tag(), silent = TRUE)
      if (!inherits(t, "try-error") && !is.null(t) && nzchar(t)) return(t)
      
      # fall back to de_sel ref/test if present
      if (!is.null(de_sel) && is.function(de_sel$ref_level) && is.function(de_sel$test_level)) {
        ref  <- try(de_sel$ref_level(),  silent = TRUE)
        test <- try(de_sel$test_level(), silent = TRUE)
        if (!inherits(ref, "try-error") && !inherits(test, "try-error") &&
            !is.null(ref) && !is.null(test) && nzchar(ref) && nzchar(test)) {
          return(paste0(.safe_tag(test), "_vs_", .safe_tag(ref)))
        }
      }
      "contrast"
    }
    
    # -------- Heatmap
    output$pathheatmapPlot <- shiny::renderPlot({
      shiny::req(pathway_result_rv(), geneList_rv())
      pr <- pathway_result_rv()
      gl <- geneList_rv()
      df <- as.data.frame(pr)
      if (is.null(df) || nrow(df) == 0) {
        shiny::showNotification("No enrichment terms available for heatmap.", type = "warning")
        return(NULL)
      }
      enrichplot::heatplot(pr, foldChange = gl, showCategory = 5)
    })
    
    # -------- Tree plot (guard against missing termsim and bad k/nWords)
    output$treePlot <- shiny::renderPlot({
      shiny::req(pathway_result_rv())
      pr <- pathway_result_rv()
      df <- as.data.frame(pr)
      if (is.null(df) || nrow(df) < 2) {
        shiny::showNotification("No enrichment terms available for treeplot.", type = "warning")
        return(NULL)
      }
      
      # ensure term similarity
      if (is.null(pr@termsim) || all(is.na(pr@termsim)) || all(pr@termsim == 0)) {
        pr_try <- try(enrichplot::pairwise_termsim(pr), silent = TRUE)
        if (inherits(pr_try, "try-error")) {
          shiny::showNotification("No valid term similarity for treeplot.", type = "warning")
          return(NULL)
        } else {
          pr <- pr_try
        }
      }
      
      # safe parameters
      n_terms  <- nrow(df)
      show_cat <- min(30, max(5, n_terms))
      nCluster <- min(4, max(1, ceiling(show_cat / 10)))
      nWords   <- 4
      
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
        shiny::showNotification("treeplot failed; showing emapplot as fallback.", type = "message")
        p <- enrichplot::emapplot(pr, showCategory = min(20, n_terms))
      }
      print(p)
    })
    
    # -------- Upset plot
    output$upsetPlot <- shiny::renderPlot({
      shiny::req(pathway_result_rv())
      pr <- pathway_result_rv()
      df <- as.data.frame(pr)
      if (is.null(df) || nrow(df) == 0) {
        shiny::showNotification("No enrichment terms available for upsetplot.", type = "warning")
        return(NULL)
      }
      enrichplot::upsetplot(pr)
    })
    
    # -------- Downloads
    output$download_pathheatmap_plot <- shiny::downloadHandler(
      filename = function() paste0("Pathway_", .get_tag(), "_heatmap_plot.pdf"),
      content = function(file) {
        shiny::req(pathway_result_rv(), geneList_rv())
        pr <- pathway_result_rv(); gl <- geneList_rv()
        grDevices::pdf(file, width = 8.5, height = 11)
        print(enrichplot::heatplot(pr, foldChange = gl, showCategory = 5))
        grDevices::dev.off()
      }
    )
    
    output$download_tree_plot <- shiny::downloadHandler(
      filename = function() paste0("Pathway_", .get_tag(), "_tree_plot.pdf"),
      content = function(file) {
        shiny::req(pathway_result_rv())
        pr <- pathway_result_rv()
        df <- as.data.frame(pr)
        
        grDevices::pdf(file, width = 8.5, height = 11)
        
        if (is.null(df) || nrow(df) < 2) {
          grid::grid.newpage(); grid::grid.text("No enrichment terms for treeplot."); grDevices::dev.off(); return()
        }
        
        # ensure term similarity for PDF, same as screen
        if (is.null(pr@termsim) || all(is.na(pr@termsim)) || all(pr@termsim == 0)) {
          pr_try <- try(enrichplot::pairwise_termsim(pr), silent = TRUE)
          if (!inherits(pr_try, "try-error")) pr <- pr_try else {
            grid::grid.newpage(); grid::grid.text("No valid term similarity for treeplot."); grDevices::dev.off(); return()
          }
        }
        
        n_terms  <- nrow(df)
        show_cat <- min(30, max(5, n_terms))
        nCluster <- min(4, max(1, ceiling(show_cat / 10)))
        nWords   <- 4
        
        p <- try(
          enrichplot::treeplot(pr, showCategory = show_cat, nCluster = nCluster, nWords = nWords),
          silent = TRUE
        )
        if (inherits(p, "try-error")) {
          grid::grid.newpage()
          grid::grid.text("treeplot failed; try emapplot in the app.")
        } else {
          print(p)
        }
        grDevices::dev.off()
      }
    )
    
    output$download_upset_plot <- shiny::downloadHandler(
      filename = function() paste0("Pathway_", .get_tag(), "_upset_plot.pdf"),
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
