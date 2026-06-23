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
mod_pathway_plots_server <- function(id, pathway_result_rv, geneList_rv, tag_rv = NULL, de_sel = NULL, cmp = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # -------- tag helper (cmp preferred, then de_sel, else "contrast")
    .get_tag <- function() {
      if (!is.null(tag_rv) && is.function(tag_rv)) {
        t <- tag_rv()
        if (!is.null(t) && nzchar(t)) return(t)
      }
      
      if (!is.null(cmp) && is.function(cmp$tag)) {
        t <- cmp$tag()
        if (!is.null(t) && nzchar(t)) return(t)
      }
      
      .contrast_tag_from(de_sel)
    }
    
    .safe_show_category <- function(pr, default = 5L) {
      df <- as.data.frame(pr)
      max(1L, min(as.integer(default), nrow(df)))
    }
    .safe_treeplot <- function(pr) {
      df <- as.data.frame(pr)
      if (is.null(df) || nrow(df) < 2) return(NULL)
      
      if (is.null(pr@termsim) || all(is.na(pr@termsim)) || all(pr@termsim == 0)) {
        pr_try <- try(enrichplot::pairwise_termsim(pr), silent = TRUE)
        if (inherits(pr_try, "try-error")) return(NULL)
        pr <- pr_try
      }
      
      show_cat <- min(30L, nrow(as.data.frame(pr)))
      if (show_cat < 2L) return(NULL)
      
      p <- try(
        suppressWarnings(
          enrichplot::treeplot(
            pr,
            showCategory = show_cat,
            cluster.params = list(
              n = min(4L, max(1L, ceiling(show_cat / 10))),
              label_words_n = 4L
            )
          )
        ),
        silent = TRUE
      )
      
      if (inherits(p, "try-error")) return(NULL)
      p
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
      enrichplot::heatplot(pr, foldChange = gl, showCategory = 10)
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
      
      p <- .safe_treeplot(pr)
      
      if (is.null(p)) {
        shiny::showNotification("treeplot failed; showing emapplot instead.", type = "message")
        pr2 <- try(enrichplot::pairwise_termsim(pr), silent = TRUE)
        if (!inherits(pr2, "try-error")) pr <- pr2
        p <- enrichplot::emapplot(pr, showCategory = min(20L, nrow(as.data.frame(pr))))
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
        on.exit(grDevices::dev.off(), add = TRUE)
        
        if (is.null(df) || nrow(df) < 2) {
          grid::grid.newpage()
          grid::grid.text("No enrichment terms for treeplot.")
          return()
        }
        
        p <- .safe_treeplot(pr)
        
        if (is.null(p)) {
          pr2 <- try(enrichplot::pairwise_termsim(pr), silent = TRUE)
          if (!inherits(pr2, "try-error")) pr <- pr2
          p <- try(enrichplot::emapplot(pr, showCategory = min(20L, nrow(as.data.frame(pr)))), silent = TRUE)
        }
        
        if (inherits(p, "try-error") || is.null(p)) {
          grid::grid.newpage()
          grid::grid.text("treeplot could not be rendered.")
        } else {
          print(p)
        }
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
