#' Transcription Factor Enrichment - Server
#'
#' Uses DE results to run TF over-representation or GSEA against local GMT/TXT
#' resources under inst/extdata. Maps targets to ENTREZ IDs for enrichment.
#'
#' @param id Shiny module id
#' @param res_reactive reactiveVal or reactive with DE results (data.frame
#'        containing at least columns: log2FoldChange, padj; rownames are IDs)
#' @param filtered_data_rv reactiveValues with $species (used for OrgDb selection)
#' @importFrom shiny moduleServer req showNotification renderPlot downloadHandler
#' @importFrom shiny observeEvent
#' @importFrom DT renderDT datatable
#' @importFrom utils read.table write.csv
#' @importFrom clusterProfiler read.gmt enricher GSEA bitr
#' @importFrom enrichplot dotplot ridgeplot
#' @export
mod_tf_enrichment_server <- function(id, res_reactive, filtered_data_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Map friendly names to files in inst/extdata
    tf_file_map <- list(
      "TRANSFAC and JASPAR PWMs" = "TRANSFAC_and_JASPAR_PWMs.h.gmt.gz",
      "ENCODE and ChEA Consensus TFs from ChIP-X" = "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X.gmt.gz",
      "TRRUST_Transcription_Factors_2019" = "TRRUST_Transcription_Factors_2019h.gmt.gz",
      "TF_Perturbations_Followed_by_Expression" = "TF_Perturbations_Followed_by_Expression_Human.gmt.gz",
      "hTFtarget" = "TF-Target-information.txt.gz",
      "TFLink" = "TFLink_Homo_sapiens_interactions_All_GMT_proteinName_v1.0.gmt.gz"
    )
    
    # Helper: load TF term-to-gene mapping (SYMBOL)
    load_tf_data <- function(source_label) {
      fn <- tf_file_map[[source_label]]
      if (is.null(fn)) stop("Invalid TF data source.")
      path <- system.file("extdata", fn, package = "TranscriptoPathR")
      if (identical(path, "")) stop("TF data file not found in inst/extdata.")
      if (grepl("\\.gmt(\\.gz)?$", fn, ignore.case = TRUE)) {
        # clusterProfiler::read.gmt returns columns: term, gene
        g <- clusterProfiler::read.gmt(path)
        colnames(g) <- c("TF_name", "target_gene")
        return(g[, c("TF_name", "target_gene")])
      } else {
        # TXT table; try to find TF and Target columns
        tb <- utils::read.table(path, sep = "\t", header = TRUE, quote = "", comment.char = "",
                                stringsAsFactors = FALSE, check.names = FALSE)
        tf_col <- intersect(c("TF_name", "TF", "TF.Name", "TFname"), colnames(tb))
        tg_col <- intersect(c("target_gene", "Target", "Gene.symbol", "Gene", "symbol"), colnames(tb))
        if (!length(tf_col) || !length(tg_col)) {
          stop("Could not find TF and target columns in TXT file.")
        }
        df <- data.frame(TF_name = tb[[tf_col[1]]],
                         target_gene = tb[[tg_col[1]]],
                         stringsAsFactors = FALSE)
        df <- df[!is.na(df$TF_name) & !is.na(df$target_gene) & nzchar(df$target_gene), ]
        return(df)
      }
    }
    
    # Main analysis
    shiny::observeEvent(input$run_tf_enrichment, {
      shiny::req(res_reactive(), filtered_data_rv$species)
      species <- filtered_data_rv$species
      orgdb <- get_orgdb(species)
      
      # Load TF data (SYMBOL targets)
      tf_df <- try(load_tf_data(input$tf_data_source), silent = TRUE)
      if (inherits(tf_df, "try-error")) {
        shiny::showNotification(paste("TF data load failed:", as.character(tf_df)), type = "error")
        return()
      }
      
      # Map TF target SYMBOL -> ENTREZ
      map_tbl <- suppressMessages(clusterProfiler::bitr(
        unique(tf_df$target_gene),
        fromType = "SYMBOL",
        toType   = "ENTREZID",
        OrgDb    = orgdb
      ))
      if (!nrow(map_tbl)) {
        shiny::showNotification("Could not map TF target genes to ENTREZ IDs.", type = "error")
        return()
      }
      term2gene <- merge(tf_df, map_tbl, by.x = "target_gene", by.y = "SYMBOL", all.x = TRUE)
      term2gene <- term2gene[!is.na(term2gene$ENTREZID), c("TF_name", "ENTREZID")]
      colnames(term2gene) <- c("term", "gene")
      
      # Build geneList from DE results
      res <- res_reactive()
      res <- res[!is.na(res$log2FoldChange) & !is.na(res$padj), , drop = FALSE]
      d1 <- data.frame(log2FoldChange = res$log2FoldChange,
                       padj = res$padj,
                       gene = rownames(res),
                       stringsAsFactors = FALSE)
      
      # If rownames are Ensembl, convert to SYMBOL; else keep SYMBOL
      if (!is_symbol(d1$gene)) {
        conv <- convert_ensembl_to_symbol(d1$gene, species)
        keep <- !is.na(conv)
        d1 <- d1[keep, , drop = FALSE]
        d1$gene <- conv[keep]
      } else {
        # normalize to uppercase for robustness with human TF sets
        d1$gene <- toupper(d1$gene)
      }
      
      # Map gene SYMBOL -> ENTREZ for the DE list
      d1_ids <- suppressMessages(clusterProfiler::bitr(
        d1$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = orgdb
      ))
      d1m <- merge(d1, d1_ids, by.x = "gene", by.y = "SYMBOL")
      d1m <- d1m[!duplicated(d1m$ENTREZID), ]
      
      # Directional filter
      dir_sel <- switch(
        input$gene_direction,
        "Up"   = d1m[d1m$log2FoldChange >=  input$lfc_threshold & d1m$padj <= input$padj_threshold, ],
        "Down" = d1m[d1m$log2FoldChange <= -input$lfc_threshold & d1m$padj <= input$padj_threshold, ],
        d1m[abs(d1m$log2FoldChange) >= input$lfc_threshold & d1m$padj <= input$padj_threshold, ]
      )
      dir_sel <- dir_sel[!is.na(dir_sel$ENTREZID), ]
      if (!nrow(dir_sel)) {
        shiny::showNotification("No genes pass the selected thresholds.", type = "warning")
        return()
      }
      
      # Gene list (named by ENTREZ)
      geneList <- dir_sel$log2FoldChange
      names(geneList) <- dir_sel$ENTREZID
      geneList <- geneList[!duplicated(names(geneList))]
      geneList <- sort(geneList, decreasing = TRUE)
      if (!length(geneList)) {
        shiny::showNotification("Gene list empty after processing.", type = "warning")
        return()
      }
      
      # Enrichment
      if (identical(input$enrichment_method, "Over_representation")) {
        tf_res <- clusterProfiler::enricher(
          gene = names(geneList),
          TERM2GENE = term2gene,
          pvalueCutoff = 0.05,
          qvalueCutoff = input$tf.qval
        )
      } else {
        tf_res <- try(clusterProfiler::GSEA(
          geneList = geneList,
          TERM2GENE = term2gene,
          pvalueCutoff = input$tf.qval,
          verbose = FALSE
        ), silent = TRUE)
        if (inherits(tf_res, "try-error")) {
          shiny::showNotification(paste("TF GSEA failed:", as.character(tf_res)), type = "error")
          return()
        }
      }
      
      if (is.null(tf_res) || !nrow(as.data.frame(tf_res))) {
        shiny::showNotification("No enriched TF terms found.", type = "warning")
        return()
      }
      
      # Plots
      output$tf_dotplot <- shiny::renderPlot({
        enrichplot::dotplot(tf_res) + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6, face = "bold"))
      })
      output$download_tf_dotplot <- shiny::downloadHandler(
        filename = function() paste0("TF_Enrichment_", gsub("\\W+", "_", input$tf_data_source), "_dotplot.pdf"),
        content  = function(file) {
          p <- enrichplot::dotplot(tf_res) + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6, face = "bold"))
          ggplot2::ggsave(file, p, width = 10, height = 8, units = "in")
        }
      )
      
      output$tf_ridgeplot <- shiny::renderPlot({
        if (identical(input$enrichment_method, "GSEA")) {
          enrichplot::ridgeplot(tf_res) + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6, face = "bold"))
        }
      })
      output$download_tf_ridgeplot <- shiny::downloadHandler(
        filename = function() paste0("TF_Enrichment_", gsub("\\W+", "_", input$tf_data_source), "_ridgeplot.pdf"),
        content  = function(file) {
          if (identical(input$enrichment_method, "GSEA")) {
            p <- enrichplot::ridgeplot(tf_res) + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6, face = "bold"))
            ggplot2::ggsave(file, p, width = 10, height = 8, units = "in")
          }
        }
      )
      
      # Table
      output$tf_results_table <- DT::renderDT({
        DT::datatable(as.data.frame(tf_res@result), options = list(scrollX = TRUE, pageLength = 25))
      })
      output$download_tf_results_table <- shiny::downloadHandler(
        filename = function() paste0("TF_Enrichment_", gsub("\\W+", "_", input$tf_data_source), "_results.csv"),
        content  = function(file) utils::write.csv(as.data.frame(tf_res@result), file, row.names = FALSE)
      )
    })
  })
}