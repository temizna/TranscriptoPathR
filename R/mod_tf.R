#' TF Enrichment Server
#'
#' @param id Shiny module id
#' @param res_reactive reactiveVal with DE results
#' @param filtered_data_rv reactiveValues with at least $species
#' @param cmp Optional comparison bridge from make_cmp_bridge(); used for titles/filenames
#' @export
mod_tf_enrichment_server <- function(id, res_reactive, filtered_data_rv, cmp = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # contrast-aware helpers (fall back if cmp not supplied)
    contrast_label <- shiny::reactive({
      if (!is.null(cmp)) cmp$label() else "Test vs Ref"
    })
    contrast_tag <- shiny::reactive({
      if (!is.null(cmp)) cmp$tag() else "contrast"
    })
    
    processed_tf_data <- shiny::reactiveVal()
    
    load_tf_data <- function(tf_data_source) {
      tf_data_file_map <- list(
        "TRANSFAC and JASPAR PWMs" = "TRANSFAC_and_JASPAR_PWMs.h.gmt.gz",
        "ENCODE and ChEA Consensus TFs from ChIP-X" = "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X.gmt.gz",
        "TRRUST_Transcription_Factors_2019" = "TRRUST_Transcription_Factors_2019h.gmt.gz",
        "TF_Perturbations_Followed_by_Expression" = "TF_Perturbations_Followed_by_Expression_Human.gmt.gz",
        "hTFtarget" = "TF-Target-information.txt.gz",
        "TFLink" = "TFLink_Homo_sapiens_interactions_All_GMT_proteinName_v1.0.gmt.gz"
      )
      tf_data_file <- tf_data_file_map[[tf_data_source]]
      if (is.null(tf_data_file)) stop("Invalid TF data source selected.")
      tf_file_path <- system.file("extdata", tf_data_file, package = "TranscriptoPathR")
      if (tf_file_path == "") stop("TF data file not found.")
      
      if (grepl("\\.gmt", tf_data_file, ignore.case = TRUE)) {
        gson::read.gmt(tf_file_path)
      } else {
        vroom::vroom(tf_file_path, delim = "\t", show_col_types = FALSE)
      }
    }
    
    shiny::observeEvent(input$run_tf_enrichment, {
      shiny::req(res_reactive(), filtered_data_rv$species)
      species <- filtered_data_rv$species
      orgdb   <- get_orgdb(species)
      
      tf_data_source <- input$tf_data_source
      tf_data <- load_tf_data(tf_data_source)
      
      tf_data <- dplyr::rename(tf_data, TF_name = 1, target_gene = 2)
      tf_data_entrez <- suppressMessages(
        clusterProfiler::bitr(tf_data$target_gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = orgdb)
      )
      tf_data_gmt <- tf_data |>
        dplyr::left_join(tf_data_entrez, by = c("target_gene" = "SYMBOL")) |>
        dplyr::select(TF_name, ENTREZID) |>
        dplyr::filter(!is.na(ENTREZID))
      
      res <- res_reactive()
      res <- res[!is.na(res$log2FoldChange) & !is.na(res$padj), , drop = FALSE]
      d1 <- res[, c("log2FoldChange", "padj")]
      d1$gene <- toupper(rownames(res))
      
      if (is_symbol(d1$gene)) {
        d1_ids <- suppressMessages(clusterProfiler::bitr(d1$gene, "SYMBOL", "ENTREZID", OrgDb = orgdb))
      } else {
        d1_ids <- suppressMessages(clusterProfiler::bitr(d1$gene, "ENSEMBL", "ENTREZID", OrgDb = orgdb))
      }
      
      d1_merged <- merge(d1, d1_ids, by.x = "gene", by.y = 1)
      d1_merged <- d1_merged[!duplicated(d1_merged$ENTREZID), ]
      
      direction <- input$gene_direction
      gene_vector <- switch(
        direction,
        "Up"   = d1_merged[d1_merged$log2FoldChange >= input$lfc_threshold & d1_merged$padj <= input$padj_threshold, ],
        "Down" = d1_merged[d1_merged$log2FoldChange <= -input$lfc_threshold & d1_merged$padj <= input$padj_threshold, ],
        d1_merged[abs(d1_merged$log2FoldChange) >= input$lfc_threshold & d1_merged$padj <= input$padj_threshold, ]
      )
      gene_vector <- gene_vector[!is.na(gene_vector$log2FoldChange) & !is.na(gene_vector$ENTREZID), ]
      geneList <- gene_vector$log2FoldChange
      names(geneList) <- gene_vector$ENTREZID
      geneList <- geneList[!duplicated(names(geneList))]
      geneList <- sort(geneList, decreasing = TRUE)
      if (length(geneList) > 1000) geneList <- head(geneList, 1000)
      
      selected_genes <- names(geneList)
      tf_result <-
        if (identical(input$enrichment_method, "Over_representation")) {
          clusterProfiler::enricher(
            gene = selected_genes,
            TERM2GENE = tf_data_gmt,
            pvalueCutoff = 0.05,
            qvalueCutoff = input$tf.qval
          )
        } else {
          tryCatch({
            clusterProfiler::GSEA(
              geneList = geneList,
              TERM2GENE = tf_data_gmt,
              pvalueCutoff = input$tf.qval,
              verbose = FALSE
            )
          }, error = function(e) {
            shiny::showNotification(paste("TF:GSEA failed:", e$message), type = "error")
            NULL
          })
        }
      
      if (is.null(tf_result) || nrow(as.data.frame(tf_result)) < 1) {
        shiny::showNotification("No enriched terms found in TF enrichment.", type = "warning")
        return()
      }
      
      # Titles use contrast_label()
      output$tf_dotplot <- shiny::renderPlot({
        p <- enrichplot::dotplot(tf_result) + ggplot2::labs(title = paste0("TF Enrichment: ", contrast_label()))
        p + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6, face = "bold"))
      })
      
      output$download_tf_dotplot <- shiny::downloadHandler(
        filename = function() {
          paste0("TF_", contrast_tag(), "_", .safe_tag(input$tf_data_source), "_dotplot.pdf")
        },
        content = function(file) {
          p <- enrichplot::dotplot(tf_result) + ggplot2::labs(title = paste0("TF Enrichment: ", contrast_label()))
          ggplot2::ggsave(file, p, width = 10, height = 8, units = "in")
        }
      )
      
      if (identical(input$enrichment_method, "GSEA")) {
        output$tf_ridgeplot <- shiny::renderPlot({
          p <- enrichplot::ridgeplot(tf_result) + ggplot2::labs(title = paste0("TF GSEA: ", contrast_label()))
          p + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6, face = "bold"))
        })
        output$download_tf_ridgeplot <- shiny::downloadHandler(
          filename = function() {
            paste0("TF_", contrast_tag(), "_", .safe_tag(input$tf_data_source), "_ridgeplot.pdf")
          },
          content = function(file) {
            p <- enrichplot::ridgeplot(tf_result) + ggplot2::labs(title = paste0("TF GSEA: ", contrast_label()))
            ggplot2::ggsave(file, p, width = 10, height = 8, units = "in")
          }
        )
      }
      
      output$tf_results_table <- DT::renderDT({
        DT::datatable(as.data.frame(tf_result@result), options = list(scrollX = TRUE))
      })
      output$download_tf_results_table <- shiny::downloadHandler(
        filename = function() paste0("TF_", contrast_tag(), "_", .safe_tag(input$tf_data_source), "_results.csv"),
        content = function(file) utils::write.csv(as.data.frame(tf_result@result), file, row.names = FALSE)
      )
    })
  })
}
