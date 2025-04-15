# === Module: mod_differential_expression ===

#' Differential Expression Module
#'
#' @description Handles DESeq2-based differential expression, renders results table, and provides a clustered heatmap of top genes
#'
#' @param input Shiny input object
#' @param output Shiny output object
#' @param session Shiny session object
#' @param data ReactiveValues with counts, samples, norm_counts, species
#' @param dds_rv ReactiveVal holding DESeqDataSet
#' @param res_reactive ReactiveVal to store DESeq2 result
#' @return None (outputs rendered)
#' @import DESeq2
#' @importFrom dplyr mutate
#' @importFrom DT renderDT datatable
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation
#' @importFrom circlize colorRamp2
#' @importFrom RColorBrewer brewer.pal
#' @importFrom utils head read.csv write.csv str
#' @importFrom stats as.formula dist model.matrix prcomp quantile relevel var
#' @importFrom grDevices dev.off pdf colorRampPalette
#' @importFrom grid gpar
#' @importFrom shiny isolate req renderPlot renderImage showNotification downloadHandler renderPrint
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom shinythemes shinytheme
#' @export
mod_differential_expression <- function(input, output, session, data, dds_rv, res_reactive) {

  observe({
    req(data$samples)
    cols <- colnames(data$samples)
    updateSelectInput(session, "metadata_column", choices = cols)
  })

  observeEvent(input$metadata_column, {
    req(data$samples)
    col <- input$metadata_column
    levels <- unique(as.character(data$samples[[col]]))
    updateSelectInput(session, "reference_condition", choices = levels)
    updateSelectInput(session, "test_condition", choices = levels)
  })

  observeEvent(input$run_de, {
    req(input$metadata_column, input$reference_condition, input$test_condition)
    req(dds_rv())

    dds <- dds_rv()
    dds[[input$metadata_column]] <- relevel(factor(dds[[input$metadata_column]]), ref = input$reference_condition)
    design(dds) <- as.formula(paste("~", input$metadata_column))

    dds <- DESeq(dds)
    res <- results(dds, contrast = c(input$metadata_column, input$test_condition, input$reference_condition))

    res$symbol <- rownames(res)

    if (is_ensembl_id(rownames(res))) {
      conv <- convert_ensembl_to_symbol(rownames(res), data$species)
      res$symbol <- conv[rownames(res)]
    }

    res_df <- as.data.frame(res)
    res_df$symbol[is.na(res_df$symbol)] <- rownames(res_df)[is.na(res_df$symbol)]
    res_df <- res_df[order(res_df$padj), ]

    res_df <- res_df[!is.na(res_df$padj) & !is.na(res_df$log2FoldChange), ]
    res_df <- res_df[abs(res_df$log2FoldChange) >= input$lfc_threshold & res_df$padj <= input$padj_threshold, ]

    res_reactive(res_df)
    
    # Update the dds_rv with the final dds object (if you need to use it elsewhere)
    dds_rv(dds)  # This will update the dds_rv with the new DESeq2 results
    #print(resultsNames(dds)) 
    
  })

  output$deTable <- renderDT({
    req(res_reactive())
    datatable(res_reactive(), options = list(scrollX = TRUE))
  })

  output$download_de_table <- downloadHandler(
    filename = function() "differential_expression_results.csv",
    content = function(file) {
      write.csv(res_reactive(), file, row.names = FALSE)
    }
  )

  output$heatmapPlot <- renderPlot({
    req(res_reactive(), data$norm_counts, data$samples, input$metadata_column)
    top_n <- input$num_genes
    top_genes <- head(res_reactive()[order(res_reactive()$padj), ], top_n)
    expr <- data$norm_counts[rownames(top_genes), ]
    expr <- log2(expr + 1)

    group_values <- data$samples[[input$metadata_column]]
    group_levels <- unique(group_values)
    palette <- colorRampPalette(brewer.pal(8, "Set2"))(length(group_levels))
    ha <- ComplexHeatmap::HeatmapAnnotation(
      df = data.frame(Group = group_values),
      col = list(Group = setNames(palette, group_levels))
    )

    cluster_cols <- if (!is.null(input$cluster_columns)) input$cluster_columns else TRUE

    ComplexHeatmap::Heatmap(expr,
                            name = "log2(norm counts)",
                            top_annotation = ha,
                            cluster_rows = TRUE,
                            cluster_columns = cluster_cols,
                            show_column_names = FALSE,
                            show_row_names = FALSE,
                            column_title = paste("Top", top_n, "Differentially Expressed Genes"),
                            column_title_gp = grid::gpar(fontface = "bold"))
  })

  output$download_heatmap <- downloadHandler(
    filename = function() { "de_heatmap.pdf" },
    content = function(file) {
      pdf(file, width = 10, height = 8)
      top_n <- input$num_genes
      top_genes <- head(res_reactive()[order(res_reactive()$padj), ], top_n)
      expr <- data$norm_counts[rownames(top_genes), ]
      expr <- log2(expr + 1)

      group_values <- data$samples[[input$metadata_column]]
      group_levels <- unique(group_values)
      palette <- colorRampPalette(brewer.pal(8, "Set2"))(length(group_levels))
      ha <- ComplexHeatmap::HeatmapAnnotation(
        df = data.frame(Group = group_values),
        col = list(Group = setNames(palette, group_levels))
      )

      cluster_cols <- if (!is.null(input$cluster_columns)) input$cluster_columns else TRUE

      print(ComplexHeatmap::Heatmap(expr,
                                    name = "log2(norm counts)",
                                    top_annotation = ha,
                                    cluster_rows = TRUE,
                                    cluster_columns = cluster_cols,
                                    show_column_names = FALSE,
                                    show_row_names = FALSE,
                                    column_title = paste("Top", top_n, "Differentially Expressed Genes"),
                                    column_title_gp = grid::gpar(fontface = "bold")))
      dev.off()
    }
  )
}

# === Register in server ===
# mod_differential_expression(input, output, session, data, dds_rv, res_reactive)

