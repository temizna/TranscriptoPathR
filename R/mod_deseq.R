# === Module: mod_de_server ===
#' Differential Expression Server Module (DESeq2)
#'
#' Runs DESeq2, renders a results table, and draws a ComplexHeatmap of top genes.
#' Exposes selected contrast so other modules (e.g., GSVA) can reuse it.
#'
#' @param id Shiny module id (must match mod_de_tab_ui)
#' @param filtered_data_rv reactiveValues with $samples, $norm_counts, $species
#' @param filtered_dds_rv reactiveVal holding a DESeqDataSet (subsetted)
#' @param res_reactive reactiveVal to store DESeq2 result data.frame
#' @return list of reactives: group_var, ref_level, test_level
#'
#' @import DESeq2
#' @importFrom DT renderDT datatable
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation draw
#' @importFrom RColorBrewer brewer.pal
#' @importFrom circlize colorRamp2
#' @importFrom grid gpar
#' @importFrom utils write.csv
#' @importFrom stats as.formula relevel aggregate
#' @importFrom grDevices pdf dev.off
#' @importFrom shiny req showNotification downloadHandler renderPlot observe observeEvent
#' @export
mod_de_server <- function(id, filtered_data_rv, filtered_dds_rv, res_reactive) {
  moduleServer(id, function(input, output, session) {
    
    # Populate "Variable to test"
    observe({
      req(filtered_data_rv$samples)
      cols <- colnames(filtered_data_rv$samples)
      updateSelectInput(session, "metadata_column", choices = cols,
                        selected = if (length(cols)) cols[1] else character(0))
    })
    
    # Populate levels for Reference/Test when variable changes
    observeEvent(input$metadata_column, {
      req(filtered_data_rv$samples, input$metadata_column)
      lv <- unique(as.character(filtered_data_rv$samples[[input$metadata_column]]))
      updateSelectInput(session, "reference_condition", choices = lv,
                        selected = if (length(lv)) lv[1] else character(0))
      updateSelectInput(session, "test_condition", choices = lv,
                        selected = if (length(lv) >= 2) lv[2] else if (length(lv)) lv[1] else character(0))
    }, ignoreInit = TRUE)
    
    # Run DESeq2
    observeEvent(input$run_de, {
      req(input$metadata_column, input$reference_condition, input$test_condition, filtered_dds_rv())
      if (identical(input$reference_condition, input$test_condition)) {
        showNotification("Reference and Test conditions must be different.", type = "error"); return()
      }
      
      dds <- filtered_dds_rv()
      dds[[input$metadata_column]] <- stats::relevel(factor(dds[[input$metadata_column]]),
                                                     ref = input$reference_condition)
      design(dds) <- as.formula(paste("~", input$metadata_column))
      
      dds <- suppressMessages(DESeq(dds))
      res <- suppressMessages(results(dds, contrast = c(input$metadata_column,
                                                        input$test_condition,
                                                        input$reference_condition)))
      
      # Gene symbol column
      res$symbol <- rownames(res)
      if (is_ensembl_id(rownames(res))) {
        conv <- convert_ensembl_to_symbol(rownames(res), filtered_data_rv$species)
        res$symbol <- conv[rownames(res)]
      }
      
      # Filter and order
      res_df <- as.data.frame(res)
      res_df$symbol[is.na(res_df$symbol)] <- rownames(res_df)[is.na(res_df$symbol)]
      res_df <- res_df[!is.na(res_df$padj) & !is.na(res_df$log2FoldChange), ]
      res_df <- res_df[abs(res_df$log2FoldChange) >= input$lfc_threshold &
                         res_df$padj <= input$padj_threshold, ]
      res_df <- res_df[order(res_df$padj), ]
      
      res_reactive(res_df)
      filtered_dds_rv(dds)  # keep updated object if needed downstream
    })
    
    # Results table
    output$deTable <- DT::renderDT({
      req(res_reactive())
      DT::datatable(res_reactive(), options = list(scrollX = TRUE, pageLength = 25))
    })
    
    # Heatmap (screen)
    output$heatmapPlot <- renderPlot({
      req(res_reactive(), filtered_data_rv$norm_counts, filtered_data_rv$samples)
      top_n <- input$num_genes
      meta_col <- input$metadata_column
      cluster_cols <- isTRUE(input$cluster_columns)
      
      top_genes <- head(res_reactive()[order(res_reactive()$padj), ], top_n)
      req(nrow(top_genes) >= 2)
      
      expr <- filtered_data_rv$norm_counts[rownames(top_genes), , drop = FALSE]
      expr <- log2(expr + 1)
      
      group_values <- as.character(filtered_data_rv$samples[[meta_col]])
      group_levels <- unique(group_values)
      
      palette <- if (length(group_levels) == 1) {
        setNames("#1f78b4", group_levels)
      } else if (length(group_levels) == 2) {
        setNames(RColorBrewer::brewer.pal(3, "Dark2")[1:2], group_levels)
      } else if (length(group_levels) <= 8) {
        setNames(RColorBrewer::brewer.pal(length(group_levels), "Dark2"), group_levels)
      } else {
        setNames(colorspace::rainbow_hcl(length(group_levels)), group_levels)
      }
      
      ha <- ComplexHeatmap::HeatmapAnnotation(
        df = data.frame(Group = factor(group_values, levels = group_levels)),
        col = list(Group = palette)
      )
      
      p1 <- ComplexHeatmap::Heatmap(
        expr,
        name = "log2(norm counts)",
        top_annotation = ha,
        cluster_rows = TRUE,
        cluster_columns = cluster_cols,
        show_column_names = FALSE,
        show_row_names = TRUE,
        row_names_gp = grid::gpar(fontsize = 6, fontface = "bold"),
        column_title = paste("Top", top_n, "Diff Genes"),
        column_title_gp = grid::gpar(fontface = "bold")
      )
      ComplexHeatmap::draw(p1)
    })
    
    # Heatmap (PDF)
    output$download_heatmap <- downloadHandler(
      filename = function() paste0("heatmap_", Sys.Date(), ".pdf"),
      content = function(file) {
        grDevices::pdf(file, width = 10, height = 8)
        req(filtered_data_rv, res_reactive())
        
        top_n <- input$num_genes
        meta_col <- input$metadata_column
        cluster_cols <- isTRUE(input$cluster_columns)
        
        top_genes <- head(res_reactive()[order(res_reactive()$padj), ], top_n)
        if (nrow(top_genes) < 2) {
          grid::grid.newpage(); grid::grid.text("Not enough DE genes for heatmap."); grDevices::dev.off(); return()
        }
        
        expr <- filtered_data_rv$norm_counts[rownames(top_genes), , drop = FALSE]
        expr <- log2(expr + 1)
        
        group_values <- as.character(filtered_data_rv$samples[[meta_col]])
        group_levels <- unique(group_values)
        
        palette <- if (length(group_levels) == 1) {
          setNames("#1f78b4", group_levels)
        } else if (length(group_levels) == 2) {
          setNames(RColorBrewer::brewer.pal(3, "Dark2")[1:2], group_levels)
        } else if (length(group_levels) <= 8) {
          setNames(RColorBrewer::brewer.pal(length(group_levels), "Dark2"), group_levels)
        } else {
          setNames(colorspace::rainbow_hcl(length(group_levels)), group_levels)
        }
        
        ha <- ComplexHeatmap::HeatmapAnnotation(
          df = data.frame(Group = factor(group_values, levels = group_levels)),
          col = list(Group = palette)
        )
        
        hp <- ComplexHeatmap::Heatmap(
          expr,
          name = "log2(norm counts)",
          top_annotation = ha,
          cluster_rows = TRUE,
          cluster_columns = cluster_cols,
          show_column_names = FALSE,
          show_row_names = TRUE,
          row_names_gp = grid::gpar(fontsize = 6, fontface = "bold"),
          column_title = paste("Top", top_n, "Diff Genes"),
          column_title_gp = grid::gpar(fontface = "bold")
        )
        ComplexHeatmap::draw(hp)
        grDevices::dev.off()
      }
    )
    
    # Results CSV
    output$download_de_table <- downloadHandler(
      filename = function() "differential_expression_results.csv",
      content = function(file) utils::write.csv(res_reactive(), file, row.names = FALSE)
    )
    
    # Expose the selected contrast to other modules
    return(list(
      group_var  = shiny::reactive(input$metadata_column),
      ref_level  = shiny::reactive(input$reference_condition),
      test_level = shiny::reactive(input$test_condition)
    ))
  })
}
