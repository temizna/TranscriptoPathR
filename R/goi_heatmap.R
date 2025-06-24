#' Genes of Interest Heatmap Module
#'
#' Allows the user to upload a gene list and visualize expression as a heatmap
#' using `filtered_data_rv$norm_counts`.
#'
#' @param input Shiny input object
#' @param output Shiny output object
#' @param session Shiny session object
#' @param filtered_data_rv A reactiveValues list containing:
#'   \itemize{
#'     \item \code{norm_counts}: normalized gene expression matrix (genes x samples)
#'     \item \code{samples}: metadata for samples
#'     \item \code{species}: species name
#'   }
#'
#' @importFrom shiny req fileInput renderPlot downloadHandler NS tagList
#' @importFrom utils read.csv write.csv
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation draw
#' @importFrom grDevices pdf dev.off
#' @importFrom grid gpar
#' @importFrom colorspace rainbow_hcl 
#' @importFrom RColorBrewer brewer.pal
#' @export
mod_goi_heatmap <- function(input, output, session, filtered_data_rv) {
    goi_expr_rv <- reactiveVal(NULL)
    goi_ha_rv <- reactiveVal(NULL)
    
    observe({
      req(filtered_data_rv$samples)
      updateSelectInput(session, "goi_group_column", choices = colnames(filtered_data_rv$samples))
    })
    
    observeEvent(input$goi_file, {
      req(input$goi_file, filtered_data_rv$norm_counts)
      gene_file <- input$goi_file$datapath
      gene_list <- tryCatch({
        ext <- tools::file_ext(input$goi_file$name)
        if (tolower(ext) == "csv") {
          read.csv(gene_file, header = TRUE)[[1]]
        } else if (tolower(ext) == "txt") {
          read.table(gene_file, header = FALSE, stringsAsFactors = FALSE)[[1]]
        } else stop("Unsupported file format.")
      }, error = function(e) {
        showNotification(paste("Failed to read gene list:", e$message), type = "error")
        return(NULL)
      })
      
      selected_genes <- trimws(gene_list)
      available_genes <- rownames(filtered_data_rv$norm_counts)
      species <- filtered_data_rv$species
      
      if (is_ensembl_id(available_genes)) {
        symbol_to_ens <- convert_symbol_to_ensembl(selected_genes, species)
        found_ens <- unname(symbol_to_ens)
        matched_symbols <- names(symbol_to_ens)
        
        valid_idx <- !is.na(found_ens) & found_ens %in% available_genes
        found_genes <- found_ens[valid_idx]
        matched_symbols <- matched_symbols[valid_idx]
        names(found_genes) <- matched_symbols
        
        expr <- log2(filtered_data_rv$norm_counts[found_genes, , drop = FALSE] + 1)
        rownames(expr) <- names(found_genes)
      } else {
        found_genes <- selected_genes[selected_genes %in% available_genes]
        expr <- log2(filtered_data_rv$norm_counts[found_genes, , drop = FALSE] + 1)
        rownames(expr) <- found_genes
      }
      
      group_col <- input$goi_group_column
      group_values <- as.character(filtered_data_rv$samples[[group_col]])
      group_levels <- unique(group_values)
      if (length(group_levels) == 1) {
        palette <- setNames("#1f78b4", group_levels)
      } else if (length(group_levels) == 2) {
        palette <- setNames(RColorBrewer::brewer.pal(3, "Dark2")[1:2], group_levels)
      } else if (length(group_levels) <= 8) {
        palette <- setNames(RColorBrewer::brewer.pal(length(group_levels), "Dark2"), group_levels)
      } else {
        palette <- setNames(colorspace::rainbow_hcl(length(group_levels)), group_levels)
      }
      
      names(palette) <- group_levels
      
      ha <- ComplexHeatmap::HeatmapAnnotation(
        df = data.frame(Group = factor(group_values, levels = group_levels)),
        col = list(Group = palette)
      )
      
      goi_expr_rv(expr)
      goi_ha_rv(ha)
    })
    
    output$goi_heatmap <- renderPlot({
      req(goi_expr_rv(), goi_ha_rv())
      ComplexHeatmap::draw(ComplexHeatmap::Heatmap(
        goi_expr_rv(),
        name = "log2(norm counts)",
        top_annotation = goi_ha_rv(),
        cluster_rows = TRUE,
        cluster_columns = input$cluster_columns,
        show_column_names = FALSE,
        show_row_names = TRUE,
        row_names_gp = grid::gpar(fontsize = 6, fontface = "bold"),
        column_title = "Genes of Interest",
        column_title_gp = grid::gpar(fontface = "bold", fontsize = 6)
      ))
    })
    
    output$download_goi_heatmap <- downloadHandler(
      filename = function() { paste0("GOI_heatmap_", Sys.Date(), ".pdf") },
      content = function(file) {
        req(goi_expr_rv(), goi_ha_rv())
        pdf(file)
        ComplexHeatmap::draw(ComplexHeatmap::Heatmap(
          goi_expr_rv(),
          name = "log2(norm counts)",
          top_annotation = goi_ha_rv(),
          cluster_rows = TRUE,
          cluster_columns = input$cluster_columns,
          show_column_names = FALSE,
          show_row_names = TRUE,
          row_names_gp = grid::gpar(fontsize = 4, fontface = "bold"),
          column_title = "Genes of Interest",
          column_title_gp = grid::gpar(fontface = "bold", fontsize = 6)
        ))
        dev.off()
      }
    )
    observeEvent(input$goi_group_column, {
      req(goi_expr_rv(), filtered_data_rv$samples)
      
      group_col <- input$goi_group_column
      group_values <- as.character(filtered_data_rv$samples[[group_col]])
      group_levels <- unique(group_values)
      
      if (length(group_levels) == 1) {
        palette <- setNames("#1f78b4", group_levels)
      } else if (length(group_levels) == 2) {
        palette <- setNames(RColorBrewer::brewer.pal(3, "Dark2")[1:2], group_levels)
      } else if (length(group_levels) <= 8) {
        palette <- setNames(RColorBrewer::brewer.pal(length(group_levels), "Dark2"), group_levels)
      } else {
        palette <- setNames(colorspace::rainbow_hcl(length(group_levels)), group_levels)
      }
      
      names(palette) <- group_levels
      
      ha <- ComplexHeatmap::HeatmapAnnotation(
        df = data.frame(Group = factor(group_values, levels = group_levels)),
        col = list(Group = palette)
      )
      
      goi_ha_rv(ha)
    })
    
  }
  