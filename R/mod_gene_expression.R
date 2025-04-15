# === Module: mod_gene_expression_plot ===

#' Module to visualize individual gene expression across sample groups
#'
#' This module creates a boxplot (with jittered points) of log2 normalized gene expression
#' values for a user-input gene or genes, grouped by a selected metadata variable.
#'
#' @param input Shiny input object
#' @param output Shiny output object
#' @param session Shiny session object
#' @param data ReactiveValues containing counts, samples, norm_counts, and species
#' @return None. Outputs are rendered to the UI
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_jitter labs theme_minimal scale_color_brewer scale_y_continuous
#' @importFrom stringr str_split
#' @importFrom DT datatable
#' @importFrom utils head read.csv write.csv str
#' @importFrom stats as.formula dist model.matrix prcomp quantile relevel var na.omit cor
#' @importFrom grDevices dev.off pdf colorRampPalette
#' @importFrom grid gpar
#' @importFrom shiny isolate req renderPlot renderImage showNotification downloadHandler renderPrint
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom shinythemes shinytheme
#' @importFrom tidyr pivot_longer
#' @importFrom tidyselect all_of
#' @importFrom shiny validate need
#' @export
mod_gene_expression_plot <- function(input, output, session, data) {

  # Update group selections based on sample metadata
  observe({
    req(data$samples)
    choices <- colnames(data$samples)
    updateSelectInput(session, "group_select_geneexpr", choices = choices)
  })

  # Render gene expression plot
  output$geneExpressionPlot <- renderPlot({
    req(input$gene_select, input$group_select_geneexpr)
    req(data$norm_counts, data$samples)
# Parse and clean user input (space/comma-separated gene names)
genes_raw <- input$gene_select
# DEBUGGING: Check input string for any issues
#print("Raw input gene string:")
#print(genes_raw)

# Split the input string by spaces or commas
selected_genes <- unlist(stringr::str_split(genes_raw, "[\\s,]+"))

# DEBUGGING: Check after splitting the string
#print("Selected genes after splitting:")
#print(selected_genes)

# Remove any extra spaces (trim) and empty entries
selected_genes <- trimws(selected_genes)  # Remove leading and trailing spaces
selected_genes <- selected_genes[selected_genes != ""]  # Remove empty strings

# DEBUGGING: Print the final list of selected genes
#print("Final list of selected genes:")
#print(selected_genes)
# Available genes from data$norm_counts
available_genes <- rownames(data$norm_counts)
#print("Available genes:")
#print(head(available_genes))

# Check if available genes are Ensembl IDs and selected genes are symbols
if (is_ensembl_id(available_genes) && !is_ensembl_id(selected_genes)) {
  # If available genes are Ensembl and selected genes are symbols, convert symbols to Ensembl IDs
  symbol_to_ens <- convert_ensembl_to_symbol(available_genes, data$species)
  
  # Map selected gene symbols to Ensembl IDs
  symbol_map <- names(symbol_to_ens)[match(selected_genes, symbol_to_ens)]
  
  # Find the matched Ensembl genes in available genes
  found_genes <- symbol_map[!is.na(symbol_map) & symbol_map %in% available_genes]
  
  # Create labels for the found genes
  gene_labels <- selected_genes[!is.na(symbol_map) & symbol_map %in% available_genes]
  names(gene_labels) <- found_genes
} else {
  # If both are in symbol format, directly match the symbols
  found_genes <- selected_genes[selected_genes %in% available_genes]
  
  # Debugging: Check the found genes
  #print("Found genes:")
  #print(found_genes)
  
  # If no valid genes are found, show an error and stop execution
  if (length(found_genes) == 0) {
    showNotification("None of the entered genes matched the dataset.", type = "error")
    return(NULL)
  }
  
  # Assign labels to the found genes
  gene_labels <- setNames(found_genes, found_genes)
}

# Prepare expression data for plotting
  df <- data.frame(t(data$norm_counts[found_genes, , drop = FALSE]))
  df$Sample <- rownames(df)
  df <- merge(df, data$samples, by.x = "Sample", by.y = "row.names")

# Convert to long format
   long_df <- tidyr::pivot_longer(
    df, cols = all_of(found_genes), names_to = "Gene", values_to = "Expression"
   )

# Annotate
   long_df$Group <- df[[input$group_select_geneexpr]][match(long_df$Sample, df$Sample)]
   long_df$Expression <- as.numeric(long_df$Expression)
   long_df$log2_Expression <- log2(long_df$Expression + 1)

# Label genes
   long_df$Gene <- gene_labels[long_df$Gene]
   long_df$Gene <- factor(long_df$Gene, levels = unique(long_df$Gene))

# Plot
   grid::grid.newpage()
   grid::pushViewport(grid::viewport())

    ggplot(long_df, aes(x = Group, y = log2_Expression, color = Group)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.25, alpha = 0.6) +
      facet_wrap(~Gene, scales = "free_y") +
      theme_minimal() +
      theme(axis.title = element_text(face = "bold")) +
      labs(title = "Individual Gene Expression by Group", x = "Group", y = "Log2 Normalized Expression") +
      scale_color_brewer(palette = "Dark2")
  })

  # Download gene expression plot
  output$download_gene_plot <- downloadHandler(
    filename = function() { paste0("gene_expression_plot.pdf") },
    content = function(file) {
      pdf(file)
      p <- isolate(output$geneExpressionPlot())
      if (!is.null(p)) print(p)
      dev.off()
    } 
 )
}

# === Register in server ===
# mod_gene_expression_plot(input, output, session, data)

