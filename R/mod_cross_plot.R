# === Module: mod_cross_plot ===

#' Cross Plot and Venn Diagram Module
#'
#' @param input Shiny input object.
#' @param output Shiny output object.
#' @param session Shiny session object.
#' @param dds_rv A reactive value containing the DESeq2 object.
#' @param data Data containing counts, samples, norm_counts, and species
#'
#' @return None. Updates Cross Plot and Venn Diagram and download handlers.
#' @importFrom ggplot2 ggplot aes geom_point geom_text geom_vline geom_hline theme_minimal theme
#' @importFrom ggrepel geom_label_repel
#' @importFrom dplyr inner_join distinct
#' @importFrom VennDiagram draw.pairwise.venn
#' @importFrom utils head read.csv write.csv str
#' @importFrom stats as.formula dist model.matrix prcomp quantile relevel var cor
#' @importFrom grDevices dev.off pdf colorRampPalette
#' @importFrom grid gpar
#' @importFrom shiny isolate req renderPlot renderImage showNotification downloadHandler renderPrint
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom shinythemes shinytheme
#' @importFrom gridExtra grid.arrange
#' @export

mod_cross_plot <- function(input, output, session, data, dds_rv) {
  
  # Reactive expression to store the calculated contrast data for cross plot
  crossplot_data <- reactiveVal()
  
  # Ensure that the dropdowns are populated once the data is available
  observe({
    req(data$samples)
    cols <- colnames(data$samples)
    updateSelectInput(session, "metadata_column_x", choices = cols)
    updateSelectInput(session, "metadata_column_y", choices = cols)
  })
  
  # Update the test and reference conditions based on the selected metadata columns
  observeEvent(input$metadata_column_x, {
    req(data$samples)
    col <- input$metadata_column_x
    levels <- unique(as.character(data$samples[[col]]))
    updateSelectInput(session, "reference_condition_x", choices = levels)
    updateSelectInput(session, "test_condition_x", choices = levels)
  })
  
  observeEvent(input$metadata_column_y, {
    req(data$samples)
    col <- input$metadata_column_y
    levels <- unique(as.character(data$samples[[col]]))
    updateSelectInput(session, "reference_condition_y", choices = levels)
    updateSelectInput(session, "test_condition_y", choices = levels)
  })
  # Observing when the cross plot button is clicked
  observeEvent(input$run_crossplot, {
    req(dds_rv())  # Ensure the DESeq2 results object is available
    # Ensure that required inputs are selected
    req(input$metadata_column_x, input$test_condition_x, input$reference_condition_x,
        input$metadata_column_y, input$test_condition_y, input$reference_condition_y)
    
    # Create contrasts for X and Y axes based on user inputs
    x_comp <- c(input$metadata_column_x, input$test_condition_x, input$reference_condition_x)
    y_comp <- c(input$metadata_column_y, input$test_condition_y, input$reference_condition_y)
    
    # Check if results for the X and Y comparisons already exist
    contrast_exists_x <- paste(input$metadata_column_x, input$test_condition_x, input$reference_condition_x, sep = "_") %in% resultsNames(dds_rv())
    contrast_exists_y <- paste(input$metadata_column_y, input$test_condition_y, input$reference_condition_y, sep = "_") %in% resultsNames(dds_rv())
    
    # If either contrast doesn't exist, run DESeq2 analysis for the respective contrast
    if (!contrast_exists_x) {
      # Run DESeq2 for the X contrast
      design_formula_x <- as.formula(paste("~", input$metadata_column_x))
      dds_x <- DESeq2::DESeqDataSetFromMatrix(countData = data$counts, colData = data$samples, design = design_formula_x)
      dds_x <- DESeq2::DESeq(dds_x)  # Run DESeq2 analysis for X contrast
    } else { 
      dds_x<-dds_rv()
      }
    if (!contrast_exists_y) {
      # Run DESeq2 for the Y contrast
      design_formula_y <- as.formula(paste("~", input$metadata_column_y))
      dds_y <- DESeq2::DESeqDataSetFromMatrix(countData = data$counts, colData = data$samples, design = design_formula_y)
      dds_y <- DESeq2::DESeq(dds_y)  # Run DESeq2 analysis for Y contrast
    } else {
      dds_y<-dds_rv()
    }
    
    # Extract DESeq2 results for X and Y comparisons
    res_x_data <- results(dds_x, contrast = x_comp)
    res_y_data <- results(dds_y, contrast = y_comp)
    
    # Store results in reactive values
    crossplot_data(data.frame(
      gene = rownames(res_x_data),
      log2FoldChange_x = res_x_data$log2FoldChange,
      padj_x = res_x_data$padj,
      log2FoldChange_y = res_y_data$log2FoldChange,
      padj_y = res_y_data$padj
    ))  # Store merged data for later plotting
      
    df <- crossplot_data()
    # Top N significant genes by combined adjusted p-values
    df$combined_padj <- pmax(df$padj_x, df$padj_y, na.rm = TRUE)
    top_genes <- head(df[order(df$combined_padj), "gene"], input$crossplot_topgenes)
    #print(top_genes)
    # Process highlight genes from input
    highlight_genes<-input$crossplot_gene_label
    highlight_genes <- if (nchar(trimws(highlight_genes)) > 0) {
      unlist(stringr::str_split(highlight_genes, "[\\s,]+"))
    } else {
      character(0)
    }
    if (length(highlight_genes) == 1 && highlight_genes == "") highlight_genes <- character(0)
    # Combine top genes and highlighted genes
    label_genes <- union(top_genes, highlight_genes)
    
    # Label genes based on this combined list
    df$label <- ifelse(df$gene %in% label_genes, df$gene, NA)
    
    # Update crossplot data
    crossplot_data(df)  # Update the reactive value with the new data
    print("Updated Data with Labels:")
    print(head(df))
   
  })
  
  # Observe the changes and trigger updates in the reactive data
  #observeEvent(
  #  list(input$crossplot_topgenes, input$crossplot_gene_label), 
  #  {
  #    
  #  })
  
  # Crossplot rendering
 
  output$crossPlot <- renderPlot({
    req(crossplot_data())  # Ensure data is available
    df <- crossplot_data()
    # Top N significant genes by combined adjusted p-values
    df$combined_padj <- pmax(df$padj_x, df$padj_y, na.rm = TRUE)
    # Define categories based on log2FoldChange for X and Y comparisons
    df$category <- "Other"
    df$category[df$log2FoldChange_x >= 1 & df$log2FoldChange_y >= 1] <- "Up-Up"
    df$category[df$log2FoldChange_x <= -1 & df$log2FoldChange_y <= -1] <- "Down-Down"
    df$category[df$log2FoldChange_x > 1 & abs(df$log2FoldChange_y) <= 1] <- "Comp1-only"
    df$category[abs(df$log2FoldChange_x) <= 1 & df$log2FoldChange_y >= 1] <- "Comp2-only"
    
    # Limit the dataset to the selected number of genes for plotting
    df <- head(df, input$crossplot_gene_count)
    #print(head(df))
    # Calculate correlation values (Pearson and Spearman)
    pearson_r <- round(cor(df$log2FoldChange_x, df$log2FoldChange_y, method = 'pearson', use = 'complete.obs'), 3)
    spearman_rho <- round(cor(df$log2FoldChange_x, df$log2FoldChange_y, method = 'spearman', use = 'complete.obs'), 3)
    
    # Create the cross plot
    ggplot(df, aes(x = log2FoldChange_x, y = log2FoldChange_y, color = category)) +
      annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -1.1,
               label = paste0("Pearson r = ", pearson_r, "   Spearman r = ", spearman_rho),
               size = 4) +
      geom_point(alpha = 0.7) +
      geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "black") +
      geom_vline(xintercept = c(-1, 1), linetype = "solid", color = "blue") +
      geom_hline(yintercept = c(-1, 1), linetype = "solid", color = "blue") +
      ggrepel::geom_text_repel(data = subset(df, !is.na(label)), aes(label = label), max.overlaps = Inf, size = 5) +
      scale_color_manual(values = c(
        "Up-Up" = "firebrick",
        "Down-Down" = "royalblue",
        "Comp1-only" = "darkorange",
        "Comp2-only" = "darkgreen",
        "Other" = "gray"
      )) +
      theme_minimal() +
      labs(title = "Cross Plot: log2FC vs log2FC",
           x = paste("log2FC:", input$test_condition_x, "vs", input$reference_condition_x),
           y = paste("log2FC:", input$test_condition_y, "vs", input$reference_condition_y),
           color = "Regulation")
  })
  
  output$crosspathplot <-renderPlot({
    req(crossplot_data(), data$species)  # Ensure data is available
    df <- crossplot_data()
    species <- data$species
    orgdb <- get_orgdb(species)
    # Top N significant genes by combined adjusted p-values
    
    # Define categories based on log2FoldChange for X and Y comparisons
    df$combined_padj <- pmax(df$padj_x, df$padj_y, na.rm = TRUE)
    df$category <- "Other"
    df$category[df$log2FoldChange_x >= 1 & df$log2FoldChange_y >= 1] <- "Up-Up"
    df$category[df$log2FoldChange_x <= -1 & df$log2FoldChange_y <= -1] <- "Down-Down"
    df$category[df$log2FoldChange_x > 1 & abs(df$log2FoldChange_y) <= 1] <- "Comp1-only"
    df$category[abs(df$log2FoldChange_x) <= 1 & df$log2FoldChange_y >= 1] <- "Comp2-only"
    
    # Convert gene symbols to ENTREZ IDs
    if (is_symbol(df$gene)) {
      df_ids <- bitr(df$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = orgdb)
    } else {
      df_ids <- bitr(df$gene, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = orgdb)
    }
    
    # Merge the data with the ENTREZ IDs
    df_merged <- merge(df, df_ids, by.x = "gene", by.y = 1)
    
    # Remove duplicate rows based on ENTREZID
    df_merged <- dplyr::distinct(df_merged, ENTREZID, .keep_all = TRUE)
   
    # Filter rows where log2FoldChange_x is greater than or equal to 1
    df_merged <- df_merged[abs(df_merged$log2FoldChange_x) >= 1, ]
    
    # Select relevant columns
    df_merged2 <- df_merged[, c("ENTREZID", "log2FoldChange_x", "category")]
    
    # Remove rows with NA values in log2FoldChange_x or ENTREZID
    df_merged2 <- df_merged2 %>%
      filter(!is.na(log2FoldChange_x) & !is.na(ENTREZID))
    
    #print(head(df_merged2))
    
    # Perform pathway analysis using compareCluster
    formula_res <- compareCluster(ENTREZID ~ category, data = df_merged2, fun = "enrichKEGG")
    
    # Check if formula_res is NULL or empty
    if (is.null(formula_res) || nrow(formula_res) == 0) {
      showNotification("No pathway enrichment results available.", type = "error")
      return(NULL)
    } else {
      enrichplot::dotplot(formula_res, x = "category")
    }
    
  })
  
  # Cross Venn Plot
  output$crossVennPlot <- renderPlot({
    req(crossplot_data())
    df <- crossplot_data()
    
    # Set the cutoffs for log2FoldChange and adjusted p-value
    lfc_cutoff <- 1
    padj_cutoff <- 0.05
    
    # Identify upregulated and downregulated genes for X and Y comparisons
    up_x <- rownames(df[df$log2FoldChange_x > lfc_cutoff & df$padj_x < padj_cutoff, ])
    up_y <- rownames(df[df$log2FoldChange_y > lfc_cutoff & df$padj_y < padj_cutoff, ])
    down_x <- rownames(df[df$log2FoldChange_x < -lfc_cutoff & df$padj_x < padj_cutoff, ])
    down_y <- rownames(df[df$log2FoldChange_y < -lfc_cutoff & df$padj_y < padj_cutoff, ])
    
    # Create the upregulated genes Venn plot
    up_venn <- VennDiagram::venn.diagram(
      x = list(X_Up = up_x, Y_Up = up_y),
      filename = NULL,
      fill = c("darkorange", "darkgreen"),
      alpha = 0.5,
      main = "Upregulated Genes Venn"
    )
    
    # Create the downregulated genes Venn plot
    down_venn <- VennDiagram::venn.diagram(
      x = list(X_Down = down_x, Y_Down = down_y),
      filename = NULL,
      fill = c("royalblue", "purple"),
      alpha = 0.5,
      main = "Downregulated Genes Venn"
    )
    
    # Turn off the logging (no filename output)
    gridExtra::grid.arrange(
      grid::grid.grabExpr(grid::grid.draw(up_venn)),
      grid::grid.grabExpr(grid::grid.draw(down_venn)),
      ncol = 2
    )
  })
  
  # Download handlers for both plots
  output$download_cross_plot <- downloadHandler(
    filename = function() { "cross_plot.pdf" },
    content = function(file) {
      pdf(file)
      print(isolate(output$crossPlot()))
      dev.off()
    }
  )
  
  output$download_crossVennPlot <- downloadHandler(
    filename = function() { "cross_venn_plot.pdf" },
    content = function(file) {
      pdf(file)
      print(isolate(output$crossVennPlot()))
      dev.off()
    }
  )
}
#mod_cross_plot(input, output, session, dds_rv)

