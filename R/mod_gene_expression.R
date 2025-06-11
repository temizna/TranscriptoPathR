# === Module: mod_gene_expression_plot ===

#' Module to visualize individual gene expression across sample groups
#'
#' This module creates a boxplot (with jittered points) of log2 normalized gene expression
#' values for a user-input gene or genes, grouped by a selected metadata variable.
#'
#' @param input Shiny input object
#' @param output Shiny output object
#' @param session Shiny session object
#' @param filtered_data_rv reactive list containing counts, samples, norm_counts, and species
#' @return None. Outputs are rendered to the UI
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_jitter labs theme_minimal scale_color_brewer scale_y_continuous
#' @importFrom stringr str_split
#' @importFrom DT datatable
#' @importFrom utils head read.csv write.csv str
#' @importFrom stats as.formula dist model.matrix prcomp quantile relevel var na.omit cor
#' @importFrom grDevices dev.off pdf colorRampPalette
#' @importFrom grid gpar
#' @importFrom shiny isolate req renderPlot updateSelectInput renderImage showNotification downloadHandler renderPrint
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom shinythemes shinytheme
#' @importFrom tidyr pivot_longer
#' @importFrom tidyselect all_of
#' @importFrom shiny validate need
#' @importFrom utils install.packages
#' @export
mod_gene_expression_plot <- function(input, output, session, filtered_data_rv) {
  
  # Update group selections based on sample metadata
  observe({
    req(filtered_data_rv$samples)
    choices <- colnames(filtered_data_rv$samples)
    updateSelectInput(session, "group_select_geneexpr", choices = choices)
  })
  generate_gene_expression_plot <- function(filtered_data, selected_genes, group_col, species) {
    available_genes <- rownames(filtered_data$norm_counts)
    
    if (is_ensembl_id(available_genes)) {
      # Convert user-given symbols to Ensembl
      symbol_to_ens <- convert_symbol_to_ensembl(selected_genes, species)
      matched_symbols <- names(symbol_to_ens)
      found_genes <- unname(symbol_to_ens)
      found_genes <- found_genes[!is.na(found_genes) & found_genes %in% available_genes]
      matched_symbols <- matched_symbols[!is.na(found_genes)]
      names(matched_symbols) <- found_genes  # Ensembl -> Symbol
      gene_labels <- matched_symbols
    } else {
      found_genes <- selected_genes[selected_genes %in% available_genes]
      gene_labels <- setNames(found_genes, found_genes)  # Symbol -> Symbol
    }
    expr_mat <- filtered_data$norm_counts[found_genes, , drop = FALSE]
    expr_df <- as.data.frame(t(expr_mat))
    expr_df$Sample <- rownames(expr_df)
    
    # Join with sample metadata
    df <- cbind(expr_df, filtered_data$samples[rownames(expr_df), , drop = FALSE])
    # Convert long format
    long_df <- tidyr::pivot_longer(
      df, cols = all_of(found_genes), names_to = "Gene", values_to = "Expression"
    )
    
    # Add group and clean labels
    long_df$Group <- df[[group_col]][match(long_df$Sample, df$Sample)]
    long_df$Expression <- as.numeric(long_df$Expression)
    long_df$log2_Expression <- log2(long_df$Expression + 1)
    
    # Replace Gene IDs with user input labels if possible
    long_df$Gene <- gene_labels[long_df$Gene]
    long_df$Gene[is.na(long_df$Gene)] <- long_df$Gene[is.na(long_df$Gene)]
    long_df$Gene <- factor(long_df$Gene, levels = unique(long_df$Gene))
    
    
    p <- ggplot(long_df, aes(x = Group, y = log2_Expression, color = Group)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.25, alpha = 0.6) +
      facet_wrap(~Gene, scales = "fixed") +
      theme_minimal() +
      theme(
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(face = "bold", angle = 90, hjust = 1)
      ) +
      labs(
        title = "Individual Gene Expression by Group",
        x = "Group", y = "Log2 Normalized Expression"
      ) +
      scale_color_brewer(palette = "Paired")
    
    return(p)
  }
  
  compute_anova_table <- function(filtered_data, selected_genes, group_col, species) {
    available_genes <- rownames(filtered_data$norm_counts)
    
    if (is_ensembl_id(available_genes)) {
      symbol_to_ens <- convert_symbol_to_ensembl(selected_genes, species)
      matched_symbols <- names(symbol_to_ens)
      found_genes <- unname(symbol_to_ens)
      found_genes <- found_genes[!is.na(found_genes) & found_genes %in% available_genes]
      matched_symbols <- matched_symbols[!is.na(found_genes)]
      names(matched_symbols) <- found_genes
      gene_labels <- matched_symbols
    } else {
      found_genes <- selected_genes[selected_genes %in% available_genes]
      gene_labels <- setNames(found_genes, found_genes)
    }
    
    expr_mat <- filtered_data$norm_counts[found_genes, , drop = FALSE]
    expr_df <- as.data.frame(t(expr_mat))
    expr_df$Sample <- rownames(expr_df)
    df <- cbind(expr_df, filtered_data$samples[rownames(expr_df), , drop = FALSE])
    
    long_df <- tidyr::pivot_longer(
      df, cols = all_of(found_genes), names_to = "Gene", values_to = "Expression"
    )
    long_df$Group <- df[[group_col]][match(long_df$Sample, df$Sample)]
    long_df$Expression <- as.numeric(long_df$Expression)
    long_df$log2_Expression <- log2(long_df$Expression + 1)
    long_df$Gene <- gene_labels[long_df$Gene]
    long_df$Gene[is.na(long_df$Gene)] <- long_df$Gene[is.na(long_df$Gene)]
    
    all_results <- lapply(unique(long_df$Gene), function(gene) {
      gene_data <- subset(long_df, Gene == gene)
      fit <- aov(log2_Expression ~ Group, data = gene_data)
      anova_p <- summary(fit)[[1]]["Group", "Pr(>F)"]
      tukey <- tryCatch(TukeyHSD(fit, "Group"), error = function(e) NULL)
      if (is.null(tukey)) return(NULL)
      
      pairwise <- as.data.frame(tukey$Group)
      pairwise$Comparison <- rownames(pairwise)
      pairwise$Gene <- gene
      pairwise$FDR <- p.adjust(pairwise$`p adj`, method = "fdr")
      pairwise[, c("Gene", "Comparison", "diff", "lwr", "upr", "p adj", "FDR")]
    })
    
    result_df <- do.call(rbind, all_results)
    colnames(result_df) <- c("Gene", "Comparison", "Difference", "Lower CI", "Upper CI", "P.Value", "FDR")
    result_df$P.Value <- signif(result_df$P.Value, 4)
    result_df$FDR <- signif(result_df$FDR, 4)
    result_df
  }
  
  
  
  output$geneExpressionPlot <- renderPlot({
    req(input$gene_select, input$group_select_geneexpr)
    req(filtered_data_rv$species, filtered_data_rv)
    
    # Parse and clean user input (space/comma-separated gene names)
    genes_raw <- input$gene_select
    selected_genes <- unlist(stringr::str_split(genes_raw, "[\\s,]+"))
    selected_genes <- trimws(selected_genes)
    selected_genes <- selected_genes[selected_genes != ""]
    
    plot <- generate_gene_expression_plot(
      filtered_data = filtered_data_rv,
      selected_genes = selected_genes,
      group_col = input$group_select_geneexpr,
      species = filtered_data_rv$species
    )
    if (is.null(plot)) return(NULL)
    plot
  })
  output$geneExpressionStats <- DT::renderDT({
    req(input$gene_select, input$group_select_geneexpr, filtered_data_rv$species)
    
    genes_raw <- input$gene_select
    selected_genes <- unlist(stringr::str_split(genes_raw, "[\\s,]+"))
    selected_genes <- trimws(selected_genes)
    selected_genes <- selected_genes[selected_genes != ""]
    
    df <- compute_anova_table(
      filtered_data = filtered_data_rv,
      selected_genes = selected_genes,
      group_col = input$group_select_geneexpr,
      species = filtered_data_rv$species
    )
    
    # Format all numeric columns for readability
    numeric_cols <- sapply(df, is.numeric)
    df[numeric_cols] <- lapply(df[numeric_cols], function(col) {
      ifelse(abs(col) < 0.001, formatC(col, format = "e", digits = 2), round(col, 3))
    })
    
    DT::datatable(
      df,
      options = list(pageLength = 10, autoWidth = TRUE),
      rownames = FALSE
    ) %>%
      DT::formatStyle(
        "FDR",
        backgroundColor = DT::styleInterval(0.05, c("#ffcccc", NA)),
        fontWeight = DT::styleInterval(0.05, c("bold", "normal"))
      ) %>%
      DT::formatStyle(
        "Gene",
        fontWeight = "bold"
      )
  })
  
  output$download_gene_stats <- downloadHandler(
    filename = function() {
      paste0("gene_expression_anova_", Sys.Date(), ".csv")
    },
    content = function(file) {
      genes_raw <- input$gene_select
      selected_genes <- unlist(stringr::str_split(genes_raw, "[\\s,]+"))
      selected_genes <- trimws(selected_genes)
      selected_genes <- selected_genes[selected_genes != ""]
      
      df <- compute_anova_table(
        filtered_data = filtered_data_rv,
        selected_genes = selected_genes,
        group_col = input$group_select_geneexpr,
        species = filtered_data_rv$species
      )
      
      # Apply same formatting as DT table
      numeric_cols <- sapply(df, is.numeric)
      df[numeric_cols] <- lapply(df[numeric_cols], function(col) {
        ifelse(abs(col) < 0.001, formatC(col, format = "e", digits = 2), round(col, 3))
      })
      
      write.csv(df, file, row.names = FALSE)
    }
  )
  
  
  # Download gene expression plot
  output$download_gene_plot <- downloadHandler(
    filename = function() {
      paste0("gene_expression.plot.pdf")  # You can change the file extension as needed (e.g., ".pdf")
    },
    content = function(file) {
      # Generate the plot using the generate_gene_expression_plot function
    
      genes_raw <- input$gene_select
      selected_genes <- unlist(stringr::str_split(genes_raw, "[\\s,]+"))
      selected_genes <- trimws(selected_genes)
      selected_genes <- selected_genes[selected_genes != ""]
      
      plot <- generate_gene_expression_plot(
        filtered_data = filtered_data_rv,
        selected_genes = selected_genes,
        group_col = input$group_select_geneexpr,
        species = filtered_data_rv$species
      )
      if (is.null(plot)) {
        showNotification("No valid genes found for the plot.", type = "error")
        return(NULL)
      }
      ggsave(file, plot, device = "pdf", width = 10, height = 8) # Save as PNG (change format as needed)
    }
  )
  
  
  
}

# === Register in server ===
# mod_gene_expression_plot(input, output, session, data)

