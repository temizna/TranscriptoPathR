# === Module: mod_qc_plots ===

#' Quality Control Plot Module
#'
#' Generates various QC plots such as PCA, sample distance heatmap, mean-variance plots,
#' and gene expression variance histograms from normalized RNA-seq count data.
#'
#' @param input Shiny input object
#' @param output Shiny output object
#' @param session Shiny session object
#' @param data ReactiveValues containing counts, samples, norm_counts, and species
#' @return None. Outputs are rendered to UI.
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_jitter labs theme_minimal scale_color_brewer scale_y_continuous
#' @importFrom stringr str_split
#' @importFrom DT datatable
#' @importFrom utils head read.csv write.csv str
#' @importFrom stats as.formula dist model.matrix prcomp quantile relevel var cor na.omit
#' @importFrom grDevices dev.off pdf colorRampPalette
#' @importFrom grid gpar
#' @importFrom shiny isolate req renderPlot renderImage showNotification reactiveVal downloadHandler
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom shinythemes shinytheme
#' @importFrom org.Mm.eg.db org.Mm.eg.db
#' @export
mod_qc_plot <- function(input, output, session, data) {

  # Update QC group selection dropdown based on design metadata
  observe({
    req(data$samples)
    updateSelectInput(session, "group_select_qc", choices = colnames(data$samples))
  })

  output$qcPlot <- renderPlot({
    req(input$qc_plot_type, data$samples, data$norm_counts)
    group_col <- input$group_select_qc

    if (!group_col %in% colnames(data$samples)) {
      showNotification("Please select a valid grouping column for QC.", type = "error")
      return(NULL)
    }

    group_factor <- factor(data$samples[[group_col]])

    if (input$qc_plot_type == "PCA") {
      expr <- log2(data$norm_counts + 1)
      expr <- expr[apply(expr, 1, function(x) var(x, na.rm = TRUE) > 0), , drop = FALSE]
      pca <- prcomp(t(expr), scale. = TRUE)
      percentVar <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)
      df <- data.frame(
        PC1 = pca$x[, 1],
        PC2 = pca$x[, 2],
        Group = group_factor,
        Sample = rownames(data$samples)
      )
      
      ggplot(df, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
        geom_point(size = 3) +
        #geom_text_repel(show.legend = FALSE, max.overlaps = if (input$show_labels) 100 else 0) +
        xlab(paste0("PC1 (", percentVar[1], "%)")) +
        ylab(paste0("PC2 (", percentVar[2], "%)")) +
        theme_minimal() +
        theme(axis.title = element_text(face = "bold")) +
        scale_color_brewer(palette = "Set2")
    } else if (input$qc_plot_type == "Sample Distance") {
      dist_matrix <- dist(t(log2(data$norm_counts + 1)))
      mat <- as.matrix(dist_matrix)
      rownames(mat) <- colnames(data$norm_counts)
      colnames(mat) <- colnames(data$norm_counts)
      ComplexHeatmap::Heatmap(mat, name = "Distance",  column_names_gp = grid::gpar(fontsize = 4),
                              width = unit(20, "cm"), height = unit(20, "cm"))

    } else if (input$qc_plot_type == "Mean-Variance") {
      means <- rowMeans(data$norm_counts)
      vars <- apply(data$norm_counts, 1, var)
      df <- data.frame(Mean = means, Variance = vars)
      ggplot(df, aes(x = Mean, y = Variance)) +
        geom_point(alpha = 0.5) +
        scale_x_log10() +
        scale_y_log10() +
        theme_minimal() +
        theme(axis.title = element_text(face = "bold")) +
        labs(title = "Mean-Variance Plot", x = "Mean Expression", y = "Variance")

    } else if (input$qc_plot_type == "Variance Histogram") {
      gene_vars <- apply(data$norm_counts, 1, var)
      df <- data.frame(Variance = gene_vars, Group = rep(group_factor, each = length(gene_vars)))
      ggplot(df, aes(x = Variance, fill = Group)) +
        geom_histogram(bins = 100, position = "identity", alpha = 0.6, color = "black") +
        xlim(0, quantile(gene_vars, 0.99)) +
        theme_minimal() +
        theme(axis.title = element_text(face = "bold")) +
        labs(title = "Variance Histogram by Group", x = "Variance", y = "Gene Count") +
        scale_fill_brewer(palette = "Set1")
    }
  })

  # Download handler for QC plot
  output$download_qc_plot <- downloadHandler(
    filename = function() paste0(input$qc_plot_filename),
    content = function(file) {
      pdf(file)
      p <- isolate(output$qcPlot())
      if (!is.null(p)) print(p)
      dev.off()
    }
  )
}

# === Register in server ===
# mod_qc_plots(input, output, session, data)

