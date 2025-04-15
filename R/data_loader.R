# === Module: mod_data_upload_design ===
#' Data upload and design
#' @description Load and process uploaded count and design files into DESeq2 objects
#'
#' @param input Shiny input object.
#' @param output Shiny output object.
#' @param session Shiny session object.
#' @param data ReactiveValues to store processed data (counts, samples, norm_counts).
#' @param dds_rv A reactiveVal container for storing the DESeq2 dataset.
#' @return A list with elements: counts, samples, norm_counts, dds
#' @importFrom DESeq2 DESeqDataSetFromMatrix estimateSizeFactors counts
#' @importFrom readxl read_excel
#' @importFrom tibble column_to_rownames
#' @importFrom DT renderDT datatable
#' @importFrom utils head read.csv write.csv str
#' @importFrom stats as.formula dist model.matrix prcomp quantile relevel var
#' @importFrom grDevices dev.off pdf colorRampPalette
#' @importFrom grid gpar
#' @importFrom shiny downloadHandler renderPrint
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom shinythemes shinytheme
#' @importFrom utils install.packages
#' @export
mod_data_upload_design <- function(input, output, session, data, dds_rv) {
load_data_file <- function(file_input) {
  ext <- tools::file_ext(file_input$name)
  if (ext == "csv") {
    read.csv(file_input$datapath, row.names = 1, check.names = FALSE)
  } else if (ext == "xlsx") {
    readxl::read_excel(file_input$datapath, col_names = TRUE) |> tibble::column_to_rownames(1)
  } else {
    stop("Unsupported file format: ", ext)
  }
}
  load_uploaded_data <- function(counts_file, design_file, design_formula) {
    counts_df <- load_data_file(counts_file)
    design_df <- load_data_file(design_file)

    common_samples <- intersect(colnames(counts_df), rownames(design_df))
    if (length(common_samples) < 2) stop("Counts and design files must have matching sample names.")
    counts_df <- counts_df[, common_samples]
    design_df <- design_df[common_samples, ]

    rn <- rownames(design_df)
    design_df <- as.data.frame(lapply(design_df, function(col) {
      if (is.character(col) || is.logical(col)) factor(col) else col
    }))
    rownames(design_df) <- rn

    counts_df <- counts_df[rowSums(counts_df) >= 10, ]
    counts_df <- counts_df[apply(counts_df, 1, var) > 0.1, ]

    formula <- as.formula(design_formula)
    if (!all(all.vars(formula)[-1] %in% colnames(design_df))) {
      stop("Design formula references unknown columns in metadata.")
    }

    dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts_df, colData = design_df, design = formula)
    dds <- DESeq2::estimateSizeFactors(dds)
    norm_counts <- counts(dds, normalized = TRUE)

    list(counts = counts_df, samples = design_df, norm_counts = norm_counts, dds = dds)
  }

  # Observe loading of demo data
  observe({
    if (input$data_mode == "demo") {
      demo_counts <- data.frame(
        row.names = c("GeneA", "GeneB", "GeneC"),
        Sample1 = c(100, 150, 200),
        Sample2 = c(120, 160, 210)
      )
      demo_design <- data.frame(
        row.names = c("Sample1", "Sample2"),
        condition = c("treated", "control")
      )
      data$counts <- demo_counts
      data$samples <- demo_design
      dds <- DESeqDataSetFromMatrix(countData = demo_counts, colData = demo_design, design = ~ condition)
      dds <- estimateSizeFactors(dds)
      dds_rv(dds)
      data$norm_counts <- counts(dds, normalized = TRUE)
    }
  })

  # Observe event for loading user-uploaded data
  observeEvent(input$load_data, {
    req(input$counts_file, input$design_file, input$design_formula)
    tryCatch({
      loaded <- load_uploaded_data(input$counts_file, input$design_file, input$design_formula)
      data$counts <- loaded$counts
      data$samples <- loaded$samples
      data$norm_counts <- loaded$norm_counts
      dds_rv(loaded$dds)
      showNotification("Data successfully loaded!", type = "message")
    }, error = function(e) {
      showNotification(paste("Error loading files:", e$message), type = "error")
    })
  })

  # Provide example template for counts
  output$download_counts_template <- downloadHandler(
    filename = function() { "counts_template.csv" },
    content = function(file) {
      example_counts <- data.frame(
        Gene = c("Gene1", "Gene2"),
        Sample1 = c(100, 200),
        Sample2 = c(150, 180)
      )
      write.csv(example_counts, file, row.names = FALSE)
    }
  )
  # Provide example template for design
  output$download_design_template <- downloadHandler(
    filename = function() { "design_template.csv" },
    content = function(file) {
      example_design <- data.frame(
        Sample = c("Sample1", "Sample2"),
        condition = c("control", "treated")
      )
      write.csv(example_design, file, row.names = FALSE)
    }
  )
  # Render design matrix
  output$design_matrix <- renderDT({
    req(data$samples)
    tryCatch({
      design_formula <- as.formula(input$design_formula)
      mm <- model.matrix(design_formula, data$samples)
      data$design_matrix <- mm
      datatable(mm, options = list(scrollX = TRUE))
    }, error = function(e) {
      showNotification(paste("Failed to build design matrix:", e$message), type = "error")
      NULL
    })
  })

  # Render diagnostic for design matrix
  output$design_diag <- renderPrint({
    req(data$design_matrix)
    diag <- tryCatch({
      list(
        Rank = qr(data$design_matrix)$rank,
        Collinearity = if (qr(data$design_matrix)$rank < ncol(data$design_matrix)) "Yes" else "No"
      )
    }, error = function(e) {
      paste("Error checking design matrix:", e$message)
    })
    str(diag)
  })

  # Show uploaded counts
  output$uploaded_counts <- renderDT({ req(data$counts); datatable(data$counts) })
  # Show uploaded design metadata
  output$uploaded_design <- renderDT({ req(data$samples); datatable(data$samples) })
}

# === Register in server ===
# mod_data_upload_design(input, output, session, data, dds_rv)

