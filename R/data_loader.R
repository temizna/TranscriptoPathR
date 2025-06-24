# === Module: mod_data_upload_design ===
#' Data upload and design
#' @description Load and process uploaded count and design files into DESeq2 objects
#'
#' @param input Shiny input object.
#' @param output Shiny output object.
#' @param session Shiny session object.
#' @param loaded_data_rv ReactiveValues to store processed data (counts, samples, norm_counts).
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
#' @importFrom shiny downloadHandler renderPrint observeEvent
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom shinythemes shinytheme
#' @importFrom utils install.packages
#' @export
mod_data_upload_design <- function(input, output, session, loaded_data_rv, dds_rv) {
  
  load_data_file <- function(file_input) {
    ext <- tolower(tools::file_ext(file_input$name))
    
    df <- switch(ext,
                 "csv" = read.csv(file_input$datapath, check.names = FALSE),
                 "tsv" = read.delim(file_input$datapath, check.names = FALSE),
                 "txt" = read.delim(file_input$datapath, check.names = FALSE),
                 "xlsx" = as.data.frame(readxl::read_excel(file_input$datapath, col_names = TRUE)),
                 stop("Unsupported file format: ", ext)
    )
      row_ids <- df[[1]]
      if (any(duplicated(row_ids))) {
        warning("Duplicate row identifiers found. Applying make.unique().")
      }
      rownames(df) <- make.unique(trimws(as.character(row_ids)))
      df <- df[, -1, drop = FALSE]
    return(df)
  }
  
  load_design_file <- function(file_input) {
      ext <- tools::file_ext(file_input$name)
      if (ext == "csv") {
        read.csv(file_input$datapath, row.names = 1, check.names = FALSE)
      } else if   (ext == "tsv")  {read.delim(file_input$datapath,row.names = 1, check.names = FALSE)
        } else if (ext == "txt") { read.delim(file_input$datapath, row.names = 1,check.names = FALSE)
          } else if (ext == "xlsx") {
        readxl::read_excel(file_input$datapath, col_names = TRUE) |> tibble::column_to_rownames(1)
      } else {
        stop("Unsupported file format: ", ext)
      }
    }
    
  load_uploaded_data <- function(counts_file, design_file, design_formula) {
    counts_df <- load_data_file(counts_file)
    design_df <- load_design_file(design_file)
    
    valid_sample_idx <- rownames(design_df) %in% colnames(counts_df)
    design_df <- design_df[valid_sample_idx, , drop = FALSE]
    counts_df <- counts_df[, rownames(design_df), drop = FALSE]
    
    excluded_design <- setdiff(rownames(design_df), colnames(counts_df))
    excluded_counts <- setdiff(colnames(counts_df), rownames(design_df))
    
    if (length(excluded_design) > 0)
      warning("Some samples in design were not found in counts: ", paste(excluded_design, collapse = ", "))
    if (length(excluded_counts) > 0)
      warning("Some samples in counts were not found in design: ", paste(excluded_counts, collapse = ", "))
    
    rn <- rownames(design_df)
    design_df <- as.data.frame(lapply(design_df, function(col) {
      if (is.character(col) || is.logical(col)) factor(col) else col
    }))
    rownames(design_df) <- rn
    
    counts_df <- counts_df[rowSums(counts_df) >= 10, , drop = FALSE]
    counts_df <- counts_df[apply(counts_df, 1, var) > 0.1, , drop = FALSE]
    
    if (nrow(counts_df) < 10) stop("Too few genes left after filtering. Check your input files.")
    
    # Sanitize and collapse rows with same gene ID per sample
    clean_ids <- sanitize_ensembl_ids(rownames(counts_df))
    counts_df$gene_id <- clean_ids
    
    counts_df <- as.data.frame(counts_df)
    counts_df <- aggregate(. ~ gene_id, data = counts_df, FUN = sum)
    rownames(counts_df) <- counts_df$gene_id
    counts_df$gene_id <- NULL
    
    formula <- as.formula(design_formula)
    required_vars <- all.vars(formula)[-1]
    if (!all(required_vars %in% colnames(design_df))) {
      stop("Design formula references unknown columns in metadata: ", 
           paste(setdiff(required_vars, colnames(design_df)), collapse = ", "))
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
      
      dds <- DESeq2::DESeqDataSetFromMatrix(countData = demo_counts, colData = demo_design, design = ~ condition)
      dds <- DESeq2::estimateSizeFactors(dds)
      
      loaded_data_rv$counts <- demo_counts
      loaded_data_rv$samples <- demo_design
      loaded_data_rv$norm_counts <- counts(dds, normalized = TRUE)
      loaded_data_rv$species <- "Homo sapiens"
      
      dds_rv(dds)
    }
  })
  
  # Observe event for loading user-uploaded data
  observeEvent(input$load_data, {
    req(input$counts_file, input$design_file, input$design_formula)
    tryCatch({
      loaded <- load_uploaded_data(input$counts_file, input$design_file, input$design_formula)
      
      loaded_data_rv$counts <- loaded$counts
      loaded_data_rv$samples <- loaded$samples
      loaded_data_rv$norm_counts <- loaded$norm_counts
      loaded_data_rv$species <- input$species
      
      dds_rv(loaded$dds)
      showNotification("Data successfully loaded!", type = "message")
    }, error = function(e) {
      showNotification(paste("Error loading files:", e$message), type = "error")
    })
  })
  
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
  
  output$design_matrix <- renderDT({
    req(loaded_data_rv$samples)
    tryCatch({
      design_formula <- as.formula(input$design_formula)
      mm <- model.matrix(design_formula, loaded_data_rv$samples)
      loaded_data_rv$design_matrix <- mm
      datatable(mm, options = list(scrollX = TRUE))
    }, error = function(e) {
      showNotification(paste("Failed to build design matrix:", e$message), type = "error")
      NULL
    })
  })
  
  output$design_diag <- renderPrint({
    req(loaded_data_rv$design_matrix)
    diag <- tryCatch({
      list(
        Rank = qr(loaded_data_rv$design_matrix)$rank,
        Collinearity = if (qr(loaded_data_rv$design_matrix)$rank < ncol(loaded_data_rv$design_matrix)) "Yes" else "No"
      )
    }, error = function(e) {
      paste("Error checking design matrix:", e$message)
    })
    str(diag)
  })
  
  output$uploaded_counts <- renderDT({
    req(loaded_data_rv$counts)
    datatable(loaded_data_rv$counts)
  })
  
  output$uploaded_design <- renderDT({
    req(loaded_data_rv$samples)
    datatable(loaded_data_rv$samples)
  })
}
