# === Module: mod_sample_select_dynamic ===

#' Sample Selection Module (Dynamic)
#'
#' This module dynamically creates sample filtering UI elements based on the metadata 
#' and allows filtering of samples for downstream analyses.
#'
#' @param input Shiny input object.
#' @param output Shiny output object.
#' @param session Shiny session object.
#' @param dds_rv reactiveVal containing the full DESeq2 object.
#' @param loaded_data_rv reactiveValues containing loaded raw data.
#' @param filtered_data_rv reactiveVal for filtered dataset.
#' @param filtered_dds_rv reactiveVal for filtered DESeq2 object.
#' @importFrom DT renderDT datatable
#' @importFrom shiny NS selectInput tagList
#' @export
mod_sample_select_dynamic <- function(input, output, session, dds_rv, loaded_data_rv, filtered_data_rv, filtered_dds_rv) {
  filter_input_ids <- reactiveVal(NULL)
  
  # Dynamically generate UI filters based on loaded sample metadata
  output$dynamic_filters <- renderUI({
    req(loaded_data_rv$samples)
    samples <- loaded_data_rv$samples
    
    ns <- session$ns  # get the module namespace
    
    filter_inputs <- lapply(names(samples), function(col_name) {
      selectInput(
        inputId = ns(paste0("filter_", col_name)),  # <- namespaced here
        label = paste("Filter by", col_name),
        choices = unique(samples[[col_name]]),
        selected = NULL,
        multiple = TRUE
      )
    })
    
    do.call(tagList, filter_inputs)
  })
  
  # Apply filtering when user clicks "Apply Filters"
  observeEvent(input$run_filter, {
    req(loaded_data_rv$samples, dds_rv())
    
    full_data <- loaded_data_rv
    samples <- full_data$samples
    filtered_samples <- samples
    filter_ids <- grep("^filter_", names(input), value = TRUE)
    
    for (filter_id in filter_ids) {
      col_name <- sub("^filter_", "", filter_id)
      selected_vals <- input[[filter_id]]
      if (!is.null(selected_vals) && length(selected_vals) > 0 && col_name %in% colnames(filtered_samples)) {
        filtered_samples <- filtered_samples[filtered_samples[[col_name]] %in% selected_vals, , drop = FALSE]
      }
    }
    
    if (nrow(filtered_samples) == 0) {
      showNotification("No samples match the selected filters.", type = "error")
      return()
    }
    
    sample_names <- rownames(filtered_samples)
   # print("Filtered sample names:")
  #  print(sample_names)
    #print("Filtered samples:")
    #print(head(filtered_samples))
    
    # Replace all fields at once
    filtered_data_rv$counts <- full_data$counts[, sample_names, drop = FALSE]
    filtered_data_rv$norm_counts <- full_data$norm_counts[, sample_names, drop = FALSE]
    filtered_data_rv$samples <- filtered_samples
    filtered_data_rv$species <- full_data$species
    
    # Replace DESeq2 object
    filtered_dds_rv(dds_rv()[, sample_names])
  })
  
  # Deselect all filters and reset to full dataset
  observeEvent(input$deselect_all, {
    req(loaded_data_rv, dds_rv())
    filtered_data_rv$counts <- loaded_data_rv$counts
    filtered_data_rv$norm_counts <- loaded_data_rv$norm_counts
    filtered_data_rv$samples <- loaded_data_rv$samples
    filtered_data_rv$species<- loaded_data_rv$species
    filtered_dds_rv(dds_rv())
  })
  
  # Select all samples again
  observeEvent(input$select_all, {
    req(loaded_data_rv, dds_rv())
    filtered_data_rv$counts <- loaded_data_rv$counts
    filtered_data_rv$norm_counts <- loaded_data_rv$norm_counts
    filtered_data_rv$samples <- loaded_data_rv$samples
    filtered_data_rv$species<- loaded_data_rv$species
    filtered_dds_rv(dds_rv())
    #print(filtered_data_rv$species)
  })
  
  # Render the filtered sample table
  output$sampleTable <- renderDT({
    req(filtered_data_rv$samples)
    datatable(filtered_data_rv$samples)
  })
}
