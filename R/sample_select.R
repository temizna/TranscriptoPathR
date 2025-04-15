# === Module: mod_sample_select ===
#' Sample Selection Module
#'
#' @description This module provides functionality to select or deselect specific samples from the loaded metadata.
#' It updates the UI to reflect the available samples and allows subsetting them interactively.
#'
#' @param input Shiny input object
#' @param output Shiny output object
#' @param session Shiny session object
#' @param data ReactiveValues object containing at least `samples` (data frame of sample metadata)
#' 
#' @return None. Updates UI elements and renders selected sample table.
#' 
#' @importFrom shiny observe observeEvent updateSelectInput renderUI
#' @importFrom DT renderDT datatable
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom shinythemes shinytheme
#' @importFrom org.Mm.eg.db org.Mm.eg.db
#' @export
mod_sample_select <- function(input, output, session, data) {

  # Update sample choices when data is loaded
  observe({
    req(data$samples)
    updateSelectInput(session, "sample_select", choices = rownames(data$samples))
  })

  # Select all samples
  observeEvent(input$select_all, {
    req(data$samples)
    updateSelectInput(session, "sample_select", selected = rownames(data$samples))
  })

  # Deselect all samples
  observeEvent(input$deselect_all, {
    updateSelectInput(session, "sample_select", selected = character(0))
  })

  # Show table of selected samples
  output$sampleTable <- renderDT({
    req(data$samples, input$sample_select)
    datatable(data$samples[input$sample_select, , drop = FALSE], options = list(scrollX = TRUE))
  })
}
