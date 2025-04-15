# === Module: mod_sample_select ===

#' Sample Selection Module
#'
#' This module enables selection of samples from uploaded or demo data.
#' Users can view the design table and choose samples for downstream analyses.
#'
#' @param input Shiny input object.
#' @param output Shiny output object.
#' @param session Shiny session object.
#' @param data ReactiveValues storing `samples` (design table).
#'
#' @return None. Updates sample selection and displays filtered sample table.
#' @importFrom DT renderDT datatable
#' @importFrom shinythemes shinytheme
#' @export
mod_sample_select <- function(input, output, session, data) {

  # Update sample choices whenever design is updated
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

  # Render sample metadata table
  output$sampleTable <- renderDT({
    req(data$samples, input$sample_select)
    datatable(data$samples[input$sample_select, , drop = FALSE], options = list(scrollX = TRUE))
  })
}

