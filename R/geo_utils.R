# === Module: mod_geo_loader ===

#' GEO Loader Module
#'
#' @description Attempts to download raw counts and metadata from a GEO accession.
#' Checks if valid GEO ID, counts, and metadata are available.
#'
#' @param input Shiny input object
#' @param output Shiny output object
#' @param session Shiny session object
#' @param data ReactiveValues storing counts, samples, norm_counts, species
#' @import GEOquery Biobase readr
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom shinythemes shinytheme
#' @export
mod_geo_loader <- function(input, output, session, data) {
  observeEvent(input$load_geo, {
    req(input$geo_accession)

    geo_id <- trimws(input$geo_accession)
    showNotification(paste("Loading data from", geo_id), type = "message")

    tryCatch({
      gse <- GEOquery::getGEO(geo_id, GSEMatrix = TRUE)
      if (length(gse) > 1) gse <- gse[[1]]

      # Extract expression matrix
      exprs_data <- Biobase::exprs(gse)
      pdata <- Biobase::pData(gse)

      # Ensure expression values are counts (integers)
      if (!all(exprs_data == floor(exprs_data))) {
        showNotification("GEO dataset does not contain raw count data.", type = "error")
        return(NULL)
      }

      if (ncol(exprs_data) < 2 || nrow(exprs_data) < 10) {
        showNotification("Expression data too small or invalid format.", type = "error")
        return(NULL)
      }

      # Check metadata
      if (nrow(pdata) != ncol(exprs_data)) {
        showNotification("Mismatch between expression and metadata samples.", type = "error")
        return(NULL)
      }

      # Match sample names
      colnames(exprs_data) <- rownames(pdata)

      # Set species (best guess from title or metadata)
      title <- gse@header$organism_ch1[1]
      species <- ifelse(grepl("mus musculus", title, ignore.case = TRUE), "Mus musculus", "Homo sapiens")

      # Store
      data$counts <- exprs_data
      data$samples <- pdata
      data$species <- species

      showNotification("GEO data loaded successfully.", type = "message")

    }, error = function(e) {
      msg <- e$message
      if (grepl("ftp|html|cannot open", msg, ignore.case = TRUE)) {
        showNotification("GEO data could not be retrieved. Check network or GEO ID.", type = "error")
      } else {
        showNotification(paste("GEO loading failed:", msg), type = "error")
      }
    })
  })
}

