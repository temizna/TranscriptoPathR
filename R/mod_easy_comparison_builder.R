#' Easy Comparison Builder Server
#'
#' Uses filtered_data_rv$samples to define Test vs Reference comparisons
#' (simple, pooled, paired) and optional covariates to adjust.
#'
#' @param id module id
#' @param filtered_data_rv reactiveValues; must contain $samples (data.frame with rownames as sample IDs)
#' @return list of reactives (see code)
#' @importFrom shiny moduleServer reactive reactiveVal req observe observeEvent updateSelectInput
#' @importFrom shiny renderUI renderTable showNotification validate need
#' @export
mod_easy_compare_server <- function(id, filtered_data_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # ---- helpers ------------------------------------------------------------
    safe_tag <- function(x) gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
    guess_subject_col <- function(cols) {
      pat <- "(^|_)(subject|patient|donor).*|.*(_|^)(id)$"
      ix <- grep(pat, tolower(cols))
      if (length(ix)) cols[ix[1]] else cols[1]
    }
    
    # Wrap samples in a reactive so changes always invalidate downstream
    smp_df <- shiny::reactive({
      shiny::req(!is.null(filtered_data_rv$samples))
      x <- filtered_data_rv$samples
      # coerce tibbles/data.tables to plain df; keep character/factor levels as-is
      x <- as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
      shiny::validate(shiny::need(ncol(x) > 0, "No sample metadata columns found."))
      x
    })
    
    # ---- 1) Populate grouping column ---------------------------------------
    shiny::observeEvent(smp_df(), {
      cols <- colnames(smp_df())
      shiny::updateSelectInput(
        session, "group_var",
        choices  = cols,
        selected = if (length(cols)) cols[1] else character(0)
      )
    }, ignoreInit = FALSE)
    
    # ---- 2) Mode-specific controls -----------------------------------------
    output$mode_controls <- shiny::renderUI({
      smp <- smp_df()
      shiny::req(input$group_var)
      shiny::validate(shiny::need(input$group_var %in% colnames(smp),
                                  "Select a valid grouping column."))
      lv <- unique(as.character(smp[[input$group_var]]))
      lv <- lv[!is.na(lv)]
      
      if (identical(input$mode, "simple")) {
        shiny::tagList(
          shiny::selectInput(ns("level_test"), "Test group", choices = lv,
                             selected = if (length(lv)) lv[1] else NULL),
          shiny::selectInput(ns("level_ref"), "Reference group", choices = lv,
                             selected = if (length(lv) >= 2) lv[2] else if (length(lv)) lv[1] else NULL)
        )
      } else if (identical(input$mode, "pooled")) {
        shiny::tagList(
          shiny::selectInput(ns("levels_test"), "Test group levels (one or more)",
                             choices = lv, multiple = TRUE),
          shiny::selectInput(ns("levels_ref"), "Reference group levels (one or more)",
                             choices = lv, multiple = TRUE),
          shiny::helpText("Optional: restrict to only the selected levels below."),
          shiny::selectInput(ns("keep_levels_only"), "Keep only selected levels",
                             choices = c("Yes" = "yes", "No" = "no"), selected = "yes")
        )
      } else {
        cols <- colnames(smp)
        subj_default <- guess_subject_col(cols)
        shiny::tagList(
          shiny::selectInput(ns("subject_col"), "Subject ID column",
                             choices = cols, selected = subj_default),
          shiny::selectInput(ns("level_test"), "Test group", choices = lv,
                             selected = if (length(lv)) lv[1] else NULL),
          shiny::selectInput(ns("level_ref"), "Reference group", choices = lv,
                             selected = if (length(lv) >= 2) lv[2] else if (length(lv)) lv[1] else NULL)
        )
      }
    })
    
    # ---- 3) Covariates ------------------------------------------------------
    output$covariate_ui <- shiny::renderUI({
      smp <- smp_df()
      shiny::req(input$group_var)
      cols <- setdiff(colnames(smp), input$group_var)
      shiny::selectInput(ns("covariates"), "Adjust for (optional)",
                         choices = cols, multiple = TRUE)
    })
    
    # ---- 4) Previews --------------------------------------------------------
    output$preview_counts_ui <- shiny::renderUI({
      shiny::req(input$group_var)
      shiny::tagList(
        shiny::helpText("Sample counts by level in the selected grouping column:"),
        shiny::tableOutput(ns("counts_by_level"))
      )
    })
    
    output$counts_by_level <- shiny::renderTable({
      smp <- smp_df()
      shiny::req(input$group_var)
      shiny::validate(shiny::need(input$group_var %in% colnames(smp), ""))
      tab <- as.data.frame(table(level = as.character(smp[[input$group_var]])))
      colnames(tab) <- c("level", "n")
      tab[order(tab$n, decreasing = TRUE), , drop = FALSE]
    }, striped = TRUE, bordered = TRUE, spacing = "xs")
    
    # ---- 5) State for exports ----------------------------------------------
    cmp_group_name <- shiny::reactive("cmp_group")
    test_label     <- shiny::reactiveVal(NULL)
    ref_label      <- shiny::reactiveVal(NULL)
    label_out      <- shiny::reactiveVal(NULL)
    tag_out        <- shiny::reactiveVal(NULL)
    included_out   <- shiny::reactiveVal(character())
    cmp_factor_out <- shiny::reactiveVal(NULL)
    
    # ---- 6) Apply -----------------------------------------------------------
    shiny::observeEvent(input$apply, {
      smp <- smp_df()
      shiny::req(input$group_var)
      gv <- input$group_var
      shiny::validate(shiny::need(gv %in% colnames(smp), "Invalid grouping column."))
      
      grp_chr <- as.character(smp[[gv]])
      
      if (identical(input$mode, "simple")) {
        shiny::req(input$level_test, input$level_ref)
        if (identical(input$level_test, input$level_ref)) {
          shiny::showNotification("Test and Reference must be different.", type = "error"); return()
        }
        role <- ifelse(grp_chr == input$level_test, "Test",
                       ifelse(grp_chr == input$level_ref, "Reference", NA))
        keep <- which(!is.na(role))
        if (!length(keep)) { shiny::showNotification("No samples match the selected levels.", type = "error"); return() }
        role    <- factor(role[keep], levels = c("Reference", "Test"))
        smp_sub <- smp[keep, , drop = FALSE]
        test_label(input$level_test); ref_label(input$level_ref)
        
      } else if (identical(input$mode, "pooled")) {
        shiny::req(input$levels_test, input$levels_ref)
        if (!length(input$levels_test) || !length(input$levels_ref)) {
          shiny::showNotification("Select at least one level for both Test and Reference.", type = "error"); return()
        }
        if (length(intersect(input$levels_test, input$levels_ref))) {
          shiny::showNotification("Test and Reference pools cannot overlap.", type = "error"); return()
        }
        role_all <- ifelse(grp_chr %in% input$levels_test, "Test",
                           ifelse(grp_chr %in% input$levels_ref, "Reference", NA))
        keep <- if (identical(input$keep_levels_only, "yes")) which(!is.na(role_all)) else seq_len(nrow(smp))
        if (!length(keep)) { shiny::showNotification("No samples selected.", type = "error"); return() }
        role    <- factor(role_all[keep], levels = c("Reference", "Test"))
        smp_sub <- smp[keep, , drop = FALSE]
        test_label(paste(input$levels_test, collapse = "+"))
        ref_label(paste(input$levels_ref, collapse = "+"))
        
      } else { # paired
        shiny::req(input$subject_col, input$level_test, input$level_ref)
        if (identical(input$level_test, input$level_ref)) {
          shiny::showNotification("Test and Reference must be different.", type = "error"); return()
        }
        mask  <- grp_chr %in% c(input$level_test, input$level_ref)
        smp2  <- smp[mask, , drop = FALSE]
        if (!nrow(smp2)) { shiny::showNotification("No samples for the selected levels.", type = "error"); return() }
        role2 <- ifelse(as.character(smp2[[gv]]) == input$level_test, "Test", "Reference")
        subj  <- as.character(smp2[[input$subject_col]])
        tab   <- table(subj, role2)
        both  <- rownames(tab)[rowSums(tab > 0) == 2]
        keep  <- which(subj %in% both)
        if (!length(keep)) { shiny::showNotification("No subjects with both Test and Reference.", type = "error"); return() }
        smp_sub <- smp2[keep, , drop = FALSE]
        role    <- factor(role2[keep], levels = c("Reference", "Test"))
        test_label(input$level_test); ref_label(input$level_ref)
      }
      
      included_out(rownames(smp_sub))
      cmp_factor_out(role)
      label_out(paste0(test_label(), " vs ", ref_label()))
      tag_out(paste0(safe_tag(make.names(test_label())), "_vs_", safe_tag(make.names(ref_label()))))
      
      output$preview_samples <- shiny::renderTable({
        data.frame(
          sample = rownames(smp_sub),
          role   = as.character(role),
          smp_sub,
          row.names  = NULL,
          check.names = FALSE
        )
      }, striped = TRUE, bordered = TRUE, spacing = "xs")
      
      output$summary <- shiny::renderTable({
        data.frame(
          grouping_column = gv,
          test_label      = test_label(),
          reference_label = ref_label(),
          group_var       = cmp_group_name(),
          test_level      = "Test",
          reference_level = "Reference",
          adjusted_for    = if (length(input$covariates)) paste(input$covariates, collapse = ", ") else "",
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
      }, striped = TRUE, bordered = TRUE, spacing = "xs")
      
      shiny::showNotification("Comparison ready. Downstream modules can now use it.", type = "message")
    })
    
    # ---- Public API ---------------------------------------------------------
    list(
      group_var        = shiny::reactive(cmp_group_name()),
      ref_level        = shiny::reactive("Reference"),
      test_level       = shiny::reactive("Test"),
      label            = shiny::reactive(label_out()),
      tag              = shiny::reactive(tag_out()),
      included_samples = shiny::reactive(included_out()),
      covariates       = shiny::reactive(input$covariates),
      cmp_factor       = shiny::reactive(cmp_factor_out())
    )
  })
}
