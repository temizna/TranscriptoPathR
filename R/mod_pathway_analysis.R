#' Pathway Analysis Server (namespaced, full logic; noisy diagnostics)
#'
#' @param id module id (must match mod_pathway_tab_ui)
#' @param filtered_data_rv reactiveValues with $counts, $norm_counts, $samples, $species
#' @param res_reactive reactiveVal holding DESeq2 results (data.frame)
#' @param geneList_rv reactiveVal to store named log2FC vector (ENTREZID -> log2FC)
#' @param d1_merged_rv reactiveVal to store merged df (gene, log2FC, padj, ENTREZID)
#' @param pathway_result_rv reactiveVal to store the enrichment result object
#' @param de_sel optional list returned by mod_de_server(), used for filenames if cmp is NULL
#' @param cmp optional comparison bridge from make_cmp_bridge(), preferred for filenames
#' @export
mod_pathway_server <- function(
    id,
    filtered_data_rv,
    res_reactive,
    geneList_rv,
    d1_merged_rv,
    pathway_result_rv,
    de_sel = NULL,
    cmp = NULL
) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    .get_tag <- function() {
      tryCatch({
        if (!is.null(cmp)) {
          t <- cmp$tag()
          if (!is.null(t) && nzchar(t)) return(t)
        }
        .contrast_tag_from(de_sel)
      }, error = function(e) .contrast_tag_from(de_sel))
    }
    
    # ---- Custom DE file loader ----------------------------------------------
    custom_de_data <- reactive({
      req(input$pathway_use_custom_file)
      file <- input$pathway_custom_file
      req(file)
      
      ext <- tools::file_ext(file$datapath)
      df <- if (ext %in% c("csv")) {
        read.csv(file$datapath, stringsAsFactors = FALSE)
      } else {
        read.delim(file$datapath, stringsAsFactors = FALSE)
      }
      
      # Validate required columns
      required_cols <- c("gene", "log2FoldChange", "padj")
      if (!all(required_cols %in% colnames(df))) {
        shiny::showNotification("Uploaded file must contain columns: 'gene', 'log2FoldChange', 'padj'.", type = "error")
        return(NULL)
      }
      
      df <- df[!is.na(df$log2FoldChange) & !is.na(df$padj), , drop = FALSE]
      df$gene <- make.unique(as.character(df$gene))
      rownames(df) <- df$gene
      df
    })
    
    # ---- Run enrichment ------------------------------------------------------
    shiny::observeEvent(input$run_pathway, {
      shiny::req(filtered_data_rv$species)
      shiny::withProgress(message = "Pathway analysis", value = 0, {
        incProgress(0.05, detail = "Preparing inputs")
        
        species <- filtered_data_rv$species
        orgdb   <- .get_orgdb(species)
        if (is.null(orgdb)) {
          shiny::showNotification("OrgDb not available for the selected species.", type = "error")
          return(NULL)
        }
        
        # Decide DE source
        res <- if (isTRUE(input$pathway_use_custom_file)) custom_de_data() else res_reactive()
        if (is.null(res) || !nrow(res)) {
          shiny::showNotification("No DE results available. Upload file or run DE analysis first.", type = "warning")
          return(NULL)
        }
        
        # Ensure needed columns exist
        need_cols <- c("log2FoldChange", "padj")
        if (!all(need_cols %in% colnames(res))) {
          shiny::showNotification("DE results missing required columns (log2FoldChange/padj).", type = "error")
          return(NULL)
        }
        
        # ID prep
        res$gene <- if ("gene" %in% colnames(res)) as.character(res$gene) else rownames(res)
        res$gene_base <- sub("\\.\\d+$", "", res$gene)
        
        # Map SYMBOL or ENSEMBL -> ENTREZ
        from_type <- if (.is_symbol(res$gene_base)) "SYMBOL" else "ENSEMBL"
        incProgress(0.15, detail = paste0("Mapping IDs (", from_type, " -> ENTREZ)"))
        d1_ids <- suppressMessages({
          clusterProfiler::bitr(res$gene_base, fromType = from_type, toType = "ENTREZID", OrgDb = orgdb)
        })
        if (is.null(d1_ids) || nrow(d1_ids) == 0) {
          shiny::showNotification("No genes could be mapped to ENTREZ IDs.", type = "error")
          return(NULL)
        }
        
        d1_merged <- merge(res, d1_ids, by.x = "gene_base", by.y = 1)
        d1_merged <- d1_merged[!duplicated(d1_merged$ENTREZID), , drop = FALSE]
        d1_merged <- d1_merged[!is.na(d1_merged$log2FoldChange) & !is.na(d1_merged$padj), , drop = FALSE]
        d1_merged_rv(d1_merged)
        
        n_in   <- nrow(res)
        n_map  <- nrow(d1_merged)
        shiny::showNotification(paste0("Mapped ", n_map, " / ", n_in, " DE rows to ENTREZ."), type = "message", duration = 4)
        
        # Filter by direction
        incProgress(0.25, detail = "Filtering by direction and cutoffs")
        dir <- input$pathway_direction
        keep <- switch(
          dir,
          "Up"   = d1_merged$log2FoldChange >=  input$lfc_threshold & d1_merged$padj <= input$padj_threshold,
          "Down" = d1_merged$log2FoldChange <= -input$lfc_threshold & d1_merged$padj <= input$padj_threshold,
          abs(d1_merged$log2FoldChange) >= input$lfc_threshold & d1_merged$padj <= input$padj_threshold
        )
        gene_vector <- d1_merged[keep & !is.na(d1_merged$ENTREZID), , drop = FALSE]
        n_keep <- nrow(gene_vector)
        shiny::showNotification(paste0("Genes after filtering: ", n_keep), type = "message", duration = 4)
        
        if (!n_keep) {
          shiny::showNotification("No genes pass the selected cutoffs.", type = "warning")
          return(NULL)
        }
        
        # Ranked vector (named by ENTREZ)
        geneList <- gene_vector$log2FoldChange
        names(geneList) <- gene_vector$ENTREZID
        geneList <- geneList[!duplicated(names(geneList))]
        geneList <- sort(geneList, decreasing = TRUE)
        max_genes <- if (!is.null(input$max_genes)) input$max_genes else 1000
        if (length(geneList) > max_genes) geneList <- head(geneList, max_genes)
        
        if (length(geneList) < 10) {
          shiny::showNotification("Too few mapped genes (<10) for enrichment.", type = "warning")
          return(NULL)
        }
        geneList_rv(geneList)
        selected_genes <- names(geneList)
        
        qcut <- input$pathway.qval
        if (is.null(qcut) || is.na(qcut)) qcut <- 0.2
        
        incProgress(0.65, detail = paste0("Running enrichment: ", input$pathway_db))
        pathway_result <- NULL
        
        # ---- GO Enrichment ----
        if (identical(input$pathway_db, "GO")) {
          pathway_result <- tryCatch({
            clusterProfiler::enrichGO(
              gene = selected_genes,
              OrgDb = orgdb,
              keyType = "ENTREZID",
              ont = "BP",
              pAdjustMethod = "BH",
              pvalueCutoff = input$padj_threshold,
              qvalueCutoff = qcut,
              readable = TRUE
            )
          }, error = function(e) {
            shiny::showNotification(paste("GO enrichment failed:", e$message), type = "error")
            NULL
          })
          
          # ---- KEGG Enrichment ----
        } else if (identical(input$pathway_db, "KEGG")) {
          kegg_code <- .get_kegg_code(species)
          if (is.null(kegg_code)) {
            shiny::showNotification("KEGG not supported for this species.", type = "error")
            return(NULL)
          }
          kegg_raw <- tryCatch({
            clusterProfiler::enrichKEGG(
              gene = selected_genes,
              organism = kegg_code,
              pvalueCutoff = input$padj_threshold,
              qvalueCutoff = qcut
            )
          }, error = function(e) {
            shiny::showNotification(paste("KEGG enrichment failed:", e$message), type = "error")
            NULL
          })
          if (!is.null(kegg_raw) && nrow(as.data.frame(kegg_raw)) > 0) {
            pathway_result <- tryCatch({
              clusterProfiler::setReadable(kegg_raw, OrgDb = orgdb, keyType = "ENTREZID")
            }, error = function(e) {
              shiny::showNotification(paste("setReadable failed:", e$message), type = "error")
              kegg_raw
            })
          }
          
          # ---- Reactome Enrichment ----
        } else if (identical(input$pathway_db, "Reactome")) {
          reactome_code <- .get_reactome_code(species)
          if (is.null(reactome_code)) {
            shiny::showNotification("Reactome supported only for human, mouse, and rat.", type = "error")
            return(NULL)
          }
          pathway_result <- tryCatch({
            ReactomePA::enrichPathway(
              gene = selected_genes,
              organism = reactome_code,
              pvalueCutoff = input$padj_threshold,
              qvalueCutoff = qcut,
              readable = TRUE
            )
          }, error = function(e) {
            shiny::showNotification(paste("Reactome enrichment failed:", e$message), type = "error")
            NULL
          })
        } else {
          shiny::showNotification("Unknown pathway database selection.", type = "error")
          return(NULL)
        }
        
        incProgress(0.85, detail = "Post-processing results")
        if (!is.null(pathway_result) && nrow(as.data.frame(pathway_result)) > 0) {
          if (!"geneSymbol" %in% colnames(as.data.frame(pathway_result))) {
            pathway_result <- tryCatch({
              clusterProfiler::setReadable(pathway_result, OrgDb = orgdb, keyType = "ENTREZID")
            }, error = function(e) pathway_result)
          }
          pathway_result <- tryCatch({
            enrichplot::pairwise_termsim(pathway_result)
          }, error = function(e) pathway_result)
          
          shiny::showNotification(
            paste0("Enrichment done. Terms: ", nrow(as.data.frame(pathway_result))),
            type = "message", duration = 4
          )
        } else {
          shiny::showNotification("No enriched pathways found with current settings.", type = "warning")
        }
        
        pathway_result_rv(pathway_result)
        incProgress(1)
      })
    })
    
    # ---- Plots & Downloads ---------------------------------------------------
    output$dotPlot <- shiny::renderPlot({
      shiny::req(pathway_result_rv())
      pr <- pathway_result_rv()
      if (is.null(pr) || nrow(as.data.frame(pr)) == 0) {
        shiny::showNotification("No enrichment terms available for dotplot.", type = "warning")
        return(NULL)
      }
      enrichplot::dotplot(pr) +
        ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6, face = "bold"))
    })
    
    output$download_dot_plot <- shiny::downloadHandler(
      filename = function() paste0("Pathway_", input$pathway_db, "_", input$pathway_direction, "_", .get_tag(), "_dot_plot.pdf"),
      content  = function(file) {
        shiny::req(pathway_result_rv())
        grDevices::pdf(file)
        print(enrichplot::dotplot(pathway_result_rv()) +
                ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6, face = "bold")))
        grDevices::dev.off()
      }
    )
    
    output$emapPlot <- shiny::renderPlot({
      shiny::req(pathway_result_rv())
      pr <- pathway_result_rv()
      df <- as.data.frame(pr)
      if (is.null(df) || nrow(df) < 2 || !"ID" %in% colnames(df) || !"geneID" %in% colnames(df)) {
        shiny::showNotification("Insufficient enrichment terms or missing columns for emapplot.", type = "warning")
        return(NULL)
      }
      if (is.null(pr@termsim) || all(is.na(pr@termsim)) || all(pr@termsim == 0)) {
        shiny::showNotification("No valid term similarity matrix for emapplot.", type = "warning")
        return(NULL)
      }
      enrichplot::emapplot(pr, showCategory = 10)
    })
    
    output$download_emap_plot <- shiny::downloadHandler(
      filename = function() paste0("Pathway_", input$pathway_db, "_", input$pathway_direction, "_", .get_tag(), "_emap_plot.pdf"),
      content  = function(file) {
        shiny::req(pathway_result_rv())
        grDevices::pdf(file)
        print(enrichplot::emapplot(pathway_result_rv(), showCategory = 10))
        grDevices::dev.off()
      }
    )
    
    output$cnetPlot <- shiny::renderPlot({
      shiny::req(pathway_result_rv())
      pr <- pathway_result_rv()
      if (is.null(pr) || nrow(as.data.frame(pr)) == 0) {
        shiny::showNotification("No enrichment terms available for cnetplot.", type = "warning")
        return(NULL)
      }
      enrichplot::cnetplot(pr, showCategory = 10)
    })
    
    output$download_cnet_plot <- shiny::downloadHandler(
      filename = function() paste0("Pathway_", input$pathway_db, "_", input$pathway_direction, "_", .get_tag(), "_cnet_plot.pdf"),
      content  = function(file) {
        shiny::req(pathway_result_rv())
        grDevices::pdf(file)
        print(enrichplot::cnetplot(pathway_result_rv(), showCategory = 10))
        grDevices::dev.off()
      }
    )
    
    output$circularPlot <- shiny::renderPlot({
      shiny::req(pathway_result_rv(), geneList_rv())
      pr <- pathway_result_rv()
      if (is.null(as.data.frame(pr)) || nrow(as.data.frame(pr)) == 0) {
        shiny::showNotification("No enriched terms to show circular plot.", type = "warning")
        return(NULL)
      }
      enrichplot::cnetplot(
        pr,
        layout = input$circular_layout,
        foldChange = geneList_rv(),
        showCategory = 5,
        circular = TRUE,
        colorEdge = TRUE
      )
    })
    
    output$download_circular_plot <- shiny::downloadHandler(
      filename = function() paste0("Pathway_", input$pathway_db, "_", input$pathway_direction, "_", input$circular_layout, "_", .get_tag(), "_circular_plot.pdf"),
      content  = function(file) {
        shiny::req(pathway_result_rv(), geneList_rv())
        grDevices::pdf(file)
        print(enrichplot::cnetplot(
          pathway_result_rv(),
          layout = input$circular_layout,
          foldChange = geneList_rv(),
          showCategory = 5,
          circular = TRUE,
          colorEdge = TRUE
        ))
        grDevices::dev.off()
      }
    )
    
    output$pathwayTable <- DT::renderDT({
      shiny::req(pathway_result_rv())
      DT::datatable(as.data.frame(pathway_result_rv()), options = list(scrollX = TRUE, pageLength = 25))
    })
    
    output$download_pathway_table <- shiny::downloadHandler(
      filename = function() paste0("Pathway_", input$pathway_db, "_", input$pathway_direction, "_", .get_tag(), "_results.csv"),
      content  = function(file) {
        shiny::req(pathway_result_rv())
        utils::write.csv(as.data.frame(pathway_result_rv()), file, row.names = FALSE)
      }
    )
  })
}


# ---- helpers (internal) ------------------------------------------------------

.is_symbol <- function(ids) {
  if (exists("is_symbol", mode = "function")) return(get("is_symbol")(ids))
  !grepl("^ENS[A-Z]*\\d+", ids)
}

.get_orgdb <- function(species) {
  if (exists("get_orgdb", mode = "function")) return(get("get_orgdb")(species))
  sp <- as.character(species)
  if (sp %in% c("Homo sapiens", "human")) {
    if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) return(org.Hs.eg.db::org.Hs.eg.db)
  } else if (sp %in% c("Mus musculus", "mouse")) {
    if (requireNamespace("org.Mm.eg.db", quietly = TRUE)) return(org.Mm.eg.db::org.Mm.eg.db)
  } else if (sp %in% c("Rattus norvegicus", "rat")) {
    if (requireNamespace("org.Rn.eg.db", quietly = TRUE)) return(org.Rn.eg.db::org.Rn.eg.db)
  }
  NULL
}

.get_kegg_code <- function(species) {
  if (exists("get_kegg_code", mode = "function")) return(get("get_kegg_code")(species))
  sp <- as.character(species)
  if (sp %in% c("Homo sapiens", "human")) return("hsa")
  if (sp %in% c("Mus musculus", "mouse")) return("mmu")
  if (sp %in% c("Rattus norvegicus", "rat")) return("rno")
  if (sp %in% c("Danio rerio")) return("dre")
  if (sp %in% c("Bos taurus")) return("bta")
  NULL
}

.get_reactome_code <- function(species) {
  if (exists("get_reactome_code", mode = "function")) return(get("get_reactome_code")(species))
  sp <- as.character(species)
  if (sp %in% c("Homo sapiens", "human")) return("human")
  if (sp %in% c("Mus musculus", "mouse")) return("mouse")
  if (sp %in% c("Rattus norvegicus", "rat")) return("rat")
  NULL
}
