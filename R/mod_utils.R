#' Convert Ensembl IDs to Gene Symbols
#'
#' @param gene_ids Character vector of Ensembl IDs
#' @param species Character string
#' @return Named character vector
#' @export
convert_ensembl_to_symbol <- function(gene_ids, species = "Homo sapiens") {
  org_db <- get_orgdb(species)
  id_type <- "ENSEMBL"
  res <- tryCatch({
    AnnotationDbi::select(org_db, keys = gene_ids, columns = "SYMBOL", keytype = id_type)
  }, error = function(e) {
    warning("Ensembl to Symbol conversion failed:", e$message)
    return(NULL)
  })
  out <- setNames(res$SYMBOL, res$ENSEMBL)
  return(out[!duplicated(names(out))])
}

#' Convert Gene Symbols to Ensembl IDs
#'
#' @param gene_symbols Character vector of gene symbols
#' @param species Character string
#' @return Named character vector
#' @export
convert_symbol_to_ensembl <- function(gene_symbols, species = "Homo sapiens") {
  org_db <- get_orgdb(species)
  id_type <- "SYMBOL"
  res <- tryCatch({
    AnnotationDbi::select(org_db, keys = gene_symbols, columns = "ENSEMBL", keytype = id_type)
  }, error = function(e) {
    warning("Symbol to Ensembl conversion failed: ", e$message)
    return(NULL)
  })
  out <- setNames(res$ENSEMBL, res$SYMBOL)
  return(out[!duplicated(names(out))])
}

utils::globalVariables(c(
  "log2FoldChange", "padj", "label", "baseMean", "color",
  "log2FoldChange_x", "log2FoldChange_y","Gene",
  "PC1", "PC2", "Sample", "Group", "Mean", "Variance"
))

#' Check if IDs are Ensembl
#' @param ids Character vector of gene symbols
#' @export
is_ensembl_id <- function(ids) {
  mean(grepl("^ENSG\\d+|^ENSMUSG\\d+|^ENSRNOG\\d+|^ENSDOG\\d+|^Y[A-P][LR]\\d+", ids)) > 0.5
}

#' Get OrgDb object from species name
#' @param species Character string
#' @export
get_orgdb <- function(species) {
  switch(species,
         "Homo sapiens" = org.Hs.eg.db::org.Hs.eg.db,
         "Mus musculus" = org.Mm.eg.db::org.Mm.eg.db,
       #  "Rattus norvegicus" = org.Rn.eg.db::org.Rn.eg.db,
      #   "Saccharomyces cerevisiae" = org.Sc.sgd.db::org.Sc.sgd.db,
     #    "Canis familiaris" = org.Cf.eg.db::org.Cf.eg.db,
         stop("Unsupported species: ", species))
}

#' Get KEGG Code
#' @param species Character string
#' @export
get_kegg_code <- function(species) {
  switch(species,
         "Homo sapiens" = "hsa",
         "Mus musculus" = "mmu",
         "Rattus norvegicus" = "rno",
         "Saccharomyces cerevisiae" = "sce",
         "Canis familiaris" = "cfa",
         stop("Unsupported species for KEGG: ", species))
}

#' Get Reactome Code
#' @param species Character string
#' @export
get_reactome_code <- function(species) {
  switch(species,
         "Homo sapiens" = "human",
         "Mus musculus" = "mouse")
         #"Rattus norvegicus" = "rat",
         #"Canis familiaris" = "dog",
         # Yeast not supported by Reactome
         #"Saccharomyces cerevisiae" = "human",  # fallback
}

#' Check if IDs are gene symbols
#' @param ids Character vector of gene symbols
#' @export
is_symbol <- function(ids) {
  mean(!grepl("^ENS", ids)) > 0.5
}
#' Sanitize Ensembl Gene IDs
#'
#' This function standardizes Ensembl gene identifiers by removing version suffixes
#' (e.g., ".1") and trimming composite IDs that include gene symbols (e.g., "ENSG00000141510_TP53").
#'
#' These modifications help ensure compatibility with annotation databases and enrichment tools.
#'
#' @param ids A character vector of Ensembl gene IDs, potentially with version suffixes
#' or appended gene symbols.
#'
#' @return A character vector of cleaned Ensembl gene IDs.
#'
#' @examples
#' sanitize_ensembl_ids(c("ENSG00000141510.1", "ENSG00000141510_TP53", "ENSG00000141510.2_TP53"))
#' # Returns: "ENSG00000141510" "ENSG00000141510" "ENSG00000141510"
#'
#' @export
sanitize_ensembl_ids <- function(ids) {
  ids <- as.character(ids)
  ids <- gsub("\\.\\d+$", "", ids)         # remove version suffix
  ids <- gsub("_.*$", "", ids)             # remove gene symbol suffix
  ids
}

