server<-function (input, output, session){
#mod_geo_loader(input, output, session, data)
utils::globalVariables(c("value", "variable", "V1", "V2"))
  
gsea_result_rv <- reactiveVal()
res_reactive <- reactiveVal()
heatmap_data_reactive <- reactiveVal()
geneList_rv <- reactiveVal()
m_list_rv <- reactiveVal()
kegg_png_rv <- reactiveVal()
pathway_input_rv <- reactiveVal()
d1_merged_rv <- reactiveVal()
kegg_pathway_results <- reactiveVal()
pathway_result_rv<- reactiveVal()
# Create shared reactive values

dds_rv <- reactiveVal()
loaded_data_rv <-reactiveValues()
filtered_data_rv <- reactiveValues()
filtered_dds_rv <- reactiveVal()
geneListU_rv <-reactiveVal()
# Register module
mod_data_load_design(input, output, session, loaded_data_rv,dds_rv)
callModule(mod_sample_select_dynamic, "sample_select",
           dds_rv = dds_rv,
           loaded_data_rv = loaded_data_rv,
           filtered_data_rv = filtered_data_rv,
           filtered_dds_rv = filtered_dds_rv)
mod_gene_expression_plot(input, output, session,filtered_data_rv)
mod_qc_plot(input, output, session, filtered_data_rv)
de_sel <- mod_de_server(
  id = "de",
  filtered_data_rv = filtered_data_rv,
  filtered_dds_rv  = filtered_dds_rv,
  res_reactive     = res_reactive
)

mod_gsva_server(
  id = "gsva",
  dds_rv = dds_rv,
  de_sel = de_sel,
  filtered_data_rv = filtered_data_rv,
  group_var_react  = de_sel$group_var,
  ref_level_react  = de_sel$ref_level,
  test_level_react = de_sel$test_level
)

mod_volcano_ma_plot(input, output, session, res_reactive, filtered_data_rv)
mod_cross_server("cross", filtered_data_rv = filtered_data_rv, filtered_dds_rv = filtered_dds_rv)
mod_pathway_server(
  id = "pathway",
  filtered_data_rv     = filtered_data_rv,
  res_reactive         = res_reactive,
  geneList_rv          = geneList_rv,
  kegg_pathway_results = kegg_pathway_results,
  d1_merged_rv         = d1_merged_rv,
  de_sel = de_sel,
  pathway_result_rv    = pathway_result_rv
)
mod_gsea_server("gsea", filtered_data_rv, res_reactive,  de_sel = de_sel)
mod_pathway_plots_server(
  id = "pathplots",
  pathway_result_rv = pathway_result_rv,
  geneList_rv       = geneList_rv,
  de_sel            = de_sel   # from mod_de_server(...)
)

mod_cancer_gene_census(input, output, session, res_reactive, filtered_data_rv)
mod_nonoverlap_server(
  id = "nonOL",
  filtered_data_rv   = filtered_data_rv,
  pathway_result_rv  = pathway_result_rv,  # from the primary pathway module
  geneList_rv        = geneList_rv,        # original ENTREZID-named vector
  geneListU_rv       = geneListU_rv        # will be filled by this module
)
mod_tf_enrichment_server("tf", res_reactive = res_reactive, filtered_data_rv = filtered_data_rv)

mod_logger_server("logger", input_all = input)
mod_pca_cluster_server(input, output, session, filtered_data_rv)
mod_goi_heatmap(input, output, session, filtered_data_rv)
            # or however you store it

onStop(function() {
    cat("Shiny app is closing cleanly...\n")
    stopApp()
  })
}
