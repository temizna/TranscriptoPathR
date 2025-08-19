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
mod_differential_expression(input, output, session, filtered_data_rv, filtered_dds_rv, res_reactive)
mod_volcano_ma_plot(input, output, session, res_reactive, filtered_data_rv)
mod_cross_plot(input, output, session, filtered_data_rv,filtered_dds_rv)
mod_pathway_analysis(input, output, session, filtered_data_rv, res_reactive, geneList_rv, kegg_pathway_results, d1_merged_rv, pathway_result_rv)
mod_gsea_analysis(input, output, session, filtered_data_rv, res_reactive)
mod_cancer_gene_census(input, output, session, res_reactive, filtered_data_rv)
mod_pathway_analysis_non_overlap(input, output, session, geneList_rv, filtered_data_rv, pathway_result_rv, res_reactive,geneListU_rv)
mod_tf_enrichment_analysis(input, output, session, res_reactive,filtered_data_rv) 
mod_logger_server("logger", input_all = input)
mod_pca_cluster_server(input, output, session, filtered_data_rv)
mod_goi_heatmap(input, output, session, filtered_data_rv)
onStop(function() {
    cat("Shiny app is closing cleanly...\n")
    stopApp()
  })
}
