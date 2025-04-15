server<-function (input, output, session){
mod_geo_loader(input, output, session, data)
gsea_result_rv <- reactiveVal()
res_reactive <- reactiveVal()
heatmap_data_reactive <- reactiveVal()
geneList_rv <- reactiveVal()
m_list_rv <- reactiveVal()
kegg_png_rv <- reactiveVal()
pathway_input_rv <- reactiveVal()
d1_merged_rv <- reactiveVal()
kegg_pathway_results <- reactiveVal()
# Create shared reactive values
data <- reactiveValues(
  counts = NULL, 
  samples = NULL, 
  norm_counts = NULL, 
  design_matrix = NULL, 
  species = "Homo sapiens"
)
dds_rv <- reactiveVal()

# Register module
mod_data_upload_design(input, output, session, data = data, dds_rv = dds_rv)
mod_sample_select(input, output, session, data)
mod_gene_expression_plot(input, output, session, data)
mod_qc_plot(input, output, session, data)

  
  # Load modules
mod_differential_expression(input, output, session, data, dds_rv, res_reactive)
mod_volcano_ma_plot(input, output, session, res_reactive, data)
mod_cross_plot(input, output, session, data,dds_rv)
mod_pathway_analysis(input, output, session, data, res_reactive, geneList_rv, pathway_input_rv, kegg_pathway_results, d1_merged_rv)
mod_gsea_analysis(input, output, session, data, res_reactive)
onStop(function() {
    cat("Shiny app is closing cleanly...\n")
    stopApp()
  })
}
