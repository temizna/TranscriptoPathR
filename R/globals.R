utils::globalVariables(c("pathway_input_rv"))
utils::globalVariables(c("dds", "category", "ENTREZID", "log2_Expression"))
 # Observe clear data button click and clear data
  observeEvent(input$clear_data, {
    # Clear all data (counts, samples, and design)
    data$counts <- NULL
    data$samples <- NULL
    data$norm_counts <- NULL
    data$design_matrix <- NULL
    dds_rv(NULL)  # Clear DESeq2 results
    pathway_data(NULL)  # Clear pathway data
    
    # Optionally, update the UI to show a message
    output$clear_status <- renderText("All data has been cleared!")
  })
  
  # Observe clear DESeq2 data button click and clear dds
  observeEvent(input$clear_dds, {
    dds_rv(NULL)  # Clear DESeq2 data
    
    output$clear_status <- renderText("DESeq2 data has been cleared!")
  })
  
  # Observe clear pathway data button click and clear pathway data
  observeEvent(input$clear_pathway, {
    pathway_data(NULL)  # Clear pathway data
    
    output$clear_status <- renderText("Pathway data has been cleared!")
  })
