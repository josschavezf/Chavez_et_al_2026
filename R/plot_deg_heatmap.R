plot_deg_heatmap <- function(gobject, 
                             features,
                             cluster_column,
                             values = "scaled",
                             return_zscores_matrix = FALSE,
                             ...) {
    # select features
    my_expression <- getExpression(gobject,
                                   values = values,
                                   output = "matrix")
    
    my_expression <- my_expression[features, ]
    
    # get gene names
    my_rownames <- rownames(my_expression)
    
    # calculate zscore
    my_zscores <- data.frame(cell_ids = colnames(my_expression))
    
    for (gene in my_rownames) {
        mean_expression_gene <- mean(my_expression[gene, ])
        sd_expression_gene <- stats::sd(my_expression[gene, ])
        
        my_zscores[[gene]] <- (my_expression[gene, ] - mean_expression_gene)/sd_expression_gene
    }
    
    my_zscores_matrix <- t(as.matrix(my_zscores[,-1]))
    colnames(my_zscores_matrix) <- my_zscores$cell_ids
    
    # calculate mean zscore per polygon
    my_zscores_mean <- data.table::data.table(feat_ID = my_rownames)
    
    my_metadata <- pDataDT(gobject)
    
    for (i in unique(my_metadata[[cluster_column]])) {
        my_cells <- my_metadata[my_metadata[[cluster_column]] == i, "cell_ID"]
        my_sub_zscores <- my_zscores_matrix[, my_cells$cell_ID]
        mean_zscores <- Matrix::rowMeans(my_sub_zscores)
        my_zscores_mean <- cbind(my_zscores_mean, mean_zscores)
    }
    
    my_zscores_mean <- as.matrix(my_zscores_mean[, -1])
    colnames(my_zscores_mean) <- unique(my_metadata[[cluster_column]])
    rownames(my_zscores_mean) <- my_rownames
    
    if(return_zscores_matrix == TRUE) {
        return(my_zscores_mean)
    } else {
        # plot heatmap
        my_heatmap <- ComplexHeatmap::Heatmap(my_zscores_mean,
                                              heatmap_legend_param = list(title = paste0(values, "\nmean z score")),
                                              show_row_dend = FALSE,
                                              show_column_dend = FALSE,
                                              ...)
        
        return(my_heatmap)
    }
    
}
