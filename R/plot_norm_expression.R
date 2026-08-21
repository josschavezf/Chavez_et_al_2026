plot_norm_expression <- function(g, features, cluster_column,
                                 ...) {
  g_sc_expression <- getExpression(g,
                                   values = "normalized",
                                   output = "matrix")
  
  cell_metadata <- pDataDT(g)
  
  df <- data.frame(genes = features)
  
  for(i in unique(cell_metadata[[cluster_column]])) {
    sub_cells <- cell_metadata[cell_metadata[[cluster_column]] == i,]$cell_ID
    g_sc_expression_sub <- g_sc_expression[features, sub_cells]
    df[[i]] <- rowMeans(as.matrix(g_sc_expression_sub))
  }
  
  m <- as.matrix(df[, -1])
  rownames(m) <- df$genes
  
  Heatmap(matrix = m,
          show_column_dend = FALSE,
          show_row_dend = FALSE,
          heatmap_legend_param = list(title = "normalized\nexpression"),
          ...)
}
