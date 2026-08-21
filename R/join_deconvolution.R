library(tidyverse)

join_deconvolution <- function(gobject_list) {
  
  result <- data.frame(cell_ID = character(),
                       cell_type = character(),
                       deconvolution_proportion = numeric())
  
  for (i in names(gobject_list)) {
    x <- getSpatialEnrichment(gobject_list[[i]],
                              name = "RCTD",
                              output = "data.table")
    
    x$cell_ID <- paste0(i, "-", x$cell_ID)
    
    x_long <- x %>% 
      pivot_longer(-cell_ID,
                   names_to = "cell_type",
                   values_to = "deconvolution_proportion")
    
    result <- rbind(result, x_long)
  }
  
  return(result)
}
