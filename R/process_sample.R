process_sample <- function(g, 
                           feat_det_in_min_cells = 100) {
    g <- filterGiotto(g,
                      feat_det_in_min_cells = feat_det_in_min_cells,
                      min_det_feats_per_cell = 500)
    
    
    g <- normalizeGiotto(g)
    g <- addStatistics(g)
    g <- calculateHVF(g)
    g <- runPCA(g)
    g <- runUMAP(g)
    g <- createNearestNetwork(g)
    g <- doLeidenCluster(g)
    
    return(g)
    
}
