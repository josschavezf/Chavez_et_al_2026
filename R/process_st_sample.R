process_st_sample <- function(g,
                              feat_det_in_min_cells = 50,
                              min_det_feats_per_cell = 100) {
    g <- filterGiotto(g,
                      feat_det_in_min_cells = feat_det_in_min_cells,
                      min_det_feats_per_cell = min_det_feats_per_cell)
    
    g <- normalizeGiotto(g)
    g <- calculateHVF(g)
    g <- runPCA(g)
    g <- runUMAP(g)
    g <- createNearestNetwork(g)
    g <- doLeidenCluster(g)
    
    return(g)
}
