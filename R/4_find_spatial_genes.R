library(Giotto)
instructions <- createGiottoInstructions()

##%######################################################%##
#                                                          #
####                    Read samples                    ####
#                                                          #
##%######################################################%##

HTA12_129 <- readRDS("data/deconvoluted_visium/brca/HTA12_129.RDS")
HTA12_135 <- readRDS("data/deconvoluted_visium/brca/HTA12_135.RDS")
HTA12_139 <- readRDS("data/deconvoluted_visium/brca/HTA12_139.RDS")
HTA12_142 <- readRDS("data/deconvoluted_visium/brca/HTA12_142.RDS")
HTA12_144 <- readRDS("data/deconvoluted_visium/brca/HTA12_144.RDS")
HTA12_145 <- readRDS("data/deconvoluted_visium/brca/HTA12_145.RDS")
HTA12_246 <- readRDS("data/deconvoluted_visium/brca/HTA12_246.RDS")
HTA12_248 <- readRDS("data/deconvoluted_visium/brca/HTA12_248.RDS")
HTA12_249 <- readRDS("data/deconvoluted_visium/brca/HTA12_249.RDS")
HTA12_250 <- readRDS("data/deconvoluted_visium/brca/HTA12_250.RDS")
HTA12_251 <- readRDS("data/deconvoluted_visium/brca/HTA12_251.RDS")
HTA12_253 <- readRDS("data/deconvoluted_visium/brca/HTA12_253.RDS")
HTA12_262 <- readRDS("data/deconvoluted_visium/brca/HTA12_262.RDS")
HTA12_266 <- readRDS("data/deconvoluted_visium/brca/HTA12_266.RDS")
HTA12_271 <- readRDS("data/deconvoluted_visium/brca/HTA12_271.RDS")

##%######################################################%##
#                                                          #
####                    Join objects                    ####
#                                                          #
##%######################################################%##

g <- joinGiottoObjects(
  gobject_list = list(HTA12_129, HTA12_135, HTA12_139, HTA12_142,
                      HTA12_144, HTA12_145, HTA12_246, HTA12_248,
                      HTA12_249, HTA12_250, HTA12_251, HTA12_253, 
                      HTA12_262, HTA12_266, HTA12_271),
  gobject_names = c("HTA12_129", "HTA12_135", "HTA12_139", "HTA12_142",
                    "HTA12_144", "HTA12_145", "HTA12_246", "HTA12_248",
                    "HTA12_249", "HTA12_250", "HTA12_251", "HTA12_253", 
                    "HTA12_262", "HTA12_266", "HTA12_271"))

##%######################################################%##
#                                                          #
####                   Normalization                    ####
#                                                          #
##%######################################################%##

g <- normalizeGiotto(g)

g <- createSpatialNetwork(g,
                          method = "kNN", 
                          k = 6,
                          maximum_distance_knn = 400,
                          name = "spatial_network")

##%######################################################%##
#                                                          #
####             Plot cell type enrichment              ####
#                                                          #
##%######################################################%##

marker_genes <- list(
  Adipocytes = c("PLIN1", "FABP4", "GLUT4", "HOXC8", "HOXC9",
                 "UCP1", "CIDEA", "PRDM16", "CD137", "TMEM26",
                 "GUCY1B1", "SLC35F1"),
  B_cell = c("CD79A", "CD79B", "CD19", "MS4A1", "IGHD", "CD22", "CD52",
             "CD20", "CD23", "CD27", "PAX5", "CD40", "CD38", "GL7",
             "RUBCNL", "COL26A1", "RGS7"),
  Basal = c("SAMD5", "FBXO32", "TP63", "RBBP8", "KLHL13"),
  Dendritic_cells = c("CADM1", "XCR1", "CLEC9A", "RAB32", "C1orf54",
                      "CD1C", "FCER1A", "CLEC10A", "CD1E",
                      "LAMP3", "CCR7", "FSCN1", "CD83", "CCL22"),
  Endothelial = c("EMCN", "FLT1", "PECAM1", "VWF", "PTPRB", "ACTA2", "ANGPT2"),
  Epithelial = c("LumSec", "ELF5", "CL28", "KRT15", "BARX2", "HS3ST4"),
  Fibroblasts = c("COL1A1", "COL3A1", "COL5A1", "LUM", "MMP2"),
  Luminal = c("ANKRD30A", "ERBB4", "AFF3", "TTC6", "ESR1", "NEK10", "XBP1"),
  Luminal_progenitor = c("K15", "S100A7", "S100A8", "S100A9", "ALDH1A3", 
                         "CD49b", "ITGA2", "CD117", "GABRP", "CCL28"),
  Macrophages = c("CX3CR1", "CD80", "CD86", "CD163", "MSR1"),
  Mast =  c("HPGD", "TPSB2", "HDC", "SLC18A2", "CPA3", "SLC8A3", "RGS13", "CPA3"),
  pDC = c("IL3RA", "BCL11A", "CLEC4C", "NRP1"),
  Plasma = c("CD38", "SDC1", "IGHG1", "IGKC", "MZB1"),
  T_cell = c("IL7R", "CD4", "CD8A", "CD8B", "CD3G", "CD3D", "CD3E",
             "FCGR3A", "GZMA", "NCAM1", "IL2RA", "CTLA4", "FOXP3", 
             "TNFRSF18", "IKZF2", "CD96", "THEMIS", "GRAP2"),
  Pericytes = c("RGS5", "PLXDC1", "FN1", "MCAM"),
  Tumor = c("KRT7", "CD24", "KRT19", "TACSTD2", "KRT18", "KRT8", "SPINT2", 
            "PERP", "TRPS1",  "RIMS2", "NFIB")
)

sign_matrix <- makeSignMatrixPAGE(sign_names = names(marker_genes),
                                  sign_list = marker_genes)

g <- runPAGEEnrich(gobject = g, 
                   sign_matrix = sign_matrix,
                   min_overlap_genes = 3,
                   name = "PAGE_marker_genes")

cell_type_enrichment <- getSpatialEnrichment(g,
                                             name = "PAGE_marker_genes",
                                             output = "data.table")

spatial_locs <- getSpatialLocations(g,
                                    output = "data.table")

cell_metadata <- pDataDT(g)
cell_type_enrichment <- merge(cell_type_enrichment, cell_metadata[, c("cell_ID", "list_ID")])
cell_type_enrichment <- merge(cell_type_enrichment, spatial_locs)

ggplot(cell_type_enrichment, aes(sdimx, sdimy, colour = Tumor)) +
  geom_point(size = 0.5) +
  facet_wrap(vars(list_ID), scales = "free") +
  scale_color_gradient2(low = "blue", high = "red") +
  labs(color = "") +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank())
ggsave("plots_paper/4_enrichment_Tumor.png", device = "png", 
       width = 5, height = 5, scale = 2)

ggplot(cell_type_enrichment, aes(sdimx, sdimy, colour = Fibroblasts)) +
  geom_point(size = 0.5) +
  facet_wrap(vars(list_ID), scales = "free") +
  scale_color_gradient2(low = "blue", high = "red") +
  labs(color = "") +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank())
ggsave("plots_paper/4_enrichment_Fibroblasts.png", device = "png", 
       width = 5, height = 5, scale = 2)

ggplot(cell_type_enrichment, aes(sdimx, sdimy, colour = Macrophages)) +
  geom_point(size = 0.5) +
  facet_wrap(vars(list_ID), scales = "free") +
  scale_color_gradient2(low = "blue", high = "red") +
  labs(color = "") +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank())
ggsave("plots_paper/4_enrichment_Macrophages.png", device = "png", 
       width = 5, height = 5, scale = 2)

ggplot(cell_type_enrichment, aes(sdimx, sdimy, colour = T_cell)) +
  geom_point(size = 0.5) +
  facet_wrap(vars(list_ID), scales = "free") +
  scale_color_gradient2(low = "blue", high = "red") +
  labs(color = "") +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank())
ggsave("plots_paper/4_enrichment_T_cell.png", device = "png", 
       width = 5, height = 5, scale = 2)

##%######################################################%##
#                                                          #
####           Find spatial genes per sample            ####
#                                                          #
##%######################################################%##

g_list <- list(HTA12_129, HTA12_135, HTA12_139, HTA12_142,
               HTA12_144, HTA12_145, HTA12_246, HTA12_248,
               HTA12_249, HTA12_250, HTA12_251, HTA12_253, 
               HTA12_262, HTA12_266, HTA12_271)

names(g_list) <- c("HTA12_129", "HTA12_135", "HTA12_139", "HTA12_142",
                   "HTA12_144", "HTA12_145", "HTA12_246", "HTA12_248",
                   "HTA12_249", "HTA12_250", "HTA12_251", "HTA12_253", 
                   "HTA12_262", "HTA12_266", "HTA12_271")

ext_spatial_genes_list <- list()
ranktest_list <- list()

for(i in names(g_list)) {
  g_list[[i]] <- createSpatialNetwork(g_list[[i]],
                                      method = "kNN", 
                                      k = 6,
                                      maximum_distance_knn = 400,
                                      name = "spatial_network")
  
  # Select the top spatial genes
  ranktest_list[[i]] <- binSpect(g_list[[i]], 
                                 bin_method = "rank",
                                 calc_hub = TRUE, 
                                 hub_min_int = 6,
                                 spatial_network_name = "spatial_network")
}

# save spatial genes for testing
saveRDS(ranktest_list, "ranktest_list.RDS")

# Take the first 20 spatial genes per sample
for(i in names(g_list)) {
  ext_spatial_genes_list[[i]] <- ranktest_list[[i]][1:20,]$feats
}

ext_spatial_genes_df <- as.data.frame(ext_spatial_genes_list)
writexl::write_xlsx(ext_spatial_genes_df, "supplementary_table_1.xlsx") # supp_table_1, sheet 2
ext_spatial_genes <- unique(unlist(ext_spatial_genes_df))

# Remove MT, RP, and non-coding genes
ext_spatial_genes <- ext_spatial_genes[!ext_spatial_genes %in% grep("^MT-", ext_spatial_genes, value = TRUE)]
ext_spatial_genes <- ext_spatial_genes[!ext_spatial_genes %in% grep("^RP", ext_spatial_genes, value = TRUE)]
ext_spatial_genes <- ext_spatial_genes[!ext_spatial_genes %in% grep("\\.", ext_spatial_genes, value = TRUE)]

##%######################################################%##
#                                                          #
####    Cluster the top spatial genes into clusters     ####
#                                                          #
##%######################################################%##

# Calculate pairwise distances between genes
spat_cor_netw_DT <- detectSpatialCorFeats(
  g,
  method = "network",
  spatial_network_name = "spatial_network",
  subset_feats = ext_spatial_genes)

# Create 4 clusters
spat_cor_netw_DT <- clusterSpatialCorFeats(spat_cor_netw_DT, 
                                           name = "spat_netw_clus", 
                                           k = 4)

heatmSpatialCorFeats(g,
                     spatCorObject = spat_cor_netw_DT,
                     use_clus_name = "spat_netw_clus",
                     heatmap_legend_param = list(title = NULL))

##%######################################################%##
#                                                          #
####        Calculate PAGE enrichment per module        ####
#                                                          #
##%######################################################%##

coexpr_dt <- data.table::data.table(
  genes = names(spat_cor_netw_DT$cor_clusters$spat_netw_clus),
  cluster = spat_cor_netw_DT$cor_clusters$spat_netw_clus)

data.table::setorder(coexpr_dt, cluster)

sign_matrix <- makeSignMatrixPAGE(sign_names = paste0("cluster_", 1:4),
                                  sign_list = list(coexpr_dt[cluster == 1]$genes,
                                                   coexpr_dt[cluster == 2]$genes,
                                                   coexpr_dt[cluster == 3]$genes,
                                                   coexpr_dt[cluster == 4]$genes))

g <- runPAGEEnrich(gobject = g, 
                   sign_matrix = sign_matrix,
                   name = "PAGE_spatial_genes")

##%######################################################%##
#                                                          #
####    Heatmap sp clusters vs cell type enrichment     ####
#                                                          #
##%######################################################%##

page_metagene_enrichment <- getSpatialEnrichment(g,
                                                 name = "PAGE_spatial_genes",
                                                 output = "data.table")

cell_type_enrichment <- getSpatialEnrichment(g,
                                             name = "PAGE_marker_genes",
                                             output = "data.table")

joined_enrichment <- merge(page_metagene_enrichment, cell_type_enrichment)

lm_coefficients <- data.frame()
for(cluster in c("cluster_1", "cluster_2", "cluster_3", "cluster_4")) {
  for(cell_type in c("Adipocytes", "B_cell", "Basal", "Dendritic_cells", 
                     "Endothelial", "Epithelial", "Fibroblasts", "Luminal", 
                     "Luminal_progenitor", "Macrophages", "Mast",
                     "Pericytes", "Plasma", "pDC", "T_cell", "Tumor")) {
    
    cor_value <- cor(x = joined_enrichment[[cell_type]],
                     y = joined_enrichment[[cluster]],
                     method = "pearson")
    
    lm_result <- lm(joined_enrichment[[cluster]]~joined_enrichment[[cell_type]])
    sum_lm_result <- summary(lm_result)
    
    lm_coefficients <- rbind(lm_coefficients,
                             data.frame("cluster" = cluster, 
                                        "cell_type" = cell_type,
                                        "cor_value" = cor_value,
                                        "r2" = round(sum_lm_result$r.squared, 8)))
  }
}

lm_coefficients_cor <- lm_coefficients %>% 
  dplyr::select(cluster, cell_type, cor_value) %>% 
  pivot_wider(names_from = cluster, values_from = cor_value)

m <- as.matrix(lm_coefficients_cor[, -1])
rownames(m) <- lm_coefficients_cor$cell_type

Heatmap(m,
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        heatmap_legend_param = list(title = ""))

pdf("plots_paper/4_sp_celltype_cor_genes.pdf", width = 5, height = 5)
Heatmap(m,
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        heatmap_legend_param = list(title = ""))
dev.off()

##%######################################################%##
#                                                          #
####          Plot spatial cluster enrichment           ####
#                                                          #
##%######################################################%##

page_metagene_enrichment <- getSpatialEnrichment(g,
                                                 name = "PAGE_spatial_genes",
                                                 output = "data.table")

spatial_locs <- getSpatialLocations(g,
                                    output = "data.table")

cell_metadata <- pDataDT(g)

page_metagene_enrichment <- merge(page_metagene_enrichment, 
                                  spatial_locs)

page_metagene_enrichment <- merge(page_metagene_enrichment,
                                  cell_metadata)

ggplot(page_metagene_enrichment, aes(x = sdimx, y = sdimy, color = cluster_1)) +
  geom_point(size = 0.5) +
  facet_wrap(vars(list_ID), scales = "free") +
  labs(title = "Cluster 1",
       color = "",
       x = "", y = "") +
  scale_color_gradient2(low = "blue", high = "red") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank())
ggsave("plots_paper/4_spatial_cluster_1.pdf", device = "pdf",
       width = 5, height = 5, scale = 2)

##%######################################################%##
#                                                          #
####               Plot spatial clusters                ####
####           enrichment in scRNAseq samples           ####
#                                                          #
##%######################################################%##

spatial_genes <- readxl::read_xlsx("supplementary_table_1.xlsx", sheet = 3)

sign_matrix <- makeSignMatrixPAGE(
  sign_names = paste0("cluster_", 1:4),
  sign_list = list(spatial_genes[spatial_genes$cluster == 1,]$genes,
                   spatial_genes[spatial_genes$cluster == 2,]$genes,
                   spatial_genes[spatial_genes$cluster == 3,]$genes,
                   spatial_genes[spatial_genes$cluster == 4,]$genes))

g_sc <- runPAGEEnrich(gobject = g_sc, 
                   sign_matrix = sign_matrix,
                   name = "PAGE_spatial_genes")

page_enrichment_sc <- getSpatialEnrichment(g_sc,
                                           name = "PAGE_spatial_genes",
                                           output = "data.table")

cell_metadata <- pDataDT(g_sc)
cell_metadata <- cell_metadata[cell_metadata$cell_types_final == "Tumor", ]

sp_enrichment_sc <- data.frame()
for (i in sort(unique(cell_metadata$list_ID))) {
  cell_ids <- cell_metadata[cell_metadata$list_ID == i,]$cell_ID
  page_enrichment_sc_small <- page_enrichment_sc[page_enrichment_sc$cell_ID %in% cell_ids,]
  mean_enrichment_sc <- colMeans(page_enrichment_sc_small[,-1])
  sp_enrichment_sc <- rbind(sp_enrichment_sc, mean_enrichment_sc)
}
colnames(sp_enrichment_sc) <- colnames(page_enrichment_sc)[-1]
rownames(sp_enrichment_sc) <- sort(unique(cell_metadata$list_ID))

sp_enrichment_sc$cluster_2_minus_3 <- sp_enrichment_sc$cluster_2 - sp_enrichment_sc$cluster_3

m <- as.matrix(sp_enrichment_sc)

pdf("plots_paper/4_sp_in_sc.pdf", width = 6, height = 3)
Heatmap(t(m),
        row_order = c("cluster_1", "cluster_2", "cluster_3", "cluster_4", "cluster_2_minus_3"),
        column_order = c("HTA12_144", "HTA12_139", "HTA12_145", "HTA12_248",
                         "HTA12_135", "HTA12_253", "HTA12_249", "HTA12_246",
                         "HTA12_266", "HTA12_251", "HTA12_250", "HTA12_129",
                         "HTA12_262", "HTA12_271", "HTA12_142"),
        heatmap_legend_param = list(title = "Cluster enrichment"))
dev.off()

##%######################################################%##
#                                                          #
####                  EnrichR pathways                  ####
#                                                          #
##%######################################################%##

pathways_1 <- read.delim("MSigDB_Hallmark_2020_table_module_1.txt")
pathways_2 <- read.delim("MSigDB_Hallmark_2020_table_module_2.txt")
pathways_3 <- read.delim("MSigDB_Hallmark_2020_table_module_3.txt")
pathways_4 <- read.delim("MSigDB_Hallmark_2020_table_module_4.txt")

pathways_1$module <- "M1"
pathways_2$module <- "M2"
pathways_3$module <- "M3"
pathways_4$module <- "M4"

all_pathways <- rbind(pathways_1[1:5, c("Term", "Adjusted.P.value", "module", "Odds.Ratio")],
                      pathways_2[1:5, c("Term", "Adjusted.P.value", "module", "Odds.Ratio")],
                      pathways_3[1:5, c("Term", "Adjusted.P.value", "module", "Odds.Ratio")],
                      pathways_4[1:5, c("Term", "Adjusted.P.value", "module", "Odds.Ratio")])

ggplot(all_pathways, aes(x = module, y = Term, size = -log10(Adjusted.P.value),
                         color = Odds.Ratio)) +
  geom_point() +
  scale_y_discrete(limits = rev(unique(all_pathways$Term))) +
  theme_bw() +
  labs(x = "",
       y = "")
ggsave("plots_paper/4_pathways.pdf", device = "pdf",
       width = 5.5, height = 4)
