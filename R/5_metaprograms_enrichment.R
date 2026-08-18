library(Giotto)
library(cluster)
library(factoextra)

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
                      HTA12_249, HTA12_250, HTA12_251, HTA12_253, HTA12_262,
                      HTA12_266, HTA12_271),
  gobject_names = c("HTA12_129", "HTA12_135", "HTA12_139", "HTA12_142",
                    "HTA12_144", "HTA12_145", "HTA12_246", "HTA12_248",
                    "HTA12_249", "HTA12_250", "HTA12_251", "HTA12_253", "HTA12_262",
                    "HTA12_266", "HTA12_271"))

##%######################################################%##
#                                                          #
####          Calculate metaprogram enrichment          ####
#                                                          #
##%######################################################%##

metaprograms <- readxl::read_xlsx("mp_table.xlsx")
metaprograms <- janitor::clean_names(metaprograms)
metaprograms <- metaprograms[-41]

sign_matrix <- makeSignMatrixPAGE(sign_names = colnames(metaprograms), 
                                  sign_list = as.list(metaprograms))

g <- runPAGEEnrich(g,
                   spat_unit = "cell",
                   feat_type = "rna",
                   sign_matrix)

mp_enrichment <- getSpatialEnrichment(g,
                                      name = "PAGE",
                                      output = "data.table")

##%######################################################%##
#                                                          #
####                Cluster metaprograms                ####
#                                                          #
##%######################################################%##

df <- as.data.frame(t(mp_enrichment[,-1]))

## Remove non-BRCA-associated metaprograms
remove_mp <- c("mp37_hemato_related_ii", 
               "mp27_oligo_progenitor",
               "mp35_hemato_related_i",
               "mp16_mes_glioma",
               "mp29_npc_opc",
               "mp28_oligo_normal",
               "mp25_astrocytes",
               "mp31_alveolar",
               "mp32_skin_pigmentation",
               "mp38_glutathione",
               "mp22_secreted_i",
               "mp39_metal_response",
               "mp40_pdac_related",
               "mp34_platelet_activation",
               "mp26_npc_glioma",
               "mp24_cilia",
               "mp23_secreted_ii",
               "mp33_rb_cs",
               "mp30_pdac_classical",
               "mp15_emt_iv")
df <- df[!rownames(df) %in% remove_mp,]

# Perform K-means clustering
km_result <- kmeans(df, centers = 4, nstart = 25)

# Calculate the distance matrix
dist_matrix <- dist(df)

# Compute silhouette coefficients for each observation
sil_info <- silhouette(km_result$cluster, dist_matrix)

# Visualization using factoextra for a ggplot2-based plot
fviz_silhouette(sil_info)
fviz_nbclust(df, kmeans, method = "silhouette")

s_results <- as.data.frame(sil_info)
rownames(s_results) <- rownames(df)
View(s_results[s_results$cluster == 1,])
View(s_results[s_results$cluster == 2,])
View(s_results[s_results$cluster == 3,])
View(s_results[s_results$cluster == 4,])

##%######################################################%##
#                                                          #
####    Plot mean enrichment per metaprogram cluster    ####
#                                                          #
##%######################################################%##

## using clusters from k = 4

### cluster 1
mp_means_1 <-  mp_enrichment$mp12_emt_i

### cluster 2
mp_means_2 <- rowMeans(mp_enrichment[,c("mp17_interferon_mhc_ii_i", 
                                     "mp18_interferon_mhc_ii_ii", 
                                     "mp36_ig")])

### cluster 3
mp_means_3 <- rowMeans(mp_enrichment[,c("mp21_respiration", 
                                     "mp3_cell_cylce_hmg_rich", 
                                     "mp14_emt_iii")])

### cluster 4
mp_means_4 <- rowMeans(mp_enrichment[,c("mp2_cell_cycle_g1_s", 
                                     "mp7_stress_in_vitro", 
                                     "mp1_cell_cycle_g2_m", 
                                     "mp6_hypoxia",
                                     "mp4_chromatin",
                                     "mp11_translation_initiation",
                                     "mp9_unfolded_protein_response",
                                     "mp19_epithelial_senescence",
                                     "mp20_myc",
                                     "mp5_stress",
                                     "mp13_emt_ii",
                                     "mp8_proteasomal_degradation",
                                     "mp10_protein_maturation")])

spatial_locs <- as.data.frame(getSpatialLocations(g,
                                                  output = "data.table"))

cell_metadata <- pDataDT(g)

spatial_locs$mp_means <- mp_means_4 # repeat for each cluster
spatial_locs$list_ID <- cell_metadata$list_ID

ggplot(spatial_locs, aes(sdimx, sdimy, colour = mp_means)) +
  geom_point(size = 0.5) +
  facet_wrap(vars(list_ID), scales = "free") +
  labs(x = "", y = "", color = "") +
  scale_color_gradient2(low = "blue", high = "red") +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank())
ggsave("plots_paper/5_mp_cluster4.png", width = 5, height = 5, scale = 2)

##%######################################################%##
#                                                          #
####           Join deconvolution proportions           ####
#                                                          #
##%######################################################%##

g_list <- list(HTA12_129, HTA12_135, HTA12_139,
               HTA12_142, HTA12_144, HTA12_145, HTA12_246, HTA12_248,
               HTA12_249, HTA12_250, HTA12_251, HTA12_253, HTA12_262,
               HTA12_266, HTA12_271)

names(g_list) <- c("HTA12_129", "HTA12_135", "HTA12_139",
                   "HTA12_142", "HTA12_144", "HTA12_145", "HTA12_246", "HTA12_248",
                   "HTA12_249", "HTA12_250", "HTA12_251", "HTA12_253", "HTA12_262",
                   "HTA12_266", "HTA12_271")

source("R/join_deconvolution.R")

g_deconvolution <- join_deconvolution(g_list)

##%######################################################%##
#                                                          #
####             Linear regression model of             ####
####        clusters enrichment vs deconvolution        ####
#                                                          #
##%######################################################%##

mp_clusters_enrichment <- data.frame(cell_ID  = mp_enrichment$cell_ID,
                                     cluster1 = mp_means_1$mp12_emt_i,
                                     cluster2 = mp_means_2,
                                     cluster3 = mp_means_3,
                                     cluster4 = mp_means_4)

# Join metaprograms enrichment and deconvolution
g_deconvolution <- g_deconvolution %>% 
  pivot_wider(names_from = cell_type, values_from = deconvolution_proportion)

mp_clusters_enrichment <- merge(mp_clusters_enrichment,
                                g_deconvolution,
                                by = "cell_ID")
mp_clusters_enrichment[is.na(mp_clusters_enrichment)] <- 0

# Calculate correlations between individual cell type and cluster 
lm_coefficients <- data.frame()

for(cluster in c("cluster1", "cluster2", "cluster3", "cluster4")) {
  for(cell_type in c("B_cell", "Basal", "Dendritic_cells", "Endothelial", 
                     "Epithelial", "Fibroblasts", "Luminal", "Macrophages",
                     "Mast", "pDC", "Plasma", "T_cell", "Tumor")) {
    
    cor_value <- cor(x = mp_clusters_enrichment[[cell_type]],
                     y = mp_clusters_enrichment[[cluster]],
                     method = "pearson")
    
    lm_result <- lm(mp_clusters_enrichment[[cluster]]~mp_clusters_enrichment[[cell_type]])
    sum_lm_result <- summary(lm_result)
    
    lm_coefficients <- rbind(lm_coefficients,
                             data.frame("cluster" = cluster, 
                                        "cell_type" = cell_type,
                                        "cor_value" = cor_value,
                                        "r2" = round(sum_lm_result$r.squared, 8)))
  }
}

# Save results
writexl::write_xlsx(lm_coefficients, "supplementary_table_3.xlsx")
