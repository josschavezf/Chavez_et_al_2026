library(Giotto)
library(tidyverse)
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
####           Join deconvoluted proportions            ####
#                                                          #
##%######################################################%##

source("R/join_deconvolution.R")

g_list <- list(HTA12_129, HTA12_135, HTA12_139,
               HTA12_142, HTA12_144, HTA12_145, HTA12_246, 
               HTA12_248, HTA12_249, HTA12_250, HTA12_251, 
               HTA12_253, HTA12_262, HTA12_266, HTA12_271)

names(g_list) <- c("HTA12_129", "HTA12_135", "HTA12_139", "HTA12_142",
                   "HTA12_144", "HTA12_145", "HTA12_246", "HTA12_248",
                   "HTA12_249", "HTA12_250", "HTA12_251", "HTA12_253", 
                   "HTA12_262", "HTA12_266", "HTA12_271")

deconvoluted_ct_composition <- join_deconvolution(g_list)
deconvoluted_ct_composition <- deconvoluted_ct_composition %>% 
  pivot_wider(names_from = cell_type, values_from = deconvolution_proportion)
deconvoluted_ct_composition[is.na(deconvoluted_ct_composition)] <- 0

##%######################################################%##
#                                                          #
####            Calculate nearest neighbors             ####
#                                                          #
##%######################################################%##

g <- createSpatialNetwork(g,
                          method = "kNN", 
                          k = 6,
                          name = "kNN_6")

kNN_6 <- getSpatialNetwork(g,
                           name = "kNN_6",
                           output = "networkDT")

##%######################################################%##
#                                                          #
####         Calculate neighborhood composition         ####
#                                                          #
##%######################################################%##

cell_neighborhood <- data.frame()

for (i in unique(kNN_6$from)) {
  cells_to <- kNN_6[kNN_6$from == i,]$to
  deconvoluted_ct_composition_small <- deconvoluted_ct_composition[deconvoluted_ct_composition$cell_ID %in% cells_to,-1]
  mean_deconvoluted_ct_composition_small <- colMeans(deconvoluted_ct_composition_small)
  cell_neighborhood <- rbind(cell_neighborhood, mean_deconvoluted_ct_composition_small)
}

cell_metadata <- pDataDT(g)
all_cell_ids <- cell_metadata$cell_ID
cell_ids_from <- unique(kNN_6$from)
cell_ids_missing <- all_cell_ids[!all_cell_ids %in% cell_ids_from]

cell_neighborhood_2 <- data.frame()
for (i in cell_ids_missing) {
  cells_from <- kNN_6[kNN_6$to == i,]$from
  deconvoluted_ct_composition_small <- deconvoluted_ct_composition[deconvoluted_ct_composition$cell_ID %in% cells_from,-1]
  mean_deconvoluted_ct_composition_small <- colMeans(deconvoluted_ct_composition_small)
  cell_neighborhood_2 <- rbind(cell_neighborhood_2, mean_deconvoluted_ct_composition_small)
}

rownames(cell_neighborhood) <- unique(kNN_6$from)
colnames(cell_neighborhood) <- colnames(deconvoluted_ct_composition)[-1]
rownames(cell_neighborhood_2) <- cell_ids_missing
colnames(cell_neighborhood_2) <- colnames(deconvoluted_ct_composition)[-1]

cell_neighborhood <- rbind(cell_neighborhood, cell_neighborhood_2)

##%######################################################%##
#                                                          #
####             Perform K-means clustering             ####
#                                                          #
##%######################################################%##

km_result <- kmeans(cell_neighborhood, centers = 6, nstart = 25)

#### Calculate the distance matrix
dist_matrix <- dist(cell_neighborhood)

#### Compute silhouette coefficients for each observation
sil_info <- cluster::silhouette(km_result$cluster, dist_matrix)
s_results <- as.data.frame(sil_info)

##%######################################################%##
#                                                          #
####               Spatial clusters plot                ####
#                                                          #
##%######################################################%##

spatial_locs <- getSpatialLocations(g,
                                    output = "data.table")
cell_metadata <- pDataDT(g)

s_results$cell_ID <- rownames(cell_neighborhood)
s_results <- merge(s_results, spatial_locs)
s_results <- merge(s_results, cell_metadata[, c("cell_ID", "list_ID")])

ggplot(s_results, aes(sdimx, sdimy, color = factor(cluster))) +
  geom_point(size = 0.5) +
  facet_wrap(vars(list_ID), scales = "free") +
  labs(color = "") +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank()) +
  scale_color_manual(values = c("1" = "#4287f5",
                                "2" = "#2cb816",
                                "3" = "#f5b942",
                                "4" = "#f55442",
                                "5" = "#b042f5",
                                "6" = "#f542d7")) +
  guides(color = guide_legend(override.aes = list(size = 4)))
ggsave("plots_paper/3_cellneighborhood_k6_6clusters.png", 
       width = 10, height = 9)

s_results %>% 
  filter(list_ID == "HTA12_145") %>% 
  ggplot(aes(sdimx, sdimy, color = factor(cluster))) +
  geom_point(size = 3) +
  labs(color = "") +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank()) +
  scale_color_manual(values = c("1" = "#4287f5",
                                "2" = "#2cb816",
                                "3" = "#f5b942",
                                "4" = "#f55442",
                                "5" = "#b042f5",
                                "6" = "#f542d7")) +
  guides(color = guide_legend(override.aes = list(size = 4)))
ggsave("plots_paper/3_cellneighborhood_k6_6clusters_HTA_145.png", 
       width = 8, height = 8)

s_results %>% 
  filter(list_ID == "HTA12_250") %>% 
  ggplot(aes(sdimx, sdimy, color = factor(cluster))) +
  geom_point(size = 3) +
  labs(color = "") +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank()) +
  scale_color_manual(values = c("1" = "#4287f5",
                                "2" = "#2cb816",
                                "3" = "#f5b942",
                                "4" = "#f55442",
                                "5" = "#b042f5",
                                "6" = "#f542d7")) +
  guides(color = guide_legend(override.aes = list(size = 4)))
ggsave("plots_paper/3_cellneighborhood_k6_6clusters_HTA_250.png", 
       width = 8, height = 8)

##%######################################################%##
#                                                          #
####       Heatmap with neighborhood composition        ####
#                                                          #
##%######################################################%##

m <- as.matrix(km_result$centers)
rownames(m) <- paste0("cluster", 1:6)

col_fun_prop <- colorRamp2::colorRamp2(
  c(0, 0.2, 0.4, 0.6, 0.8, 1), 
  c("lightgray", "purple", "green", "orange", "pink", "red" ))

pdf("plots_paper/3_cellneighborhoodheatmap_k6_6clusters.pdf", 
    width = 5, height = 6)
Heatmap(t(m),
        heatmap_legend_param = list(title = "",
                                    at = c(0, 0.2, 0.4, 0.6, 0.8, 1)),
        col = col_fun_prop)
dev.off()

# Barplots per sample
cluster_proportions <- s_results %>% 
  group_by(list_ID) %>% 
  count(cluster) %>% 
  mutate(total = sum(n),
         percentage = n/total)

ggplot(cluster_proportions, aes(x = percentage, y = list_ID, fill = factor(cluster))) +
  geom_col() +
  scale_fill_manual(values = c("1" = "#4287f5",
                               "2" = "#2cb816",
                               "3" = "#f5b942",
                               "4" = "#f55442",
                               "5" = "#b042f5",
                               "6" = "#f542d7")) +
  theme_bw() +
  labs(x = "",
       y = "",
       fill = "") +
  scale_x_continuous(labels = scales::label_percent()) +
  scale_y_discrete(limits = rev(c("HTA12_142", "HTA12_271", "HTA12_262", "HTA12_145",
                                  "HTA12_248", "HTA12_144", "HTA12_139", "HTA12_253", 
                                  "HTA12_135", "HTA12_249", "HTA12_266", "HTA12_250",
                                  "HTA12_129",  "HTA12_251", "HTA12_246")))
ggsave("plots_paper/3_neighborhood_barplot.pdf", device = "pdf",
       width = 5, height = 5)

# Save results
saveRDS(km_result, "km_result_k6_6clusters.RDS")
saveRDS(s_results, "s_results_k6_6clusters.RDS")
