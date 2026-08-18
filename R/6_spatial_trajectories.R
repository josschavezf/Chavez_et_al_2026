library(Giotto)
library(tidyr)
library(ggplot2)
library(ComplexHeatmap)
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
####               Create metadata file                 ####
#                                                          #
##%######################################################%##

cell_metadata <- pDataDT(g)
spat_locs <- getSpatialLocations(g,
                                 output = "data.table")

samples_er <- c("HTA12_135", "HTA12_139", "HTA12_144", "HTA12_145", "HTA12_246",
                "HTA12_248", "HTA12_249", "HTA12_253", "HTA12_266")
samples_tnbc <- c("HTA12_129", "HTA12_142", "HTA12_250",
                  "HTA12_262", "HTA12_271")
samples_her2 <- "HTA12_251"

cell_ids_er <- cell_metadata[cell_metadata$list_ID %in% samples_er,]$cell_ID
cell_ids_tnbc <- cell_metadata[cell_metadata$list_ID %in% samples_tnbc,]$cell_ID
cell_ids_her2 <- cell_metadata[cell_metadata$list_ID %in% samples_her2,]$cell_ID

cell_metadata <- cell_metadata[cell_metadata$cell_ID %in% cell_ids_er,]
spat_locs <- spat_locs[spat_locs$cell_ID %in% cell_ids_er,]

cell_metadata <- cell_metadata[cell_metadata$cell_ID %in% cell_ids_tnbc,]
spat_locs <- spat_locs[spat_locs$cell_ID %in% cell_ids_tnbc,]

cell_metadata <- cell_metadata[cell_metadata$cell_ID %in% cell_ids_her2,]
spat_locs <- spat_locs[spat_locs$cell_ID %in% cell_ids_her2,]

meta_input <- data.frame(Spot_ID = cell_metadata$cell_ID,
                         x = spat_locs$sdimx,
                         y = spat_locs$sdimy,
                         Sample = cell_metadata$list_ID)

dir.create("ontrac_files/input_files/brca_er")
write.csv(meta_input, "ontrac_files/input_files/brca_er/meta_input.csv",
          quote = FALSE, row.names = FALSE)

dir.create("ontrac_files/input_files/brca_tnbc")
write.csv(meta_input, "ontrac_files/input_files/brca_tnbc/meta_input.csv",
          quote = FALSE, row.names = FALSE)

dir.create("ontrac_files/input_files/brca_her2")
write.csv(meta_input, "ontrac_files/input_files/brca_her2/meta_input.csv",
          quote = FALSE, row.names = FALSE)

##%######################################################%##
#                                                          #
####             Join cell type composition             ####
#                                                          #
##%######################################################%##

source("R/join_deconvolution.R")

### ER
g_list <- list(HTA12_135, HTA12_139, HTA12_144, HTA12_145, HTA12_246,
               HTA12_248, HTA12_249, HTA12_253, HTA12_266)
names(g_list) <- c("HTA12_135", "HTA12_139", "HTA12_144", "HTA12_145", 
                   "HTA12_246", "HTA12_248", "HTA12_249", "HTA12_253", 
                   "HTA12_266")

### TNBC
g_list <- list(HTA12_129, HTA12_142, HTA12_250,
               HTA12_262, HTA12_271)
names(g_list) <- c("HTA12_129", "HTA12_142", "HTA12_250", 
                   "HTA12_262", "HTA12_271")

### HER2
g_list <- list(HTA12_251)
names(g_list) <- "HTA12_251"

deconvoluted_ct_composition <- join_deconvolution(g_list)

## ER+
deconvoluted_ct_composition <- deconvoluted_ct_composition %>% 
  pivot_wider(names_from = cell_type, values_from = deconvolution_proportion) %>% 
  relocate(cell_ID, Adipocytes, B_cell, Basal, Dendritic_cells, Endothelial,
           Epithelial, Fibroblasts, Luminal, Luminal_progenitor, Macrophages, 
           Mast, pDC, Plasma, T_cell, Tumor)

## TNBC
deconvoluted_ct_composition <- deconvoluted_ct_composition %>% 
  pivot_wider(names_from = cell_type, values_from = deconvolution_proportion) %>% 
  relocate(cell_ID, B_cell, Basal, Dendritic_cells, Endothelial,
           Epithelial, Fibroblasts, Luminal, Macrophages, 
           Mast, pDC, Pericytes, Plasma, T_cell, Tumor)

## HER2
deconvoluted_ct_composition <- deconvoluted_ct_composition %>% 
  pivot_wider(names_from = cell_type, values_from = deconvolution_proportion) %>% 
  relocate(cell_ID, B_cell, Fibroblasts, Luminal, Macrophages, pDC, Plasma, 
           T_cell, Tumor)

deconvoluted_ct_composition[is.na(deconvoluted_ct_composition)] <- 0

write.csv(deconvoluted_ct_composition, "ontrac_files/input_files/brca_er/ct_composition.csv",
          row.names = FALSE, quote = FALSE)
write.csv(deconvoluted_ct_composition, "ontrac_files/input_files/brca_tnbc/ct_composition.csv",
          row.names = FALSE, quote = FALSE)
write.csv(deconvoluted_ct_composition, "ontrac_files/input_files/brca_her2/ct_composition.csv",
          row.names = FALSE, quote = FALSE)

##%######################################################%##
#                                                          #
####            Join deconvoluted expression            ####
#                                                          #
##%######################################################%##

HTA12_129 <- readRDS("data/annotated_scRNAseq/brca/HTA12_129.RDS")
HTA12_135 <- readRDS("data/annotated_scRNAseq/brca/HTA12_135.RDS")
HTA12_139 <- readRDS("data/annotated_scRNAseq/brca/HTA12_139.RDS")
HTA12_142 <- readRDS("data/annotated_scRNAseq/brca/HTA12_142.RDS")
HTA12_144 <- readRDS("data/annotated_scRNAseq/brca/HTA12_144.RDS")
HTA12_145 <- readRDS("data/annotated_scRNAseq/brca/HTA12_145.RDS")
HTA12_246 <- readRDS("data/annotated_scRNAseq/brca/HTA12_246.RDS")
HTA12_248 <- readRDS("data/annotated_scRNAseq/brca/HTA12_248.RDS")
HTA12_249 <- readRDS("data/annotated_scRNAseq/brca/HTA12_249.RDS")
HTA12_250 <- readRDS("data/annotated_scRNAseq/brca/HTA12_250.RDS")
HTA12_251 <- readRDS("data/annotated_scRNAseq/brca/HTA12_251.RDS")
HTA12_253 <- readRDS("data/annotated_scRNAseq/brca/HTA12_253.RDS")
HTA12_262 <- readRDS("data/annotated_scRNAseq/brca/HTA12_262.RDS")
HTA12_266 <- readRDS("data/annotated_scRNAseq/brca/HTA12_266.RDS")
HTA12_271 <- readRDS("data/annotated_scRNAseq/brca/HTA12_271.RDS")

### ER
g_sc <- joinGiottoObjects(
  gobject_list = list(HTA12_135, HTA12_139, HTA12_144, HTA12_145, 
                      HTA12_246, HTA12_248, HTA12_249, HTA12_253,
                      HTA12_266),
  gobject_names = c("HTA12_135", "HTA12_139", "HTA12_144", "HTA12_145", 
                    "HTA12_246", "HTA12_248", "HTA12_249", "HTA12_253", 
                    "HTA12_266"))

### TNBC
g_sc <- joinGiottoObjects(
  gobject_list = list(HTA12_129, HTA12_142, HTA12_250, 
                      HTA12_262, HTA12_271),
  gobject_names = c("HTA12_129", "HTA12_142", "HTA12_250", 
                    "HTA12_262", "HTA12_271"))

### HER2
g_sc <- HTA12_251

## Prepare data
g_sc <- normalizeGiotto(g_sc)

cell_metadata <- pDataDT(g_sc)

df <- data.frame()
m <- getExpression(g_sc,
                   values = "normalized",
                   output = "matrix")
m <- as.matrix(m)

for (i in sort(unique(cell_metadata$cell_types_final))) {
  cell_ids <- cell_metadata[cell_metadata$cell_types_final == i,]$cell_ID
  x_small <- m[, cell_ids]
  mean_expression <- rowMeans(x_small)
  df <- rbind(df, mean_expression)
}

colnames(df) <- rownames(m)
df <- cbind(Cell_Type = sort(unique(cell_metadata$cell_types_final)),
            df)

write.csv(df, "ontrac_files/input_files/brca_er/deconvoluted_expression.csv",
          row.names = FALSE, quote = FALSE)
write.csv(df, "ontrac_files/input_files/brca_tnbc/deconvoluted_expression.csv",
          row.names = FALSE, quote = FALSE)
write.csv(df, "ontrac_files/input_files/brca_her2/deconvoluted_expression.csv",
          row.names = FALSE, quote = FALSE)

##%######################################################%##
#                                                          #
####                     Run ONTraC                     ####
#                                                          #
##%######################################################%##

## Run ONTraC on a bash terminal, ER
ONTraC --meta-input ontrac_files/input_files/brca_er/meta_input.csv \
--NN-dir ontrac_files/results/brca_er/nn_dir/ \
--GNN-dir ontrac_files/results/brca_er/gnn_dir/ \
--NT-dir ontrac_files/results/brca_er/nt_dir/ \
--deconvoluted-ct-composition ontrac_files/input_files/brca_er/ct_composition.csv \
--deconvoluted-exp-input ontrac_files/input_files/brca_er/deconvoluted_expression.csv \
-k 6 \
--n-neighbors 6 \
--regularization-loss-weight 8 \
--modularity-loss-weight 1 \
-s 1991 \
--sigma 1 \
--embedding-adjust

ONTraC_analysis --meta-input ontrac_files/input_files/brca_er/meta_input.csv \
--NN-dir ontrac_files/results/brca_er/nn_dir/ \
--GNN-dir ontrac_files/results/brca_er/gnn_dir/ \
--NT-dir ontrac_files/results/brca_er/nt_dir/ \
-o ontrac_files/results/brca_er/plots \
--embedding-adjust \
--suppress-niche-cluster-loadings \
--suppress-cell-type-composition \
--sigma 1

## Run ONTraC on a bash terminal, TNBC
ONTraC --meta-input ontrac_files/input_files/brca_tnbc/meta_input.csv \
--NN-dir ontrac_files/results/brca_tnbc/nn_dir/ \
--GNN-dir ontrac_files/results/brca_tnbc/gnn_dir/ \
--NT-dir ontrac_files/results/brca_tnbc/nt_dir/ \
--deconvoluted-ct-composition ontrac_files/input_files/brca_tnbc/ct_composition.csv \
--deconvoluted-exp-input ontrac_files/input_files/brca_tnbc/deconvoluted_expression.csv \
-k 6 \
--n-neighbors 6 \
--regularization-loss-weight 8 \
--modularity-loss-weight 1 \
-s 1991 \
--sigma 1 \
--embedding-adjust

ONTraC_analysis --meta-input ontrac_files/input_files/brca_tnbc/meta_input.csv \
--NN-dir ontrac_files/results/brca_tnbc/nn_dir/ \
--GNN-dir ontrac_files/results/brca_tnbc/gnn_dir/ \
--NT-dir ontrac_files/results/brca_tnbc/nt_dir/ \
-o ontrac_files/results/brca_tnbc/plots \
--embedding-adjust \
--suppress-niche-cluster-loadings \
--suppress-cell-type-composition \
--sigma 1

## Run ONTraC on a bash terminal, HER2
ONTraC --meta-input ontrac_files/input_files/brca_her2/meta_input.csv \
--NN-dir ontrac_files/results/brca_her2/nn_dir/ \
--GNN-dir ontrac_files/results/brca_her2/gnn_dir/ \
--NT-dir ontrac_files/results/brca_her2/nt_dir/ \
--deconvoluted-ct-composition ontrac_files/input_files/brca_her2/ct_composition.csv \
--deconvoluted-exp-input ontrac_files/input_files/brca_her2/deconvoluted_expression.csv \
-k 6 \
--n-neighbors 6 \
--regularization-loss-weight 8 \
--modularity-loss-weight 1 \
-s 1991 \
--sigma 1 \
--embedding-adjust

ONTraC_analysis --meta-input ontrac_files/input_files/brca_her2/meta_input.csv \
--NN-dir ontrac_files/results/brca_her2/nn_dir/ \
--GNN-dir ontrac_files/results/brca_her2/gnn_dir/ \
--NT-dir ontrac_files/results/brca_her2/nt_dir/ \
-o ontrac_files/results/brca_her2/plots \
--embedding-adjust \
--suppress-niche-cluster-loadings \
--suppress-cell-type-composition \
-r \
--sigma 1

##%######################################################%##
#                                                          #
####            Create spatial NT score plot            ####
#                                                          #
##%######################################################%##

# Run this in python using the ONTraC conda environment
import os 
import sys
import matplotlib as mpl
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'

from ONTraC. utils import write_version_info
write_version_info()

from ONTraC.analysis.data import AnaData
from ONTraC.analysis.cell_type import plot_violin_cell_type_composition_along_NT_score_from_anadata
from optparse import Values

options = Values()

options.NN_dir = 'ontrac_files/results/brca_her2/nn_dir'  
options.GNN_dir = 'ontrac_files/results/brca_her2/gnn_dir'
options.NT_dir = 'ontrac_files/results/brca_her2/nt_dir'

options. embedding_adjust = True
options.sigma = 1
options.reverse = True  
options.output = None  
ana_data = AnaData(options)

NT_score=ana_data.NT_score
meta_data_df=ana_data.meta_data_df
reverse=ana_data.options.reverse
output_file_path=ana_data.options.output

samples = meta_data_df['Sample'].unique().tolist()

N = len(samples)
fig, axes = plt.subplots(1, N, figsize=(3.5 * N, 3))  # control figsize here
for i, sample in enumerate(samples):
  sample_df = NT_score.loc[meta_data_df[meta_data_df['Sample'] == sample].index]
ax: plt.Axes = axes[i] if N > 1 else axes  # type: ignore
NT_score_values = sample_df['Cell_NTScore'] if not reverse else 1 - sample_df['Cell_NTScore']
scatter = ax.scatter(sample_df['x'], sample_df['y'], c=NT_score_values, cmap='rainbow', vmin=0, vmax=1, s=1) # contorl the points size using `s`
ax.set_xticks([])
ax.set_yticks([])
plt.colorbar(scatter)
ax.set_title(f"{sample} Cell-level NT Score")

fig.tight_layout()
fig.savefig('ontrac_files/results/brca_tnbc/plots/cellNTscore.pdf', transparent=True)

##%######################################################%##
#                                                          #
####                 Create violin plot                 ####
#                                                          #
##%######################################################%##

plt.rcParams['ytick.labelsize'] = 6

## ER
fig, axes = plot_violin_cell_type_composition_along_NT_score_from_anadata(
  ana_data, 
  max_violin_width = 0.6, 
  palette = {"Adipocytes": "#d6d62d", "B_cell": "#7ba7ed", "Basal": "#1168f5", 
    "Dendritic_cells": "#02239e", "Endothelial": "#edad87", "Epithelial": "#9335b8", 
    "Fibroblasts": "#f57842", "Luminal": "#5e1c58", "Luminal_progenitor": "#fc0303", 
    "Macrophages": "#fae987", "Mast": "#f272a9", "pDC": "#87faf8", 
    "Pericytes": "#e63cc9", "Plasma": "#91030f", "T_cell": "#169103", "Tumor": "#7c7d80"}, 
  figsize = [9, 6], 
  show_fill = True,
  violin_width_scale = {"Adipocytes": 20, "B_cell": 30, "Basal": 10, 
    "Dendritic_cells": 100, "Endothelial": 5, "Epithelial": 5, 
    "Fibroblasts": 2, "Luminal": 10, "Luminal_progenitor": 30, "Macrophages": 4, 
    "Mast": 40, "pDC": 20, "Pericytes": 10, "Plasma": 5, "T_cell": 10, 
    "Tumor": 1})
fig.savefig('ontrac_files/results/brca_er/plots/violin.pdf', transparent=True)

## TNBC
fig, axes = plot_violin_cell_type_composition_along_NT_score_from_anadata(
  ana_data, max_violin_width=0.5, 
  palette = {"Adipocytes": "#d6d62d", "B_cell": "#7ba7ed", "Basal": "#1168f5", 
    "Dendritic_cells": "#02239e", "Endothelial": "#edad87", "Epithelial": "#9335b8", 
    "Fibroblasts": "#f57842", "Luminal": "#5e1c58", "Luminal_progenitor": "#fc0303", 
    "Macrophages": "#fae987", "Mast": "#f272a9", "pDC": "#87faf8", 
    "Pericytes": "#e63cc9", "Plasma": "#91030f", "T_cell": "#169103", "Tumor": "#7c7d80"}, 
  figsize = [9,6], 
  show_fill = True,
  violin_width_scale = {"Adipocytes": 20, "B_cell": 30, "Basal": 10, 
    "Dendritic_cells": 100, "Endothelial": 5, "Epithelial": 5, 
    "Fibroblasts": 2, "Luminal": 10, "Luminal_progenitor": 30, "Macrophages": 4, 
    "Mast": 40, "pDC": 20, "Pericytes": 10, "Plasma": 5, "T_cell": 10, 
    "Tumor": 1})
fig.savefig('ontrac_files/results/brca_tnbc/plots/violin.pdf', transparent=True)

## HER2
fig, axes = plot_violin_cell_type_composition_along_NT_score_from_anadata(
  ana_data, max_violin_width=0.4, 
  palette = {"Adipocytes": "#d6d62d", "B_cell": "#7ba7ed", "Basal": "#1168f5", 
    "Dendritic_cells": "#02239e", "Endothelial": "#edad87", "Epithelial": "#9335b8", 
    "Fibroblasts": "#f57842", "Luminal": "#5e1c58", "Luminal_progenitor": "#fc0303", 
    "Macrophages": "#fae987", "Mast": "#f272a9", "pDC": "#87faf8", 
    "Pericytes": "#e63cc9", "Plasma": "#91030f", "T_cell": "#169103", "Tumor": "#7c7d80"}, 
  figsize = [9,6], 
  show_fill = True,
  violin_width_scale = {"Adipocytes": 20, "B_cell": 30, "Basal": 10, 
    "Dendritic_cells": 100, "Endothelial": 5, "Epithelial": 5, 
    "Fibroblasts": 2, "Luminal": 10, "Luminal_progenitor": 30, "Macrophages": 4, 
    "Mast": 40, "pDC": 20, "Pericytes": 10, "Plasma": 5, "T_cell": 10, 
    "Tumor": 1})
fig.savefig('ontrac_files/results/brca_her2/plots/violin.pdf', transparent=True)

##%######################################################%##
#                                                          #
####          Cell-type composition per niche           ####
#                                                          #
##%######################################################%##

cell_type_enrichment <- getSpatialEnrichment(g,
                                             name = "PAGE_marker_genes",
                                             output = "data.table")

# ER
niches <- read.csv("ontrac_files/results/brca_er/gnn_dir/cell_level_max_niche_cluster.csv.gz")
colnames(niches)[1] <- "cell_ID"
niches <- merge(niches, cell_type_enrichment)

niche_composition <- data.frame()
for(i in 0:5) {
  niches_small <- niches[niches$Niche_Cluster == i,-1]
  mean_composition <- colMeans(niches_small)
  niche_composition <- rbind(niche_composition, mean_composition)
}
colnames(niche_composition) <- colnames(niches)[-1]
rownames(niche_composition) <- paste0("Niche_Cluster_", niche_composition$Niche_Cluster)
niche_composition <- niche_composition[,-1]
niche_composition <- as.matrix(niche_composition)

pdf("plots_paper/6_niche_composition_er.pdf", width = 5, height = 5)
Heatmap(t(niche_composition),
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        column_order = c("Niche_Cluster_0",
                         "Niche_Cluster_2",
                         "Niche_Cluster_1",
                         "Niche_Cluster_4",
                         "Niche_Cluster_3",
                         "Niche_Cluster_5"),
        heatmap_legend_param = list(title = ""))
dev.off()

# TNBC
niches <- read.csv("ontrac_files/results/brca_tnbc/gnn_dir/cell_level_max_niche_cluster.csv.gz")
colnames(niches)[1] <- "cell_ID"
niches <- merge(niches, cell_type_enrichment)

niche_composition <- data.frame()
for(i in 0:5) {
  niches_small <- niches[niches$Niche_Cluster == i,-1]
  mean_composition <- colMeans(niches_small)
  niche_composition <- rbind(niche_composition, mean_composition)
}
colnames(niche_composition) <- colnames(niches)[-1]
rownames(niche_composition) <- paste0("Niche_Cluster_", niche_composition$Niche_Cluster)
niche_composition <- niche_composition[,-1]
niche_composition <- as.matrix(niche_composition)

pdf("plots_paper/6_niche_composition_tnbc.pdf", width = 5, height = 5)
Heatmap(t(niche_composition),
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        column_order = c("Niche_Cluster_1",
                         "Niche_Cluster_0",
                         "Niche_Cluster_4",
                         "Niche_Cluster_5",
                         "Niche_Cluster_3",
                         "Niche_Cluster_2"),
        heatmap_legend_param = list(title = ""))
dev.off()

# HER2
niches <- read.csv("ontrac_files/results/brca_her2/gnn_dir/cell_level_max_niche_cluster.csv.gz")
colnames(niches)[1] <- "cell_ID"
niches <- merge(niches, cell_type_enrichment)

niche_composition <- data.frame()
for(i in 0:5) {
  niches_small <- niches[niches$Niche_Cluster == i,-1]
  mean_composition <- colMeans(niches_small)
  niche_composition <- rbind(niche_composition, mean_composition)
}
colnames(niche_composition) <- colnames(niches)[-1]
rownames(niche_composition) <- paste0("Niche_Cluster_", niche_composition$Niche_Cluster)
niche_composition <- niche_composition[,-1]
niche_composition <- as.matrix(niche_composition)

pdf("plots_paper/6_niche_composition_her2.pdf", width = 5, height = 5)
Heatmap(t(niche_composition),
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        column_order = c("Niche_Cluster_1",
                         "Niche_Cluster_2",
                         "Niche_Cluster_4",
                         "Niche_Cluster_5",
                         "Niche_Cluster_3",
                         "Niche_Cluster_0"),
        heatmap_legend_param = list(title = ""))
dev.off()
