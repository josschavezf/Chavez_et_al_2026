library(Giotto)
library(ComplexHeatmap)
source("plot_deg_heatmap.R")
source("process_sample.R")

##%######################################################%##
#                                                          #
####                    Load samples                    ####
#                                                          #
##%######################################################%##

# HTA12_129
x <- readRDS("data/scRNAseq/brca/HT206B1-XBN1-1.rds")

# HTA12_135
x_sc <- readRDS("data/scRNAseq/brca/HT243B1.rds")

# HTA12_139
x_sc <- readRDS("data/scRNAseq/brca/HT263B1.rds")

# HTA12_142
x_sc <- readRDS("data/scRNAseq/brca/HT271B1.rds")

# HTA12_144
x_sc <- readRDS("data/scRNAseq/brca/HT297B1.rds")

# HTA12_145
x_sc <- readRDS("data/scRNAseq/brca/HT305B1.rds")

# HTA12_246
x_sc <- readRDS("data/scRNAseq/brca/TWCE-HT262B1-XBN2_1.rds")

# HTA12_248
x_sc <- readRDS("data/scRNAseq/brca/HT308B1.rds")

# HTA12_249
x_sc <- readRDS("data/scRNAseq/brca/HT323B1.rds")

# HTA12_250
x_sc <- readRDS("data/scRNAseq/brca/HT339B1.rds")

# HTA12_251
x_sc <- readRDS("data/scRNAseq/brca/HT365B1.rds")

# HTA12_253
x_sc <- readRDS("data/scRNAseq/brca/snRNA_L4__HT397B1-S1H4_combo.rds")

# HTA12_262
x_sc <- readRDS("data/scRNAseq/brca/snRNA_L4__HT339B2-S1H2_combo.rds")

# HTA12_266
x_sc <- readRDS("data/scRNAseq/brca/snRNA_L4__HT425B1-S1H1_combo.rds")

# HTA12_271
x_sc <- readRDS("data/scRNAseq/brca/snRNA_L4__HT486B1-S1H1_combo.rds")

##%######################################################%##
#                                                          #
####              Convert to Giotto object              ####
#                                                          #
##%######################################################%##

g_sc <- seuratToGiottoV4(x,
                         spatial_assay = "RNA",
                         dim_reduction = "pca")

##%######################################################%##
#                                                          #
####                   Process sample                   ####
#                                                          #
##%######################################################%##

g_sc <- process_sample(g_sc)

##%######################################################%##
#                                                          #
####          Plot expression of marker genes           ####
#                                                          #
##%######################################################%##

marker_genes <- list(
    Adipocytes = c("PLIN1", "FABP4", "GLUT4", "HOXC8", "HOXC9",
                   "UCP1", "CIDEA", "PRDM16", "CD137", "TMEM26",
                   "GUCY1B1", "SLC35F1"),
    B_cell = c("CD79A", "CD79B", "CD19", "MS4A1", "IGHD", "CD22", "CD52",
               "CD20", "CD23", "CD27", "PAX5", "CD40", "CD38", "GL7",
               "RUBCNL", "COL26A1", "RGS7"),
    dendritic = c("CADM1", "XCR1", "CLEC9A", "RAB32", "C1orf54",
                  "CD1C", "FCER1A", "CLEC10A", "CD1E",
                  "LAMP3", "CCR7", "FSCN1", "CD83", "CCL22"),
    pDC = c("IL3RA", "BCL11A", "CLEC4C", "NRP1"),
    macrophage = c("CX3CR1", "CD80", "CD86", "CD163", "MSR1"),
    mast =  c("HPGD", "TPSB2", "HDC", "SLC18A2", "CPA3", "SLC8A3",
              "RGS13", "CPA3"),
    endothelial = c("EMCN", "FLT1", "PECAM1", "VWF", "PTPRB", "ACTA2", "ANGPT2"),
    fibroblast = c("COL1A1", "COL3A1", "COL5A1", "LUM", "MMP2"),
    pericyte = c("RGS5", "PLXDC1", "FN1", "MCAM"),
    plasma = c("CD38", "SDC1", "IGHG1", "IGKC", "MZB1"),
    T_cell = c("IL7R", "CD4", "CD8A", "CD8B", "CD3G", "CD3D", "CD3E",
               "FCGR3A", "GZMA", "NCAM1", "IL2RA", "CTLA4", "FOXP3", 
               "TNFRSF18", "IKZF2", "CD96", "THEMIS", "GRAP2"),
    brca_epithelial = c("LumSec", "ELF5", "CL28", "KRT15", 
                        "BARX2", "HS3ST4"),
    Luminal = c("ANKRD30A", "ERBB4", "AFF3", "TTC6", "ESR1", "NEK10", "XBP1"),
    basal = c("SAMD5", "FBXO32", "TP63", "RBBP8", "KLHL13"),
    luminal_progenitor = c("K15", "S100A7", "S100A8", "S100A9", "ALDH1A3", 
                           "CD49b", "ITGA2", "CD117", "GABRP", "CCL28")
)

my_genes <- unlist(marker_genes)
my_genes <- my_genes[my_genes %in% rownames(g_sc)]

cell_genes <- c(rep("Adipocytes", 12),
                rep("B_cell", 17),
                rep("dendritic", 14),
                rep("pDC", 4),
                rep("macrophage", 5),
                rep("mast", 8),
                rep("endothelial", 7),
                rep("fibroblast", 5),
                rep("pericyte", 4),
                rep("plasma", 5),
                rep("T_cell", 18),
                rep("brca_epithelial", 6),
                rep("Luminal", 7),
                rep("basal", 5),
                rep("luminal_progenitor", 10))

names(cell_genes) <- unlist(marker_genes)
cell_genes <- cell_genes[my_genes]

z_score_matrix <- plot_deg_heatmap(g_sc,
                                   features = my_genes,
                                   cluster_column = "leiden_clus",
                                   values = "normalized",
                                   return_zscores_matrix = TRUE)

Heatmap(z_score_matrix,
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        row_split = cell_genes,
        row_title_rot = 0,
        row_order = my_genes,
        row_names_gp = gpar(fontsize = 3),
        row_title_gp = gpar(fontsize = 10),
        heatmap_legend_param = list(title = ""))

##%######################################################%##
#                                                          #
####        Add cell type to each Leiden cluster        ####
#                                                          #
##%######################################################%##

# HTA12_129
cell_types <- c("1" = "T_cell",
                "2" = "T_cell",
                "3" = "Macrophages",
                "4" = "T_cell",
                "5" = "Mast",
                "6" = "Macrophages",
                "7" = "T_cell",
                "8" = "Fibroblasts",
                "9" = "Epithelial",
                "10" = "Epithelial",
                "11" = "B_cell",
                "12" = "pDC",
                "13" = "Plasma",
                "14" = "Epithelial",
                "15" = "Luminal",
                "16" = "Dendritic_cells",
                "17" = "Basal",
                "18" = "Endothelial")

# HTA12_135
cell_types <- c("1" = "Macrophages",
                "2" = "Luminal",
                "3" = "Fibroblasts",
                "4" = "Luminal",
                "5" = "Luminal",
                "6" = "Luminal",
                "7" = "Epithelial",
                "8" = "T_cell",
                "9.1" = "Plasma",
                "9.2" = "B_cell",
                "9.3" = "B_cell",
                "9.4" = "Plasma",
                "10.1" = "Luminal",
                "10.2" = "Luminal",
                "10.3" = "Luminal",
                "10.4" = "Mast",
                "10.5" = "Luminal",
                "11" = "Luminal",
                "12" = "Basal",
                "13" = "Luminal",
                "14" = "Macrophages",
                "15" = "Dendritic_cells",
                "16" = "T_cell",
                "17" = "pDC")

# HTA12_139
cell_types <- c("1" = "Luminal",
                "2" = "Luminal",
                "3" = "Fibroblasts",
                "4" = "Luminal",
                "5" = "Luminal",
                "6" = "Luminal",
                "7" = "Luminal",
                "8" = "Macrophages",
                "9" = "Luminal",
                "10" = "Luminal",
                "11" = "Basal",
                "12" = "Luminal",
                "13" = "Epithelial",
                "14" = "Endothelial",
                "15" = "Adipocytes")

# HTA12_142
cell_types <- c("1" = "Tumor",
                "2" = "Tumor",
                "3" = "Tumor",
                "4" = "Tumor",
                "5" = "Tumor",
                "6" = "Tumor",
                "7" = "Tumor",
                "8" = "Tumor",
                "9.1" = "T_cell",
                "9.2" = "Macrophages",
                "9.3" = "Fibroblasts",
                "9.4" = "T_cell",
                "9.5" = "Plasma",
                "9.6" = "Luminal",
                "10" = "Tumor")

# HTA12_144
cell_types <- c("1" = "Luminal",
                "2" = "Luminal",
                "3" = "Macrophages",
                "4" = "Luminal",
                "5" = "Fibroblasts",
                "6" = "Luminal",
                "7" = "Luminal",
                "8" = "T_cell",
                "9" = "Luminal_progenitor",
                "10" = "Basal",
                "11" = "Endothelial",
                "12" = "Plasma")

# HTA12_145
cell_types <- c("1" = "Luminal",
                "2" = "Luminal",
                "3" = "Luminal",
                "4" = "Luminal",
                "5" = "Luminal",
                "6" = "Luminal",
                "7" = "Luminal",
                "8.1" = "Fibroblasts",
                "8.2" = "Fibroblasts",
                "8.3" = "Fibroblasts",
                "8.4" = "Fibroblasts",
                "8.5" = "Adipocytes",
                "8.6" = "Fibroblasts",
                "8.7" = "Fibroblasts",
                "9.1" = "Macrophages",
                "9.2" = "Macrophages",
                "9.3" = "Macrophages",
                "9.4" = "Macrophages",
                "9.5" = "T_cell",
                "9.6" = "Macrophages",
                "9.7" = "T_cell",
                "10" = "Luminal",
                "11" = "Luminal_progenitor",
                "12" = "Basal",
                "13" = "Luminal",
                "14" = "Endothelial")

# HTA12_246
cell_types <- c("1" = "Fibroblasts",
                "2" = "Fibroblasts",
                "3" = "Basal",
                "4" = "Endothelial",
                "5" = "Luminal",
                "6" = "Basal",
                "7" = "Luminal",
                "8" = "Luminal",
                "9" = "Fibroblasts",
                "10" = "Luminal",
                "11" = "T_cell",
                "12" = "Macrophages",
                "13" = "Mast",
                "14" = "Adipocytes",
                "15" = "Luminal_progenitor")

# HTA12_248
cell_types <- c("1" = "Luminal",
                "2" = "Luminal",
                "3" = "Luminal",
                "4" = "Luminal",
                "5" = "Luminal",
                "6" = "Macrophages",
                "7" = "Luminal",
                "8" = "T_cell",
                "9" = "T_cell",
                "10" = "Luminal",
                "11" = "Luminal",
                "12" = "Fibroblasts",
                "13" = "Plasma",
                "14" = "Luminal",
                "15" = "pDC")

# HTA12_249
cell_types <- c("1" = "T_cell",
                "2" = "Fibroblasts",
                "3" = "Luminal",
                "4" = "Macrophages",
                "5" = "Luminal",
                "6" = "Macrophages",
                "7" = "Luminal",
                "8" = "Luminal",
                "9" = "Luminal",
                "10" = "Basal",
                "11" = "Luminal",
                "12" = "B_cell",
                "13" = "Mast",
                "14" = "Luminal_progenitor",
                "15" = "Endothelial")

# HTA12_250
cell_types <- c("1" = "Macrophages",
                "2" = "Epithelial",
                "3" = "Fibroblasts",
                "4" = "T_cell",
                "5" = "Epithelial",
                "6" = "Luminal",
                "7" = "Basal",
                "8" = "T_cell",
                "9" = "Macrophages",
                "10" = "Epithelial",
                "11" = "Epithelial",
                "12" = "pDC",
                "13" = "Epithelial",
                "14" = "T_cell",
                "15" = "Fibroblasts",
                "16" = "Mast")

# HTA12_251
cell_types <- c("1" = "T_cell",
                "2" = "T_cell",
                "3" = "Macrophages",
                "4" = "Macrophages",
                "5" = "Luminal",
                "6" = "Luminal",
                "7" = "T_cell",
                "8" = "Macrophages",
                "9" = "Macrophages",
                "10" = "T_cell",
                "11" = "Fibroblasts",
                "12" = "Fibroblasts",
                "13" = "Plasma",
                "14" = "Macrophages",
                "15" = "B_cell",
                "16" = "Luminal",
                "17" = "Luminal",
                "18" = "pDC")

# HTA12_253
cell_types <- c("1" = "Macrophages",
                "2" = "Luminal",
                "3" = "T_cell",
                "4" = "Luminal",
                "5" = "Luminal",
                "6" = "Luminal",
                "7.1" = "Fibroblasts",
                "7.2" = "Fibroblasts",
                "7.3" = "Fibroblasts",
                "7.4" = "Fibroblasts",
                "7.5" = "Fibroblasts",
                "7.6" = "Fibroblasts",
                "7.7" = "Fibroblasts",
                "7.8" = "Fibroblasts",
                "8" = "Luminal",
                "9" = "Basal",
                "10.1" = "Plasma",
                "10.2" = "Plasma",
                "10.3" = "Plasma",
                "10.4" = "Plasma",
                "10.5" = "Plasma",
                "10.6" = "Plasma",
                "11" = "B_cell",
                "12" = "pDC",
                "13" = "Endothelial",
                "14" = "Endothelial")

# HTA12_262
cell_types <- c("1" = "Epithelial",
                "2" = "Epithelial",
                "3" = "Epithelial",
                "4" = "Epithelial",
                "5" = "Epithelial",
                "6" = "Epithelial",
                "7" = "Epithelial",
                "8" = "Epithelial",
                "9" = "Epithelial",
                "10.1" = "Fibroblasts",
                "10.2" = "Pericytes",
                "10.3" = "Pericytes",
                "10.4" = "Fibroblasts",
                "10.5" = "Fibroblasts",
                "10.6" = "Pericytes",
                "11" = "Epithelial",
                "12" = "Epithelial",
                "13" = "Macrophages",
                "14" = "Endothelial",
                "15" = "Basal")

# HTA12_266
cell_types <- c("1" = "Macrophages",
                "2" = "Luminal",
                "3" = "Fibroblasts",
                "4" = "Luminal",
                "5" = "Luminal",
                "6" = "Fibroblasts",
                "7" = "Macrophages",
                "8" = "Fibroblasts",
                "9" = "Luminal",
                "10" = "T_cell",
                "11" = "Fibroblasts",
                "12" = "T_cell",
                "13" = "Endothelial",
                "14" = "Plasma",
                "15" = "Macrophages")

# HTA12_271
cell_types <- c("1" = "Fibroblasts",
                "2" = "Fibroblasts",
                "3" = "Fibroblasts",
                "4" = "Fibroblasts",
                "5" = "Fibroblasts",
                "6" = "Fibroblasts",
                "7" = "Fibroblasts",
                "8.1" = "Macrophages",
                "8.2" = "Macrophages",
                "8.3" = "Macrophages",
                "8.4" = "Macrophages",
                "9" = "Fibroblasts",
                "10.1" = "Pericytes",
                "10.2" = "Pericytes",
                "10.3" = "Pericytes",
                "10.4" = "Pericytes",
                "10.5" = "Pericytes",
                "10.6" = "Pericytes",
                "11" = "Fibroblasts",
                "12" = "Endothelial",
                "13" = "Fibroblasts")


##%######################################################%##
#                                                          #
####              Add annotation to object              ####
#                                                          #
##%######################################################%##

g_sc <- annotateGiotto(g_sc,
                       annotation_vector = cell_types,
                       cluster_column = "leiden_clus")

##%######################################################%##
#                                                          #
####          Label-transferring of tumor cells         ####
#                                                          #
##%######################################################%##

cell_metadata <- pDataDT(g_sc)
cell_metadata$cell_types_final <- cell_metadata$cell_types
cell_metadata$cell_types_final[cell_metadata$cell_type_broad == "Tumor"] <- "Tumor"

cell_metadata <- cell_metadata[, c("cell_ID", "cell_types", "cell_types_final")]

g_sc <- addCellMetadata(g_sc,
                        new_metadata = cell_metadata,
                        by_column = TRUE)

##%######################################################%##
#                                                          #
####                     Plot UMAP                      ####
#                                                          #
##%######################################################%##

plotUMAP(g_sc,
         point_size = 0.5,
         cell_color = "cell_types_final",
         cell_color_code = c("Adipocytes" = "#d6d62d",
                             "B_cell" = "#7ba7ed", 
                             "Basal" = "#1168f5", 
                             "Dendritic_cells" = "#02239e",
                             "Endothelial" = "#edad87",
                             "Epithelial" = "#9335b8",
                             "Fibroblasts" = "#f57842",
                             "Luminal" = "#5e1c58",
                             "Luminal_progenitor" = "#fc0303",
                             "Macrophages" = "#fae987",
                             "Mast" = "#f272a9",
                             "pDC" = "#87faf8", 
                             "Pericytes" = "#e63cc9",
                             "Plasma" = "#91030f",
                             "T_cell" = "#169103",
                             "Tumor" = "#7c7d80"
         ))

##%######################################################%##
#                                                          #
####               Save annotated object                ####
#                                                          #
##%######################################################%##

saveRDS(g_sc, "data/annotated_scRNAseq/brca/HTA12_129.RDS")
saveRDS(g_sc, "data/annotated_scRNAseq/brca/HTA12_135.RDS")
saveRDS(g_sc, "data/annotated_scRNAseq/brca/HTA12_139.RDS")
saveRDS(g_sc, "data/annotated_scRNAseq/brca/HTA12_142.RDS")
saveRDS(g_sc, "data/annotated_scRNAseq/brca/HTA12_144.RDS")
saveRDS(g_sc, "data/annotated_scRNAseq/brca/HTA12_145.RDS")
saveRDS(g_sc, "data/annotated_scRNAseq/brca/HTA12_246.RDS")
saveRDS(g_sc, "data/annotated_scRNAseq/brca/HTA12_248.RDS")
saveRDS(g_sc, "data/annotated_scRNAseq/brca/HTA12_249.RDS")
saveRDS(g_sc, "data/annotated_scRNAseq/brca/HTA12_250.RDS")
saveRDS(g_sc, "data/annotated_scRNAseq/brca/HTA12_251.RDS")
saveRDS(g_sc, "data/annotated_scRNAseq/brca/HTA12_253.RDS")
saveRDS(g_sc, "data/annotated_scRNAseq/brca/HTA12_262.RDS")
saveRDS(g_sc, "data/annotated_scRNAseq/brca/HTA12_266.RDS")
saveRDS(g_sc, "data/annotated_scRNAseq/brca/HTA12_271.RDS")
