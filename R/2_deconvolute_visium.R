library(Giotto)
library(spacexr)
library(SummarizedExperiment)
library(SpatialExperiment)
source("R/process_st_sample.R")

##%######################################################%##
#                                                          #
####                     HTA12_129                      ####
#                                                          #
##%######################################################%##

# Read Visium sample
x <- readRDS("data/visium/brca/HT206B1-S1Fc1U2Z1B1-SeuratObj.rds")
g <- seuratToGiottoV4(x)

# Read annotated scRNAseq
g_sc <- readRDS("data/annotated_scRNAseq/brca/HTA12_129.RDS")

##%######################################################%##
#                                                          #
####                     HTA12_135                      ####
#                                                          #
##%######################################################%##

# Read Visium sample
x <- readRDS("data/visium/brca/HT243B1H3A2-S1Fc1U4Z1B1-SeuratObj.rds")
g <- seuratToGiottoV4(x)

# Read annotated scRNAseq
g_sc <- readRDS("data/annotated_scRNAseq/brca/HTA12_135.RDS")

##%######################################################%##
#                                                          #
####                     HTA12_139                      ####
#                                                          #
##%######################################################%##

# Read Visium sample
x <- readRDS("data/visium/brca/HT263B1-S1H1Fc2U1Z1Bs1-SeuratObj.rds")
g <- seuratToGiottoV4(x)

# Read annotated scRNAseq
g_sc <- readRDS("data/annotated_scRNAseq/brca/HTA12_139.RDS")

##%######################################################%##
#                                                          #
####                     HTA12_142                      ####
#                                                          #
##%######################################################%##

# Read Visium sample
x <- readRDS("data/visium/brca/HT271B1-S1H3Fc2U1Z1Bs1-SeuratObj.rds")
g <- seuratToGiottoV4(x)

# Read annotated scRNAseq
g_sc <- readRDS("data/annotated_scRNAseq/brca/HTA12_142.RDS")

##%######################################################%##
#                                                          #
####                     HTA12_144                      ####
#                                                          #
##%######################################################%##

# Read Visium sample
x <- readRDS("data/visium/brca/HT297B1-S1H1Fc2U1Z1Bs1-SeuratObj.rds")
g <- seuratToGiottoV4(x)

# Read annotated scRNAseq
g_sc <- readRDS("data/annotated_scRNAseq/brca/HTA12_144.RDS")

##%######################################################%##
#                                                          #
####                     HTA12_145                      ####
#                                                          #
##%######################################################%##

# Read Visium sample
x <- readRDS("data/visium/brca/HT305B1-S1H1Fc2U1Z1Bs1-SeuratObj.rds")
g <- seuratToGiottoV4(x)

# Read annotated scRNAseq
g_sc <- readRDS("data/annotated_scRNAseq/brca/HTA12_145.RDS")

##%######################################################%##
#                                                          #
####                     HTA12_246                      ####
#                                                          #
##%######################################################%##

# Read Visium sample
x <- readRDS("data/visium/brca/HT262B1-S1H2Fc2U1Z1Bs1-SeuratObj.rds")
g <- seuratToGiottoV4(x)

# Read annotated scRNAseq
g_sc <- readRDS("data/annotated_scRNAseq/brca/HTA12_246.RDS")

##%######################################################%##
#                                                          #
####                     HTA12_248                      ####
#                                                          #
##%######################################################%##

# Read Visium sample
x <- readRDS("data/visium/brca/HT308B1-S1H4Fc2U1Z1Bs1-SeuratObj.rds")
g <- seuratToGiottoV4(x)

# Read annotated scRNAseq
g_sc <- readRDS("data/annotated_scRNAseq/brca/HTA12_248.RDS")

##%######################################################%##
#                                                          #
####                     HTA12_249                      ####
#                                                          #
##%######################################################%##

# Read Visium sample
x <- readRDS("data/visium/brca/HT323B1-S1H3Fc2U1Z1Bs1-SeuratObj.rds")
g <- seuratToGiottoV4(x)

# Read annotated scRNAseq
g_sc <- readRDS("data/annotated_scRNAseq/brca/HTA12_249.RDS")

##%######################################################%##
#                                                          #
####                     HTA12_250                      ####
#                                                          #
##%######################################################%##

# Read Visium sample
x <- readRDS("data/visium/brca/HT339B1-S1H3Fc2U2Bs2-SeuratObj.rds")
g <- seuratToGiottoV4(x)

# Read annotated scRNAseq
g_sc <- readRDS("data/annotated_scRNAseq/brca/HTA12_250.RDS")

##%######################################################%##
#                                                          #
####                     HTA12_251                      ####
#                                                          #
##%######################################################%##

# Read Visium sample
x <- readRDS("data/visium/brca/HT365B1-S1H1Fc2U1Z1Bs1-SeuratObj.rds")
g <- seuratToGiottoV4(x)

# Read annotated scRNAseq
g_sc <- readRDS("data/annotated_scRNAseq/brca/HTA12_251.RDS")

##%######################################################%##
#                                                          #
####                     HTA12_253                      ####
#                                                          #
##%######################################################%##

# Read Visium sample
x <- readRDS("data/visium/brca/HT397B1-S1H3A1U21Bp1-SeuratObj.rds")
g <- seuratToGiottoV4(x)

# Read annotated scRNAseq
g_sc <- readRDS("data/annotated_scRNAseq/brca/HTA12_253.RDS")

##%######################################################%##
#                                                          #
####                     HTA12_262                      ####
#                                                          #
##%######################################################%##

# Read Visium sample
x <- readRDS("data/visium/brca/HT339B2-H1Fp1U1Bp1-SeuratObj.rds")
g <- seuratToGiottoV4(x)

# Read annotated scRNAseq
g_sc <- readRDS("data/annotated_scRNAseq/brca/HTA12_262.RDS")

##%######################################################%##
#                                                          #
####                     HTA12_266                      ####
#                                                          #
##%######################################################%##

# Read Visium sample
x <- readRDS("data/visium/brca/HT425B1-S1H2Fs1U1Bp1-SeuratObj.rds")
g <- seuratToGiottoV4(x)

# Read annotated scRNAseq
g_sc <- readRDS("data/annotated_scRNAseq/brca/HTA12_266.RDS")

##%######################################################%##
#                                                          #
####                     HTA12_271                      ####
#                                                          #
##%######################################################%##

# Read Visium sample
x <- readRDS("data/visium/brca/HT486B1-S1H1Fs1U1Bp1-SeuratObj.rds")
g <- seuratToGiottoV4(x)

# Read annotated scRNAseq
g_sc <- readRDS("data/annotated_scRNAseq/brca/HTA12_271.RDS")

##%######################################################%##
#                                                          #
####  Process Visium sample (filtering, normalization,  ####
####    dimension reduction, and Leiden clustering)     ####
#                                                          #
##%######################################################%##

g <- process_st_sample(g)

##%######################################################%##
#                                                          #
####               Run RCTD deconvolution               ####
#                                                          #
##%######################################################%##

### Create the Reference object
sc_counts <- getExpression(g_sc,
                           values = "raw",
                           output = "matrix")

cell_metadata <- pDataDT(g_sc)
sc_cell_types <- factor(cell_metadata$cell_types_final)
names(sc_cell_types) <- cell_metadata$cell_ID

reference_se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(counts = sc_counts),
  colData = data.frame(cell_type = sc_cell_types)
)

# Create SpatialExperiment

spatial_locs <- getSpatialLocations(g,
                                    output = "data.table")

st_coords <- data.frame(x = spatial_locs$sdimx,
                        y = spatial_locs$sdimy)
rownames(st_coords) <- spatial_locs$cell_ID

st_counts <- getExpression(g,
                           values = "raw",
                           output = "matrix")

spatial_spe <- SpatialExperiment::SpatialExperiment(
  assay = st_counts,
  spatialCoords = as.matrix(st_coords)
)

# Create and run RCTD
myRCTD <- createRctd(spatial_experiment = spatial_spe, 
                     reference_experiment = reference_se
)

results_spe <- runRctd(myRCTD, rctd_mode = "full", max_cores = 1)

# Get results
results_full <- results_spe@assays@data$weights

# Save deconvoluted_expression for ONTraC analysis
deconvoluted_expression <- results_spe@metadata$cell_type_info$info[[1]]
deconvoluted_expression <- t(deconvoluted_expression)
deconvoluted_expression <- cbind(data.frame(Spot_ID = rownames(deconvoluted_expression)),
                                 deconvoluted_expression)
write.csv(deconvoluted_expression,
          file = "ontrac_files/input_files/brca/sample_name/deconvoluted_expression.csv",
          quote = FALSE,
          row.names = FALSE)

# Normalize weights
x <- as.data.frame(t(as.matrix(results_full)))
norm_weights <- sweep(x, 1, rowSums(x), '/') 
norm_weights <- cbind(data.frame(cell_ID = rownames(results_spe@colData)),
                      norm_weights)

## Add RCTD to the Giotto object
spatEnrObj <- createSpatEnrObj(enrichment_data = norm_weights,
                               name = "RCTD",
                               method = "RCTD")

g <- subsetGiotto(g,
                  cell_ids = norm_weights$cell_ID)

g <- setSpatialEnrichment(gobject = g,
                          x = spatEnrObj,
                          name = "RCTD")

##%######################################################%##
#                                                          #
####             Plot deconvoluted sample               ####
#                                                          #
##%######################################################%##

spatDeconvPlot(g, 
               deconv_name = "RCTD",
               show_image = FALSE,
               pie_scale = 0.40,
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
                                   "Tumor" = "#7c7d80"))
