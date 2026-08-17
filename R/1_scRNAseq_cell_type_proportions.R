library(Giotto)
library(tidyverse)
instructions <- createGiottoInstructions()

##%######################################################%##
#                                                          #
####                    Read samples                    ####
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

##%######################################################%##
#                                                          #
####                    Join objects                    ####
#                                                          #
##%######################################################%##

g_sc <- joinGiottoObjects(
  gobject_list = list(HTA12_129, HTA12_135, HTA12_139,
                      HTA12_142, HTA12_144, HTA12_145, HTA12_246, 
                      HTA12_248, HTA12_249, HTA12_250, HTA12_251, 
                      HTA12_253, HTA12_262, HTA12_266, HTA12_271),
  gobject_names = c("HTA12_129", "HTA12_135", "HTA12_139",
                    "HTA12_142", "HTA12_144", "HTA12_145", "HTA12_246", 
                    "HTA12_248", "HTA12_249", "HTA12_250", "HTA12_251", 
                    "HTA12_253", "HTA12_262", "HTA12_266", "HTA12_271"))

##%######################################################%##
#                                                          #
####            Global cell type proportion             ####
#                                                          #
##%######################################################%##

cell_metadata <- pDataDT(g_sc)

global_cell_type_percentage <- cell_metadata %>% 
  count(cell_types_final) %>% 
  mutate(total = sum(n),
         percentage = n/total)

global_cell_type_percentage <- global_cell_type_percentage %>% 
  mutate(cell_types_final = fct_reorder(cell_types_final, desc(percentage))) 

ggplot(global_cell_type_percentage, 
       aes(x = "", y = percentage, fill = cell_types_final)) +
  geom_bar(stat = "identity", width = 1, color="white") +
  coord_polar("y", start = 0) +
  theme_void() +
  scale_fill_manual(values = c("Adipocytes" = "#d6d62d",
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
                               "Tumor" = "#7c7d80")) +
  theme(legend.title = element_blank())
ggsave(filename = "plots_paper/1_global_celltype_percentage.svg",
       device = "svg", scale = 2)

##%######################################################%##
#                                                          #
####          Cell type proportion per sample           ####
#                                                          #
##%######################################################%##

cell_metadata <- pDataDT(g_sc)

cell_type_percentage <- cell_metadata %>% 
  group_by(list_ID) %>% 
  count(cell_types_final) %>% 
  mutate(total = sum(n),
         percentage = n/total)

##%######################################################%##
#                                                          #
####            Tumor proportion per sample             ####
#                                                          #
##%######################################################%##

cell_type_percentage %>% 
  filter(cell_types_final == "Tumor") %>% 
  ggplot(aes(x = percentage, y = reorder(list_ID, percentage), fill = cell_types_final)) +
  geom_col(show.legend = FALSE) +
  scale_fill_manual(values = c("Tumor" = "#7c7d80")) +
  theme_bw() +
  labs(x = "", y = "") +
  scale_x_continuous(labels = scales::label_percent()) 
ggsave(filename = "plots_paper/1_sample_tumor_percentage.pdf",
       device = "pdf", scale = 1)

##%######################################################%##
#                                                          #
####     Non tumor cell types proportion per sample     ####
#                                                          #
##%######################################################%##

# Normalized to non-tumor cells per sample
cell_metadata %>% 
  filter(!cell_types_final == "Tumor") %>% 
  group_by(list_ID) %>% 
  count(cell_types_final) %>% 
  mutate(total = sum(n),
         percentage = n/total) %>% 
  ggplot(aes(x = percentage, y = list_ID, fill = cell_types_final)) +
  geom_col(show.legend = FALSE) +
  scale_fill_manual(values = c("Adipocytes" = "#d6d62d",
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
                               "T_cell" = "#169103")) +
  theme_bw() +
  labs(x = "", y = "") +
  scale_x_continuous(labels = scales::label_percent()) +
  scale_y_discrete(limits = rev(c("HTA12_142", "HTA12_271", "HTA12_262", "HTA12_145",
                                  "HTA12_248", "HTA12_144", "HTA12_139", "HTA12_253", 
                                  "HTA12_135", "HTA12_249", "HTA12_266", "HTA12_250",
                                  "HTA12_129",  "HTA12_251", "HTA12_246")))
ggsave(filename = "plots_paper/1_sample_celltypes_percentage.pdf",
       device = "pdf", scale = 1)


##%######################################################%##
#                                                          #
####            Add clinical classification             ####
#                                                          #
##%######################################################%##

# ESR1: ER+
# PGR: PR+
# ERBB2: HER2-Positive

pdf("plots_paper/1_clinical.pdf", width = 5, height = 3)
plot_norm_expression(g_sc,
                     features = c("ESR1", "PGR", "ERBB2"),
                     cluster_column = "list_ID")
dev.off()


library(genefu)
library(AnnotationDbi)
library(org.Hs.eg.db)

data(pam50.robust)

## Create data matrix for genefu
g_sc_expression <- as.matrix(getExpression(g_sc,
                                           output = "matrix",
                                           values = "raw"))

cell_metadata <- pDataDT(g_sc)
cell_metadata <- cell_metadata[cell_metadata$cell_types_final == "Tumor",]
g_sc_expression <- g_sc_expression[, cell_metadata$cell_ID]

ddata <- data.frame()
for(i in sort(unique(cell_metadata$list_ID))) {
  cell_ids <- cell_metadata[cell_metadata$list_ID == i,]$cell_ID
  g_sc_expression_small <- g_sc_expression[, cell_ids]
  ddata <- rbind(ddata, rowSums(g_sc_expression_small))
}

colnames(ddata) <- rownames(ddata)
rownames(ddata) <- sort(unique(cell_metadata$list_ID))

## Create gene annotation for genefu
gene_symbols <- colnames(sample_expression)
entrez_ids <- mapIds(
  x = org.Hs.eg.db,           # The annotation database
  keys = gene_symbols,        # Your input gene names
  keytype = "SYMBOL",         # The type of your input keys
  column = "ENTREZID",        # The type of ID you want to retrieve
  multiVals = "first"         # What to do if a gene maps to multiple IDs
)

dannot <- data.frame(EntrezGene.ID = entrez_ids,
                     Gene.symbol = gene_symbols,
                     probe = gene_symbols)

SubtypePredictions <- molecular.subtyping(sbt.model = "pam50", data = ddata,
                                          annot = dannot, do.mapping = FALSE,
                                          verbose = TRUE)

## Print results
SubtypePredictions$subtype

