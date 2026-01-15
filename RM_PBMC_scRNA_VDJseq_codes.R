# RM_PBMC_scRNA_VDJseq_codes
# Author: Andy Kwan Pui Chan 
# Date: 1 Dec 2025
# Description: This R document contained all the codes employed in the production of all the main figures of the manuscript. Raw and processed files are uploaded to NCBI GEO accordingly. All RDS objects will be available upon request.

# Library -----------
library("writexl")
library("readxl")
library("dplyr")
library("Seurat") 
library("SeuratData")
library("patchwork")
library("purrr")
library("tidyr")
library("ggpubr")
library("spatstat.utils")
library("tibble")
library("cowplot")
library("ggplot2")
library("fields")
library("ROCR")
library("KernSmooth")
library("Matrix")
library("parallel")
library("clustree")
library("DoubletFinder")
library('magrittr')
library('pheatmap')
library("EnhancedVolcano")
library("metap")
library("stringr")
library("stringi")
library("SeuratWrappers")
library("grid")
library("presto")
library("BPCells")
library('SingleCellExperiment')
library('escape')
library('dittoSeq')
library('DESeq2')
library("Azimuth")
library("scRepertoire")
library("gghighlight")
library("RColorBrewer")
library("seqinr")
library("gghighlight")
library("RColorBrewer")
library("ggpubr")
library("rstatix")
library("gridExtra")
library("clusterProfiler")

# Fig1 -------------
# Fig1A and B

Fig1A <- DimPlot(integrated_MeV, reduction = "umap.rpca", group.by = "predicted.celltype.l1", label = FALSE) +
  ggtitle ("Predicted PBMC\nSeurat clusters")+ 
  theme_classic(base_size = 8)

cell_count_per_txgrp_timepoint <- function (){
  
  plot_list <<- list()
  
  datatable_subset <- integrated_MeV@meta.data %>%
    group_by(predicted.celltype.l1, txgrp, DPI, RM_ID)%>%
    summarise (cell_count = n())%>%
    ungroup()%>%
    group_by(DPI, RM_ID)%>%
    mutate(total = sum(cell_count), percentage = cell_count/total*100) %>%
    ungroup() %>%
    filter (predicted.celltype.l1 != "other") %>%
    filter (predicted.celltype.l1 != "other T")
  
  mean_data <- datatable_subset %>%
    group_by(predicted.celltype.l1, txgrp, DPI)%>%
    summarise (mean_percentage = mean (percentage, na.rm= TRUE), .groups = "drop")
  
  plot1<- ggplot(datatable_subset, aes(x = as.factor(DPI), y = percentage, fill = as.factor (DPI))) +
    geom_boxplot(aes(fill=DPI))+
    geom_jitter(alpha = 0.5, size = 1)+
    theme_bw(base_size = 7)+
    theme (legend.position= "bottom")+
    labs (fill = "DPI")+
    ylab("Cellular proportion per timepoint (%)") +
    facet_wrap(predicted.celltype.l1~txgrp, ncol = 6) +
    scale_fill_manual(values=c("#00AFBB", "#99CC66", "#E7B800", "#FC4E07"))+
    theme (axis.title.x = element_blank(), 
           axis.text.x = element_blank(), 
           axis.ticks.x = element_blank())
  
  plot_list [[1]] <<- plot1

}
cell_count_per_txgrp_timepoint()
Fig1B <-grid.arrange(grobs = plot_list[1], ncol = 1)

ggsave("Fig1A.tiff",
       plot = Fig1A,
       dpi = 600, 
       height = 2.5,
       width  = 3, 
       units = "in")   

ggsave("Fig1B.tiff",
       plot = Fig1B,
       dpi = 600, 
       height = 2.5,
       width  = 4, 
       units = "in")  

Fig1C <-DimPlot(integrated_MeV_B_cells, group.by = c("RNA_snn_res_B_0.4"), label.size = 4, label = TRUE, repel = TRUE, reduction = "umap")+
  scale_color_discrete(labels = c(
    "0 - M-MBC 1",
    "1 - CD11c+ MBC 1",
    "2 - Naive B cells 1",
    "3 - C-MBC 1",
    "4 - M-MZL BC",
    "5 - CD11c+ MBC 2",
    "6 - M-MBC 2",
    "7 - Activated B cells",
    "8 - C-MBC 2",
    "9 - Contamination",
    "10 - Naive B cells 2",
    "11 - Naive B cells 3",
    "12 - Antibody-secreting cells"
  ))+
  ggtitle ("B cell subtypes")+ 
  theme_classic(base_size = 8)

cell_count_per_txgrp_timepoint_B_cell <- function (){
  
  plot_list <<- list()
  
  datatable_subset <- integrated_MeV_B_cells@meta.data %>%
    group_by(RNA_snn_res_B_0.4, txgrp, DPI, RM_ID)%>%
    summarise (cell_count = n())%>%
    ungroup()%>%
    group_by(DPI, RM_ID)%>%
    mutate(total = sum(cell_count), percentage = cell_count/total*100) %>%
    ungroup()%>%
    mutate (RNA_snn_res_B_0.4 = factor (RNA_snn_res_B_0.4, 
                                         levels = c("0", "1", "2", "3", "4", "5", "7", "8", "12"), 
                                         labels = c("M-MBC 1", "CD11c MBC 1" ,"Naive B", "C-MBC 1", "M-MZL BC", "CD11c MBC 2", "Activated B", "C-MBC 2", "ASC")))
  
  mean_data <- datatable_subset %>%
    group_by(RNA_snn_res_B_0.4, txgrp, DPI)%>%
    summarise (mean_percentage = mean (percentage, na.rm= TRUE), .groups = "drop")
  
  plot1<- ggplot(subset (datatable_subset, RNA_snn_res_B_0.4 %in% c("M-MBC 1", "CD11c MBC 1" ,"Naive B", "C-MBC 1", "M-MZL BC", "CD11c MBC 2", "Activated B", "C-MBC 2", "ASC")), aes(x = as.factor(DPI), y = percentage, fill = as.factor(DPI))) +
    geom_boxplot ()+
    geom_jitter(alpha = 0.5, size = 1)+
    theme_bw(base_size = 7) +
    theme (legend.position= "bottom")+
    labs (fill = "DPI")+
    ylab("Cellular proportion per timepoint (%)") +
    facet_wrap(RNA_snn_res_B_0.4~txgrp, ncol = 6) +
    scale_fill_manual(values=c("#00AFBB", "#99CC66", "#E7B800", "#FC4E07")) +
    theme(axis.text.x = element_text(angle = 90))+
    theme (axis.title.x = element_blank(), 
           axis.text.x = element_blank(), 
           axis.ticks.x = element_blank())
  
  plot_list [[1]] <<- plot1
  
}
cell_count_per_txgrp_timepoint_B_cell()
Fig1D <-grid.arrange(grobs = plot_list[1], ncol = 1)


ggsave("Fig1C.tiff",
       plot = Fig1C,
       dpi = 600, 
       height = 2.5,
       width  = 4, 
       units = "in")   

ggsave("Fig1D.tiff",
       plot = Fig1D,
       dpi = 600, 
       height = 3,
       width  = 4, 
       units = "in")  



Fig1E <-DimPlot(integrated_MeV_T_cells_baseline, group.by = c("RNA_snn_res_T_0.35"), label.size = 4, label = TRUE, repel = TRUE, reduction = "umap")+
  scale_color_discrete(labels = c(
    "0 - CD4 naive T cells",
    "1 - CD8 naive T cells",
    "2 - CD4 circulating follicular helper T cells",
    "3 - CD8 cytotoxic effector memory T cells",
    "4 - CD8 EOMES+ TIGIT+ T cells",
    "5 - CD4 Th2-like helper T cells",
    "6 - BACH2+ cells",
    "7 - ZEB2+ MYO1E+ cells",
    "8 - Contamination - APC",
    "9 - XCL1+ innate-like T cells",
    "10 - Contamination - CD16+ NK cells",
    "11 - CD4 Treg cells",
    "12 - TIMD4+ T cells",
    "13 - ZBTB16+ innate-like T cells",
    "14 - CD4 Th17-like helper T cells",
    "15 - Contamination - platelets",
    "16 - Proliferating T cells",
    "17 - Unknown"
  ))+
  ggtitle ("T cell subtypes")+ 
  theme_classic(base_size = 8)


cell_count_per_txgrp_timepoint_T_cell <- function (){
  
  plot_list <<- list()
  
  datatable_subset <- integrated_MeV_T_cells@meta.data %>%
    group_by(RNA_snn_res_T_0.35, txgrp, DPI, RM_ID)%>%
    summarise (cell_count = n())%>%
    ungroup()%>%
    group_by(DPI, RM_ID)%>%
    mutate(total = sum(cell_count), percentage = cell_count/total*100) %>%
    ungroup() %>%
    mutate (RNA_snn_res_T_0.35 = factor (RNA_snn_res_T_0.35, 
                                         levels = c("0", "1", "2", "3", "4", "5"), 
                                         labels = c("CD4 naive", "CD8 naive" ,"CD4 TFH", "CD8 cytotoxic", "CD8 EOMES+", "CD4 Th2-like")))
  
  mean_data <- datatable_subset %>%
    group_by(RNA_snn_res_T_0.35, txgrp, DPI)%>%
    summarise (mean_percentage = mean (percentage, na.rm= TRUE), .groups = "drop")
  
  plot1<- ggplot(subset (datatable_subset, RNA_snn_res_T_0.35 %in% c("CD4 naive", "CD8 naive" ,"CD4 TFH", "CD8 cytotoxic", "CD8 EOMES+", "CD4 Th2-like")), aes(x = as.factor(DPI), y = percentage, fill = as.factor(DPI))) +
    geom_boxplot ()+
    geom_jitter(alpha = 0.5, size = 1)+
    theme_bw(base_size = 7) +
    theme (legend.position= "bottom")+
    labs(fill = "DPI") +
    ylab("Cellular proportion per timepoint (%)") +
    facet_wrap(RNA_snn_res_T_0.35~txgrp, ncol = 6) +
    scale_fill_manual(values=c("#00AFBB", "#99CC66", "#E7B800", "#FC4E07")) +
    theme (axis.title.x = element_blank(), 
           axis.text.x = element_blank(), 
           axis.ticks.x = element_blank())
  
  plot_list [[1]] <<- plot1

}
cell_count_per_txgrp_timepoint_T_cell()
Fig1F <-grid.arrange(grobs = plot_list[1], ncol = 1)
  


ggsave("Fig1E.tiff",
       plot = Fig1E,
       dpi = 600, 
       height = 2.5,
       width  = 4.5, 
       units = "in") 

ggsave("Fig1F_2.tiff",
       plot = Fig1F,
       dpi = 600, 
       height = 2.5,
       width  = 4, 
       units = "in") 



# T cell baseline object
genes_to_plot <- c(
  "CD3D", "TRAC","CD4", "SELL", "CCR7","IL7R", "TCF7", "SOX4", "SOCS3",
  "CXCR5","CR1", "ICOS","CTLA4", "IL6", "MAF", "BATF", "MX1", "AIM2", "ITGB1",
  "GATA3", "ADAM19", "DUSP4", "LGALS3","S100A4", "S100A6", "S100A11", "PLP2", "ANXA5", "CD82", "RUNX2", "TGFBR3", "PRDM1", "IL6R", "BHLHE40",
  "CCR6", "ADAM12", "IL18R1", "CD109", "IFNGR2", "DUSP16", "TNFRSF25","ITGB7",
  "FOXP3", "IL2RA", "TIGIT", "DUSP2",
  "CD8A", "CD8B", "PLAC8","CD7","ITGA1", "TMIGD2",
  "GZMB", "PRF1", "KLRD1", "KLRF1","TBX21", "NKG7", "CX3CR1", "CCL4L1", "CCL5", "PDCD1", "TNFRSF18",
  "EOMES", "GZMK", "GZMM", "GNLY", "SLAMF7", "NKG2D", "ITGAV", "CD74", "CD84", "CD96",  "CXCR4", "KLRB1", "CTSW",
  "XCL1", "ITGAX", "THY1",
  "TIMD4", "PRDM8", "ID3",
  "TRGV9", "ZBTB16", "KLRG1", "MKI67",
  "BACH2", "BCL2","BCL11B", "FOXO1", "ZBTB20", 
  "ZEB2", "MYO1E", "MYO1F"
)

Idents(integrated_MeV_T_cells_baseline) <- "RNA_snn_res_T_0.35"
alldata <- ScaleData(integrated_MeV_T_cells_baseline, features = unique(genes_to_plot), assay = "RNA")
alldata$new_identity <- "Unassigned"
alldata$new_identity[alldata$RNA_snn_res_T_0.35 == "0"] <- "CD4 naive T cells (Cluster 0)"
alldata$new_identity[alldata$RNA_snn_res_T_0.35 == "2"] <- "CD4 circulating follicular helper T cells (Cluster 2)"
alldata$new_identity[alldata$RNA_snn_res_T_0.35 == "5"] <- "CD4 Th2-like helper T cells (Cluster 5)"
alldata$new_identity[alldata$RNA_snn_res_T_0.35 == "14"] <- "CD4 Th17-like helper T cells (Cluster 14)"
# alldata$new_identity[alldata$RNA_snn_res_T_0.35 == "17"] <- "CD4 effector-like memory T cells (Cluster 17)"
alldata$new_identity[alldata$RNA_snn_res_T_0.35 == "11"] <- "CD4 Treg cells (Cluster 11)"
alldata$new_identity[alldata$RNA_snn_res_T_0.35 == "1"] <- "CD8 naive T cells (Cluster 1)"
alldata$new_identity[alldata$RNA_snn_res_T_0.35 == "3"] <- "CD8 cytotoxic effector memory T cells (Cluster 3)"
alldata$new_identity[alldata$RNA_snn_res_T_0.35 == "4"] <- "CD8 EOMES+ TIGIT+ T cells (Cluster 4)"
alldata$new_identity[alldata$RNA_snn_res_T_0.35 == "9"] <- "XCL1+ innate-like T cells(Cluster 9)"
alldata$new_identity[alldata$RNA_snn_res_T_0.35 == "12"] <- "TIMD4+ T cells (Cluster 12)"
alldata$new_identity[alldata$RNA_snn_res_T_0.35 == "13"] <- "ZBTB16+ innate-like T cells (Cluster 13)"
alldata$new_identity[alldata$RNA_snn_res_T_0.35 == "16"] <- "Proliferating T cells (Cluster 16)"
alldata$new_identity[alldata$RNA_snn_res_T_0.35 == "6"] <- "BACH2+ cells (Cluster 6)"
alldata$new_identity[alldata$RNA_snn_res_T_0.35 == "7"] <- "ZEB2+ MYO1E+ cells (Cluster 7)"
alldata$new_identity <- factor (alldata$new_identity, levels = rev(c("CD4 naive T cells (Cluster 0)",  "CD4 circulating follicular helper T cells (Cluster 2)", "CD4 Th2-like helper T cells (Cluster 5)",  "CD4 Th17-like helper T cells (Cluster 14)", "CD4 Treg cells (Cluster 11)", "CD8 naive T cells (Cluster 1)", "CD8 cytotoxic effector memory T cells (Cluster 3)", "CD8 EOMES+ TIGIT+ T cells (Cluster 4)","XCL1+ innate-like T cells(Cluster 9)",  "TIMD4+ T cells (Cluster 12)", "ZBTB16+ innate-like T cells (Cluster 13)", "Proliferating T cells (Cluster 16)" ,"BACH2+ cells (Cluster 6)", "ZEB2+ MYO1E+ cells (Cluster 7)", "Unassigned")))
Idents(alldata) <- "new_identity"

Fig1G <- DotPlot(subset(alldata, new_identity != "Unassigned"), features = genes_to_plot, group.by = "new_identity", assay = "RNA", dot.scale = 3) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
  scale_colour_gradient2(low = "blue", mid = "white", high = "red")+
  theme(text  = element_text(size = 6),
        axis.text  = element_text(size = 5),
        axis.text.x  = element_text(size = 5),
        legend.position = "bottom")

ggsave("Fig1G.tiff",
       plot = Fig1G,
       dpi = 600, 
       height = 3.3,
       width  = 7, 
       units = "in") 

# Fig2 -------------
Fig2A <- DotPlot (integrated_MeV, features = c("IFI6", "IFI27", "RSAD2", "EPSTI1", "RIGI", "OAS2", "OAS3", "DHX58",  "STAT1","STAT2","IRF7","IRF9", "BCL3", "SOCS3", "GADD45B", "BCL2A1","ISG15","UBE2L6",   "HERC5", "HERC6","IFIT1", "IFIT2", "IFIT3",  "MX1",  "GBP3", "GBP7", "RNF213"), group.by = "txgrp_DPI", assay = "RNA")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
  scale_colour_gradient2(low = "blue", mid = "white", high = "red")+
  theme(text  = element_text(size = 8),
        axis.text  = element_text(size = 7),
        axis.text.x  = element_text(size = 7),
        legend.position = "bottom")

ggsave("Fig2A.tiff",
       plot = Fig2A,
       dpi = 600, 
       height = 2.5,
       width  = 7, 
       units = "in")  


Fig2B <- VlnPlot(subset(integrated_MeV, txgrp_DPI == "WTMeV_D10/11"), features = c("IFI6", "IFI27", "RSAD2", "EPSTI1", "RIGI", "OAS2", "OAS3", "DHX58", "STAT1","STAT2","IRF7","IRF9", "BCL3", "SOCS3", "GADD45B", "BCL2A1","ISG15","UBE2L6", "HERC5", "HERC6","IFIT1", "IFIT2", "IFIT3", "MX1",  "GBP3", "GBP7", "RNF213"), group.by = "predicted.celltype.l1", assay = "RNA", pt.size = 0, ncol = 9) &
  theme(
    axis.text.x = element_text(size = 4, angle = 45, vjust = 0.5),    # x-axis group labels
    axis.text.y = element_text(size = 4),    # y-axis numbers
    axis.title = element_text(size = 4), 
    plot.title = element_text (size = 7, face = "bold")
  )

ggsave("Fig2B.tiff",
       plot = Fig2B,
       dpi = 600, 
       height = 3,
       width  = 7, 
       units = "in")  


library("org.Mmu.eg.db")
DefaultAssay(integrated_MeV) <- "RNA"
Idents(integrated_MeV) <- "celltype_l1_txgrp_DPI"

GO_plot_function <- function (cell_type, txgrp_var, DPI_var_1, DPI_var_2 ) {
  
  DefaultAssay(integrated_MeV) <- "RNA"
  integrated_MeV$celltype_l1_txgrp_DPI <- paste(integrated_MeV$predicted.celltype.l1, integrated_MeV$txgrp_DPI, sep = "_")
  Idents(integrated_MeV) <- "celltype_l1_txgrp_DPI"
  
  txgrp_DPI_deg_1 <- FindMarkers(integrated_MeV, ident.1 = paste0(cell_type ,"_", txgrp_var, "_", DPI_var_2), ident.2 = paste0(cell_type, "_", txgrp_var, "_", DPI_var_1), test.use="wilcox", min.pct = 0.25, logfc.threshold = log(1))
  
  top_genes_lfc_output_subtype <- txgrp_DPI_deg_1 %>%
    mutate(gene = rownames(txgrp_DPI_deg_1)) %>%      
    filter (avg_log2FC > 0.40) %>%
    filter (pct.1 > 0.25) %>%
    filter (p_val_adj < 0.05)%>%
    arrange(desc(avg_log2FC), p_val_adj) %>%
    dplyr :: slice(1:50) 
  
  genes_to_test <- top_genes_lfc_output_subtype$gene
  
  GO_results <- enrichGO(gene= genes_to_test, OrgDb = "org.Mmu.eg.db", ont = "BP", keyType = "SYMBOL") #BP, MP or CC
  as.data.frame(GO_results)  
  
  graph <- plot(barplot(GO_results, showCategory = 10))
  
  return (graph)
  
}
Fig2C_1 <- GO_plot_function("Mono", "WTMeV", "D0", "D10/11")+
  theme(
    axis.text.x = element_text(size = 5),    # x-axis group labels
    axis.text.y = element_text(size = 5),    # y-axis numbers
    axis.title = element_text(size = 5), 
    legend.text = element_text (size = 4), 
    legend.title = element_text (size = 5),
    legend.position = "bottom"
  )

Fig2C_2 <- GO_plot_function("Mono", "WTMeV", "D0", "D42/43/56")+
  theme(
    axis.text.x = element_text(size = 5),    # x-axis group labels
    axis.text.y = element_text(size = 5),    # y-axis numbers
    axis.title = element_text(size = 5), 
    legend.text = element_text (size = 4), 
    legend.title = element_text (size = 5),
    legend.position = "bottom"
  )

Fig2C_3 <- GO_plot_function("Mono", "WTMeV", "D0", "D167/168/176")+
  theme(
    axis.text.x = element_text(size = 5),    # x-axis group labels
    axis.text.y = element_text(size = 5),    # y-axis numbers
    axis.title = element_text(size = 5), 
    legend.text = element_text (size = 4), 
    legend.title = element_text (size = 5),
    legend.position = "bottom"
  )

ggsave("Fig2C_1.tiff",
       plot = Fig2C_1,
       dpi = 600, 
       height = 3,
       width  = 2.33, 
       units = "in")  

ggsave("Fig2C_2.tiff",
       plot = Fig2C_2,
       dpi = 600, 
       height = 3,
       width  = 2.33, 
       units = "in")  

ggsave("Fig2C_3.tiff",
       plot = Fig2C_3,
       dpi = 600, 
       height = 3.3,
       width  = 2.5, 
       units = "in")  

Fig2D_1 <- DotPlot(subset(integrated_MeV, predicted.celltype.l1 == "CD4 T"), features = c("BATF", "SBNO2", "VAV1"), group.by = "txgrp_DPI", assay = "RNA")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
  scale_colour_gradient2(low = "blue", mid = "white", high = "red")+
  coord_flip()+
  theme(text = element_text(size = 6),
        axis.text = element_text(size = 5),
        axis.text.x = element_text(size = 5),
        legend.position = "right")

Fig2D_2 <- DotPlot(subset(integrated_MeV, predicted.celltype.l1 == "NK"), features = c("FOS", "JUN", "JUNB", "DUSP1", "CD69", "MAMU-A", "RGS1",  "CLCN7","CCL3","RELB","KLF4", "KLF10", "FAM3C", "STMN1"), group.by = "txgrp_DPI", assay = "RNA")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
  scale_colour_gradient2(low = "blue", mid = "white", high = "red")+
  coord_flip()+
  theme(text= element_text(size = 6),
        axis.text = element_text(size = 5),
        axis.text.x = element_text(size = 5),
        legend.position = "right")

Fig2D_3 <- DotPlot(subset(integrated_MeV, predicted.celltype.l1 == "Mono"), features = c("B4GALT4", "KCNH7", "ENSMMUG00000045208", "GPX1", "ALPK1", "SELL", "NLRC5", "CSF3R","SMCHD1","TNK2", "TGFB1", "IL6", "FCN1", "PANX1","ACOD1","TEX2", "IL7R", "LRRC25","FABP5", "LPL"), group.by = "txgrp_DPI", assay = "RNA")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
  scale_colour_gradient2(low = "blue", mid = "white", high = "red")+
  coord_flip()+
  theme(text  = element_text(size = 6),
        axis.text  = element_text(size = 5),
        axis.text.x  = element_text(size = 5),
        legend.position = "right")

ggsave("Fig2D_1.tiff",
       plot = Fig2D_1,
       dpi = 600, 
       height = 3,
       width  = 2.33, 
       units = "in")  

ggsave("Fig2D_2.tiff",
       plot = Fig2D_2,
       dpi = 600, 
       height = 3,
       width  = 2.33, 
       units = "in")  

ggsave("Fig2D_3.tiff",
       plot = Fig2D_3,
       dpi = 600, 
       height = 3,
       width  = 3, 
       units = "in")  



# Fig3 --------------
Read10X_new <- function(
    data.dir,
    gene.column = 2,
    cell.column = 1,
    unique.features = TRUE,
    strip.suffix = FALSE
) {
  full.data <- list()
  for (i in seq_along(along.with = data.dir)) {
    run <- data.dir[i]
    if (!dir.exists(paths = run)) {
      stop("Directory provided does not exist")
    }
    barcode.loc <- file.path(run, 'barcodes.tsv')
    gene.loc <- file.path(run, 'genes.tsv')
    features.loc <- file.path(run, 'features.tsv.gz')
    matrix.loc <- file.path(run, 'matrix.mtx')
    # Flag to indicate if this data is from CellRanger >= 3.0
    pre_ver_3 <- file.exists(gene.loc)
    if (!pre_ver_3) {
      addgz <- function(s) {
        return(paste0(s, ".gz"))
      }
      barcode.loc <- addgz(s = barcode.loc)
      matrix.loc <- addgz(s = matrix.loc)
    }
    if (!file.exists(barcode.loc)) {
      stop("Barcode file missing. Expecting ", basename(path = barcode.loc))
    }
    if (!pre_ver_3 && !file.exists(features.loc) ) {
      stop("Gene name or features file missing. Expecting ", basename(path = features.loc))
    }
    if (!file.exists(matrix.loc)) {
      stop("Expression matrix file missing. Expecting ", basename(path = matrix.loc))
    }
    data <- readMM(file = matrix.loc)
    cell.barcodes <- read.table(file = barcode.loc, header = FALSE, sep = '\t', row.names = NULL)
    if (ncol(x = cell.barcodes) > 1) {
      cell.names <- cell.barcodes[, cell.column]
    } else {
      cell.names <- readLines(con = barcode.loc)
    }
    if (all(grepl(pattern = "\\-1$", x = cell.names)) & strip.suffix) {
      cell.names <- as.vector(x = as.character(x = sapply(
        X = cell.names,
        FUN = ExtractField,
        field = 1,
        delim = "-"
      )))
    }
    if (is.null(x = names(x = data.dir))) {
      if (length(x = data.dir) < 2) {
        colnames(x = data) <- cell.names
      } else {
        colnames(x = data) <- paste0(i, "_", cell.names)
      }
    } else {
      colnames(x = data) <- paste0(names(x = data.dir)[i], "_", cell.names)
    }
    feature.names <- read.delim(
      file = ifelse(test = pre_ver_3, yes = gene.loc, no = features.loc),
      header = FALSE,
      stringsAsFactors = FALSE
    )
    if (any(is.na(x = feature.names[, gene.column]))) {
      warning(
        'Some features names are NA. Replacing NA names with ID from the opposite column requested',
        call. = FALSE,
        immediate. = TRUE
      )
      na.features <- which(x = is.na(x = feature.names[, gene.column]))
      replacement.column <- ifelse(test = gene.column == 2, yes = 1, no = 2)
      feature.names[na.features, gene.column] <- feature.names[na.features, replacement.column]
    }
    if (unique.features) {
      fcols = ncol(x = feature.names)
      if (fcols < gene.column) {
        stop(paste0("gene.column was set to ", gene.column,
                    " but feature.tsv.gz (or genes.tsv) only has ", fcols, " columns.",
                    " Try setting the gene.column argument to a value <= to ", fcols, "."))
      }
      rownames(x = data) <- make.unique(names = feature.names[, gene.column])
    }
    # In cell ranger 3.0, a third column specifying the type of data was added
    # and we will return each type of data as a separate matrix
    if (ncol(x = feature.names) > 2) {
      data_types <- factor(x = feature.names$V3)
      lvls <- levels(x = data_types)
      if (length(x = lvls) > 1 && length(x = full.data) == 0) {
        message("10X data contains more than one type and is being returned as a list containing matrices of each type.")
      }
      expr_name <- "Gene Expression"
      if (expr_name %in% lvls) { # Return Gene Expression first
        lvls <- c(expr_name, lvls[-which(x = lvls == expr_name)])
      }
      data <- lapply(
        X = lvls,
        FUN = function(l) {
          return(data[data_types == l, , drop = FALSE])
        }
      )
      names(x = data) <- lvls
    } else{
      data <- list(data)
    }
    full.data[[length(x = full.data) + 1]] <- data
  }
  # Combine all the data from different directories into one big matrix, note this
  # assumes that all data directories essentially have the same features files
  list_of_data <- list()
  for (j in 1:length(x = full.data[[1]])) {
    list_of_data[[j]] <- do.call(cbind, lapply(X = full.data, FUN = `[[`, j))
    # Fix for Issue #913
    list_of_data[[j]] <- as(object = list_of_data[[j]], Class = "CsparseMatrix")
  }
  names(x = list_of_data) <- names(x = full.data[[1]])
  # If multiple features, will return a list, otherwise
  # a matrix.
  if (length(x = list_of_data) == 1) {
    return(list_of_data[[1]])
  } else {
    return(list_of_data)
  }
}


count_data_43F_D11 <- Read10X_new(data.dir = paste0("43F-D11-MeV", "/outs/filtered_feature_bc_matrix"))
MeV_43F_D11 <- paste0 ("43F-D11_", count_data_43F_D11@Dimnames[[2]])

count_data_83H_D11 <- Read10X_new(data.dir = paste0("83H-D11-MeV", "/outs/filtered_feature_bc_matrix"))
MeV_83H_D11 <- paste0 ("83H-D11_", count_data_83H_D11@Dimnames[[2]])

MeV_infected<- c(MeV_43F_D11, MeV_83H_D11)

all_cell_ID <- colnames (integrated_MeV)

MeV_infected_status <- ifelse (all_cell_ID %in% MeV_infected, 
                               "MeV RNA +ve", "MeV RNA -ve")

integrated_MeV <- AddMetaData(integrated_MeV, 
                              metadata = MeV_infected_status, 
                              col.name = "MeV_infection_status")


# B cell object
all_cell_ID <- colnames (integrated_MeV_B_cells)

MeV_infected_status <- ifelse (all_cell_ID %in% MeV_infected, 
                               "MeV RNA +ve", "MeV RNA -ve")

integrated_MeV_B_cells <- AddMetaData(integrated_MeV_B_cells, 
                                      metadata = MeV_infected_status, 
                                      col.name = "MeV_infection_status")

# T cell object
all_cell_ID <- colnames (integrated_MeV_T_cells)

MeV_infected_status <- ifelse (all_cell_ID %in% MeV_infected, 
                               "MeV RNA +ve", "MeV RNA -ve")

integrated_MeV_T_cells <- AddMetaData(integrated_MeV_T_cells, 
                                      metadata = MeV_infected_status, 
                                      col.name = "MeV_infection_status")

# prepare tables
pt <- table(integrated_MeV$MeV_infection_status, integrated_MeV$predicted.celltype.l1, integrated_MeV$DPI, integrated_MeV$RM_ID)
pt <- as.data.frame(pt)%>%
  filter (Var3 == "D10/11" & Var4 %in% c("43F", "83H") & Var2 != "DC" & Var2 != "other" & Var2 != "other T")

pt_summary <- pt %>%
  group_by (Var2, Var4)%>%
  mutate(total = sum(Freq), percentage = Freq/total*100)

pt_summary_2 <- pt %>%
  group_by (Var2, Var4)%>%
  mutate(total = sum(Freq), percentage = Freq/total*100)%>%
  group_by (Var1)%>%
  summarise(average_percentage = sum(Freq)/sum(total))

pt$Var1 <- as.character(pt$Var1)

Fig3A <- ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 8) +
  geom_col(position = "dodge", width = 0.5)+
  xlab("Celltypes") +
  ylab("Cell count") +
  facet_wrap (~Var4)+
  theme(legend.title = element_blank(), 
        axis.text.x = element_text(angle = 45, vjust = 0.5))+
  ggtitle ("Count of 43F and 83H MeV RNA +ve \nand -ve at D10/11")+
  
  
  ggplot(subset (pt_summary, Var1 == "MeV RNA +ve"), aes(x = Var2, y = percentage)) +
  theme_bw(base_size = 8) +
  geom_col(width = 0.5, fill = "#00BFC4") +
  xlab("Celltypes") +
  ylab("Proportion (%)") +
  facet_wrap (~Var4)+
  theme(legend.title = element_blank(), 
        axis.text.x = element_text(angle = 45, vjust = 0.5))+
  ggtitle ("Proportion of 43F and 83H MeV \nRNA +ve at D10/11")


ggsave("Fig3A.tiff",
       plot = Fig3A,
       dpi = 600, 
       height = 2,
       width  = 7, 
       units = "in")  




FindMarkers_function_4 <- function (dataset, cell_type_1, cell_type_2, avglog2fc_value= a) {
  
  baseline_deg <- FindMarkers(dataset, ident.1 = cell_type_2, ident.2 = cell_type_1, test.use="wilcox", min.pct = 0.25, logfc.threshold = 0)
  
  my_volcano_plot <- function(clus) {
    clus$levels <- "NOT SIGNIFICANT" # Default state is no change
    clus$levels[clus$avg_log2FC > avglog2fc_value & clus$p_val_adj < 1e-10] <- "UP" # Upregulated
    clus$levels[clus$avg_log2FC < -avglog2fc_value & clus$p_val_adj < 1e-10] <- "DOWN" # Downregulated
    clus$levels[clus$p_val_adj > 1e-10] <- "NO" # Not significant
    
    clus$delabel <- NA
    clus$gene_symbol <- rownames(clus)
    
    clus$delabel[clus$levels == "UP"| clus$levels == "DOWN"] <- clus$gene_symbol[clus$levels == "UP" | clus$levels == "DOWN"]
    
    num_up <- sum(clus$levels == "UP")
    num_down <- sum(clus$levels == "DOWN")
    
    clus <- clus %>%
      mutate(p_val_adj = ifelse(p_val_adj < 1e-300, 1e-300, p_val_adj))
    
    p <- ggplot(data=clus, aes(x=avg_log2FC, y=-log10(p_val_adj), color=levels, label=delabel)) +
      geom_point(alpha=0.5) + 
      theme_minimal() +
      theme(axis.line = element_line(colour = "black"))+
      theme(plot.title = element_text(size = 0), 
            axis.text=element_text(size=8),
            axis.title=element_text(size=9), 
            legend.title= element_text(size = 8), 
            legend.text = element_text(size = 8))+
      theme(legend.position="none")+
      geom_text_repel(max.overlaps = Inf, size = 2.5) +
      scale_color_manual(values=c( "DOWN" = "blue", "UP"= "red","NOT SIGNIFICANT" = "grey"),
                         breaks = c("DOWN", "UP", "NOT SIGNIFICANT"),
                         labels = c(paste0("cluster ", cell_type_1), paste0("cluster ", cell_type_2),"other genes"), 
                         name = "Highly expressed in") +
      geom_vline(xintercept=c(-avglog2fc_value, avglog2fc_value), col="black", linetype = "dashed") +
      geom_hline(yintercept=-log10(1e-10), col="black", linetype = "dashed") +
      annotate("text", x=Inf, y=Inf, label=paste("# of genes in blue:", num_down, "\n# of genes in red:", num_up), 
               hjust=1.1, vjust=2, size=0, color="black")
    
    return(p)
  }
  
  plot <- list()
  
  plot[[1]] <- my_volcano_plot (baseline_deg)+ggtitle(paste0("Clus ", cell_type_1, " vs ", cell_type_2))
  
  do.call(gridExtra::grid.arrange, c(plot, ncol=1))
  
}

DefaultAssay(integrated_MeV) <- "RNA"
Idents(integrated_MeV) <- "MeV_infection_status"

# Compare all MeV +ve and -ve cells in the entire dataset
DefaultAssay(integrated_MeV) <- "RNA"
Idents(integrated_MeV) <- "MeV_infection_status"
Fig3B <- FindMarkers_function_4 (integrated_MeV, "MeV RNA +ve", "MeV RNA -ve", 1.5) 

# Compare all MeV +ve and -ve cells in D10/11 in the WT group
DefaultAssay(integrated_MeV) <- "RNA"
Idents(integrated_MeV) <- "MeV_infection_status"

integrated_MeV$txgrp_DPI_MeV_infection_status <- paste(integrated_MeV$txgrp, integrated_MeV$DPI, integrated_MeV$MeV_infection_status, sep = "_")

Idents(integrated_MeV) <- "txgrp_DPI_MeV_infection_status"

print(unique (integrated_MeV$txgrp_DPI_MeV_infection_status))

Fig3C <- FindMarkers_function_4 (integrated_MeV, "WTMeV_D10/11_MeV RNA +ve" , "WTMeV_D10/11_MeV RNA -ve" , 0.5) 

ggsave("Fig3B.tiff",
       plot = Fig3B,
       dpi = 600, 
       height = 3.5,
       width  = 3.5, 
       units = "in")

ggsave("Fig3C.tiff",
       plot = Fig3C,
       dpi = 600, 
       height = 3.5,
       width  = 3.5, 
       units = "in")

pt <- table(integrated_MeV_B_cells$MeV_infection_status, integrated_MeV_B_cells$RNA_snn_res_B_0.4, integrated_MeV_B_cells$DPI, integrated_MeV_B_cells$RM_ID)

pt <- as.data.frame(pt)%>%
  filter (Var3 == "D10/11" & Var4 %in% c("43F", "83H") & Var2 != "9" & Var2 != "10" & Var2 != "11")

pt$Var1 <- as.character(pt$Var1)

pt_summary_B <- pt %>%
  group_by (Var2, Var4)%>%
  mutate(total = sum(Freq), percentage = Freq/total*100)


pt <- table(integrated_MeV_T_cells$MeV_infection_status, integrated_MeV_T_cells$RNA_snn_res_T_0.35, integrated_MeV_T_cells$DPI, integrated_MeV_T_cells$RM_ID)

pt <- as.data.frame(pt)%>%
  filter (Var3 == "D10/11" & Var4 %in% c("43F", "83H") & Var2 != "6" & Var2 != "7" & Var2 != "8" & Var2 != "10" & Var2 != "15" & Var2 != "17")

pt$Var1 <- as.character(pt$Var1)

pt_summary_T <- pt %>%
  group_by (Var2, Var4)%>%
  mutate(total = sum(Freq), percentage = Freq/total*100)


Fig3D <- ggplot(subset (pt_summary_B, Var1 == "MeV RNA +ve"), aes(x = Var2, y = percentage)) +
  theme_bw(base_size = 8) +
  geom_col(width = 0.5, fill = "#00BFC4") +
  xlab("B cell clusters") +
  ylab("Proportion (%)") +
  facet_wrap (~Var4)+
  theme(legend.title = element_blank())+
  ggtitle ("Proportion of 43F and 83H MeV RNA +ve at D10/11")+

ggplot(subset (pt_summary_T, Var1 == "MeV RNA +ve"), aes(x = Var2, y = percentage)) +
  theme_bw(base_size = 8) +
  geom_col(width = 0.5, fill = "#00BFC4") +
  xlab("T cell clusters") +
  ylab("Proportion (%)") +
  facet_wrap (~Var4)+
  theme(legend.title = element_blank())

ggsave("Fig3D.tiff",
       plot = Fig3D,
       dpi = 600, 
       height = 2,
       width  = 7, 
       units = "in")  

# Fig4 ----------------------
heatmap_function <- function (dataset){
  df <- as.data.frame(dataset)
  rownames(df) <- df[,1]   
  df <- df[,-1]            
  m <- as.matrix(sapply(df, as.numeric))
  rownames(m) <- rownames(df)
  pheatmap(m, cluster_rows = FALSE, 
           cluster_cols = FALSE)
}

heatmap_function_2 <- function (dataset){
  df <- as.data.frame(dataset)
  rownames(df) <- df[,1]   
  df <- df[,-1]            
  m <- as.matrix(sapply(df, as.numeric))
  pheatmap(m, cluster_rows = FALSE, 
           cluster_cols = FALSE)
}

RM1_plot <- heatmap_function(df_RM1)
RM2_plot <- heatmap_function(df_RM2)
RM3_plot <- heatmap_function(df_RM3)

RM4_plot <- heatmap_function_2(df_RM4)
RM5_plot <- heatmap_function_2(df_RM5)
RM6_plot <- heatmap_function_2(df_RM6)

ggsave("Fig4a.tiff", plot = RM1_plot, width = 2.33, height = 2.33, dpi = 600, units = "in")
ggsave("Fig4b.tiff", plot = RM2_plot, width = 2.33, height = 2.33, dpi = 600, units = "in")
ggsave("Fig4c.tiff", plot = RM3_plot, width = 2.33, height = 2.33, dpi = 600, units = "in")
ggsave("Fig4d.tiff", plot = RM4_plot, width = 1.2, height = 1.5, dpi = 600, units = "in")
ggsave("Fig4e.tiff", plot = RM5_plot, width = 1.2, height = 1.5, dpi = 600, units = "in")
ggsave("Fig4f.tiff", plot = RM6_plot, width = 1.2, height = 1.5, dpi = 600, units = "in")

# Fig5 -------
Fig5A <- quantContig(combined, cloneCall="gene+nt", scale = T)+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), text = element_text(size = 8), legend.position = "none")  +
  quantContig(combined, cloneCall="gene+nt", scale = F)+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), text = element_text(size = 8))

ggsave("Fig5A.tiff",
       plot = Fig5A,
       dpi = 600, 
       height = 1.8,
       width  = 8, 
       units = "in")  

#BCR
Fig5B_1 <- compareClonotypes(combined[1:4], 
                  numbers = 20,
                  cloneCall="aa", 
                  graph = "alluvial")+
  ylab ("Clonal proportion")+
  scale_x_discrete(
  labels = c(
    "01-13F-D0_D0"   = "D0",
    "02-13F-D11_D10/11"  = "D10/11",
    "03-13F-D43_D42/43/56" = "D42/43/56",
    "04-13F-D168_D167/168/176"= "D167/168/176"))+
  theme(legend.position = "none", 
        text = element_text(size = 7), 
        axis.text.x = element_text(angle=90))

ggsave("Fig5B_1.tiff",
       plot = Fig5B_1,
       dpi = 600, 
       height = 2.33,
       width  = 2.33, 
       units = "in") 

Fig5B_2 <-compareClonotypes(combined[9:12], 
                  numbers = 20,
                  cloneCall="aa", 
                  graph = "alluvial")+
  ylab ("Clonal proportion")+
  scale_x_discrete(
    labels = c(
      "09-47G-D0_D0"   = "D0",
      "10-47G-D10_D10/11"  = "D10/11",
      "11-47G-D42_D42/43/56" = "D42/43/56",
      "12-47G-D168_D167/168/176"= "D167/168/176"))+
  theme(legend.position = "none", 
        text = element_text(size = 7), 
        axis.text.x = element_text(angle=90))

ggsave("Fig5B_2.tiff",
       plot = Fig5B_2,
       dpi = 600, 
       height = 2.33,
       width  = 2.33, 
       units = "in") 

Fig5B_3 <-compareClonotypes(combined[13:16], 
                  numbers = 20,
                  cloneCall="aa", 
                  graph = "alluvial")+
  ylab ("Clonal proportion")+
  scale_x_discrete(
    labels = c(
      "13-69G-D0_D0"   = "D0",
      "14-69G-D10_D10/11"  = "D10/11",
      "15-69G-D42_D42/43/56" = "D42/43/56",
      "16-69G-D168_D167/168/176"= "D167/168/176"))+
  theme(legend.position = "none", 
        text = element_text(size = 7), 
        axis.text.x = element_text(angle=90))

ggsave("Fig5B_3.tiff",
       plot = Fig5B_3,
       dpi = 600, 
       height = 2.33,
       width  = 2.33, 
       units = "in") 

Fig5B_4 <- compareClonotypes(combined[5:8], 
                  numbers = 20,
                  cloneCall="aa", 
                  graph = "alluvial")+
  ylab ("Clonal proportion")+
  scale_x_discrete(
    labels = c(
      "05-43F-D0_D0"   = "D0",
      "06-43F-D11_D10/11"  = "D10/11",
      "07-43F-D56_D42/43/56" = "D42/43/56",
      "08-43F-D176_D167/168/176"= "D167/168/176"))+
  theme(legend.position = "none", 
        text = element_text(size = 7), 
        axis.text.x = element_text(angle=90))

ggsave("Fig5B_4.tiff",
       plot = Fig5B_4,
       dpi = 600, 
       height = 2.33,
       width  = 2.33, 
       units = "in") 

Fig5B_5 <- compareClonotypes(combined[17:20], 
                  numbers = 20,
                  cloneCall="aa", 
                  graph = "alluvial")+
  ylab ("Clonal proportion")+
  scale_x_discrete(
    labels = c(
      "17-83H-D0_D0"   = "D0",
      "18-83H-D11_D10/11"  = "D10/11",
      "19-83H-D42_D42/43/56" = "D42/43/56",
      "20-83H-D167_D167/168/176"= "D167/168/176"))+
  theme(legend.position = "none", 
        text = element_text(size = 7), 
        axis.text.x = element_text(angle=90))

ggsave("Fig5B_5.tiff",
       plot = Fig5B_5,
       dpi = 600, 
       height = 2.33,
       width  = 2.33, 
       units = "in") 

Fig5B_6 <- compareClonotypes(combined[21:24], 
                  numbers = 20,
                  cloneCall="aa", 
                  graph = "alluvial")+
  ylab ("Clonal proportion")+
  scale_x_discrete(
    labels = c(
      "21-84H-D0_D0"   = "D0",
      "22-84H-D11_D10/11"  = "D10/11",
      "23-84H-D42_D42/43/56" = "D42/43/56",
      "24-84H-D167_D167/168/176"= "D167/168/176"))+
  theme(legend.position = "none", 
        text = element_text(size = 7), 
        axis.text.x = element_text(angle=90))

ggsave("Fig5B_6.tiff",
       plot = Fig5B_6,
       dpi = 600, 
       height = 2.33,
       width  = 2.33, 
       units = "in") 


clonalDiversity_table <- clonalDiversity(combined, cloneCall = "gene", group.by = "sample", n.boots = 100, exportTable = TRUE)

clonalDiversity_table_dataset <- clonalDiversity_table %>%
  separate(sample, into = c("position", "RM_ID", "DPI"), sep = "-", remove=FALSE)%>%
  mutate(txgrp = ifelse (RM_ID %in% c("13F", "47G", "69G"), "LAMV", "WTMeV"))%>%
  mutate(sex = ifelse (RM_ID %in% c("13F", "43F"), "M", "F"))%>%
  mutate(DPI = ifelse (DPI == "D10", "D10/11", DPI))%>%
  mutate(DPI = ifelse (DPI == "D11", "D10/11", DPI))%>%
  mutate(DPI = ifelse (DPI == "D42", "D42/43/56", DPI))%>%
  mutate(DPI = ifelse (DPI == "D43", "D42/43/56", DPI))%>%
  mutate(DPI = ifelse (DPI == "D56", "D42/43/56", DPI))%>%
  mutate(DPI = ifelse (DPI == "D167", "D167/168/176", DPI))%>%
  mutate(DPI = ifelse (DPI == "D168", "D167/168/176", DPI))%>%
  mutate(DPI = ifelse (DPI == "D176", "D167/168/176", DPI))%>%
  mutate(RM_ID = factor(RM_ID, levels = c("13F", "47G", "69G", "43F", "83H", "84H")))%>%
  mutate(position = factor(position, levels = c("01", "02", "03", "04", "09", "10", "11" ,"12" ,"13" ,"14" ,"15" ,"16", "05", "06", "07", "08", "17", "18", "19", "20", "21" ,"22", "23", "24")))%>%
  mutate(DPI = factor (DPI, levels =c('D0', 'D10/11' ,'D42/43/56', 'D167/168/176')))

Fig5C <- ggplot(clonalDiversity_table_dataset) + aes(x= DPI, y= Shannon, color = DPI)+
  geom_boxplot()+
  geom_jitter(alpha = 0.5, width = 0.15)+
  facet_wrap(~txgrp, ncol= 5)+
  theme_bw(base_size=9)+
  xlab("DPI")+
  ylab("Shannon diversity index")+
  labs(color="DPI")+
  scale_color_manual(values=c("#00AFBB", "#99CC66", "#E7B800", "#FC4E07"))+
  theme (legend.position = "none")+
  theme (axis.title.x = element_blank(), 
         axis.text.x = element_blank(), 
         axis.ticks.x = element_blank())+
  
  ggplot(clonalDiversity_table_dataset) + aes(x= DPI, y= Chao, color = DPI)+
  geom_boxplot()+
  geom_jitter(alpha = 0.5, width = 0.15)+
  facet_wrap(~txgrp, ncol= 5)+
  theme_bw(base_size=9)+
  xlab("DPI")+
  ylab("Chao diversity index")+
  labs(color="DPI")+
  scale_color_manual(values=c("#00AFBB", "#99CC66", "#E7B800", "#FC4E07"))+
  theme (axis.title.x = element_blank(), 
         axis.text.x = element_blank(), 
         axis.ticks.x = element_blank())

ggsave("Fig5C.tiff",
       plot = Fig5C,
       dpi = 600, 
       height = 2,
       width  = 7, 
       units = "in") 

# Fig6 -------------------
vizGenes_function_1 <- function (IG_chain, gene_segment) {
  dataset <-  vizGenes(combined, gene = gene_segment, chain = IG_chain, plot = "bar", order = "variance", scale = TRUE, exportTable = TRUE) %>%
    separate(element.names, into = c("sample_id", "DPI"), sep = "_" , remove = TRUE)%>%
    separate(Var2, into = c("position", "RM_ID", "DPI_orig"), sep = "-", remove=FALSE)%>%
    mutate(txgrp = ifelse (RM_ID %in% c("13F", "47G", "69G"), "LAMV", "WTMeV"))%>%
    mutate(sex = ifelse (RM_ID %in% c("13F", "43F"), "M", "F"))%>%
    mutate(DPI = ifelse (DPI == "D11", "D10/11", DPI))%>%
    mutate(DPI = ifelse (DPI == "D43/56", "D42/43/56", DPI))%>%
    mutate(DPI = ifelse (DPI == "D168/176", "D167/168/176", DPI))%>%
    mutate(RM_ID = factor(RM_ID, levels = c("13F", "47G", "69G", "43F", "83H", "84H")))%>%
    mutate(position = factor(position, levels = c("01", "02", "03", "04", "09", "10", "11" ,"12" ,"13" ,"14" ,"15" ,"16", "05", "06", "07", "08", "17", "18", "19", "20", "21" ,"22", "23", "24")))%>%
    mutate(DPI = factor (DPI, levels =c('D0', 'D10/11' ,'D42/43/56', 'D167/168/176')))%>%
    group_by(Var1, RM_ID) %>% 
    mutate(pct_change = (n/lag(n) - 1) * 100) %>%
    mutate(pct_change_from_baseline = (((n - n[1L])/n[1L]) *100))%>%
    ungroup()
  
  return(dataset)
  
}
vizGenes_function_2 <- function (dataset){
  dataset_2 <- dataset %>%
    group_by (Var1)%>%
    filter (any(n > 0.04)) %>%
    ungroup () %>%
    mutate(RM_ID = factor(RM_ID, levels= c("13F", "47G", "69G", "43F", "83H", "84H")))%>%
    mutate(RM_ID_DPI = paste0(RM_ID, "_", DPI))
  
  return (dataset_2)
  
}
vizGenes_function_3 <- function (dataset){
  dataset_3 <-  dataset%>%
    mutate (substring_gene = substr(Var1,1,5))%>%
    group_by(sample_id,RM_ID,DPI, substring_gene) %>%
    summarise(Freq = sum(n))
  
}
vizGenes_function_3b <- function (dataset){
  dataset_3 <-  dataset%>%
    mutate (substring_gene = substr(Var1,1,6))%>%
    mutate(substring_gene = str_remove(substring_gene, "[-*].*")) %>%
    group_by(sample_id, RM_ID, DPI, substring_gene) %>%
    summarise(Freq = sum(n))
  
}
vizGenes_IGHC_original <- vizGenes_function_1("IGH", "C")  %>%
  filter (Var1 != "IGLC1")
vizGenes_IGHC_selected <- vizGenes_function_2(vizGenes_IGHC_original)
vizGenes_IGHC_substring_gene <- vizGenes_function_3(vizGenes_IGHC_original) %>%
  filter (substring_gene != "IGLC1")

Fig6A <- ggplot(vizGenes_IGHC_original, aes(x= Var1, y= n))+
  geom_bar(stat = "identity", position= position_dodge(), aes(fill = DPI))+
  theme_bw(base_size=9)+
  xlab(paste0("IGHC", " Genes"))+
  ylab("Proportions")+
  labs(fill="DPI") +
  theme(axis.text.x=element_text(angle=90))+
  scale_fill_manual(values= c("#00AFBB", "#99CC66", "#E7B800", "#FC4E07"), labels = c('D0', 'D10/11' ,'D42/43/56', 'D167/168/176'))+
  facet_wrap(~RM_ID, ncol = 3)


Fig6B <- ggplot(subset(vizGenes_IGHC_selected, Var1 != "IGHG2")) + aes(x= DPI, y= n, color = DPI)+
  geom_boxplot()+
  geom_jitter(alpha = 0.5, width = 0.15)+
  facet_wrap(txgrp~Var1, ncol= 4)+
  theme_bw(base_size=9)+
  xlab("DPI")+
  ylab("Proportions")+
  labs(color="DPI")+
  theme(axis.text.x=element_text(angle=90))+
  stat_compare_means(method = "t.test", label = "p.signif", comparison = list(c("D0", "D10/11"), c("D0", "D42/43/56"), c("D0", "D167/168/176"), c("D10/11", "D42/43/56"), c("D10/11", "D167/168/176"), c("D42/43/56", "D167/168/176")), paired= TRUE, hide.ns = TRUE, tip.length = 0)+
  scale_color_manual(values=c("#00AFBB", "#99CC66", "#E7B800", "#FC4E07"))+
  theme (axis.title.x = element_blank(), 
         axis.text.x = element_blank(), 
         axis.ticks.x = element_blank())


ggsave("Fig6A.tiff",
       plot = Fig6A,
       dpi = 600, 
       height = 3,
       width  = 7, 
       units = "in")  

ggsave("Fig6B.tiff",
       plot = Fig6B,
       dpi = 600, 
       height = 3,
       width  = 7, 
       units = "in")  



SHM_per_RM_function <- function (dataset, RM_ID_variable, chain_variable){
  stat.test <- dataset%>%
    filter(RM_ID == RM_ID_variable & chain == chain_variable)%>%
    t_test(shm ~ DPI) %>%
    adjust_pvalue(method  = "bonferroni")%>%
    add_significance()%>%
    add_xy_position(x = "DPI")
  
  plot <- ggboxplot(subset(dataset,RM_ID == RM_ID_variable & chain == chain_variable), x = "DPI", y = "shm", color = "DPI", palette =c("#00AFBB", "#99CC66", "#E7B800", "#FC4E07"), size= 0.5)+
    geom_jitter(position=position_jitter(0.15), cex=1.5, alpha=0.008)+
    stat_pvalue_manual(stat.test, label = "p.adj.signif", y.position = c(7,8,9,10,11,12), size = 2.5)+
    ggtitle(paste0(RM_ID_variable, " ", chain_variable))+
    theme_bw(base_size=7.5)+
    theme(axis.text.x=element_blank(), 
          axis.ticks.x=element_blank())+
    theme (legend.position = "none")+
    ylab ("SHM rate")
  
  return (plot)
}

rm_ids <- unique(SHM_results_combined_transformed$RM_ID)
chains <- unique(SHM_results_combined_transformed$chain)

plots_list <- list()

for (rm_id in rm_ids){
  for (chain in chains){
    plot <- SHM_per_RM_function(SHM_results_combined_transformed, rm_id, chain)
    plots_list [[paste(rm_id, chain, sep = "_")]] <- plot
  }
}

Fig6C <- grid.arrange(grobs = plots_list[c(1, 7, 10, 4, 13, 16)], ncol =3)


Fig6D <- ggboxplot((subset(SHM_results_combined_transformed_2, chain == "IGH")), x = "DPI", y = "mean", color = "DPI", palette =c("#00AFBB", "#99CC66", "#E7B800", "#FC4E07"))+
  geom_jitter(position=position_jitter(0.15), cex=1.5, alpha=0.2)+
  ggtitle("")+
  facet_wrap(~txgrp, ncol =2)+
  theme_bw(base_size=8)+
  theme (legend.position = "none")+
  stat_compare_means(method = "t.test", label = "p.signif", comparison = list(c("D0", "D10/11"), c("D0", "D42/43/56"), c("D0", "D167/168/176"), c("D10/11", "D42/43/56"), c("D10/11", "D167/168/176"), c("D42/43/56", "D167/168/176")), paired= TRUE, hide.ns = TRUE, tip.length = 0)+
  ylab ("Average SHM rate")+
  xlab ("DPI")+
  theme (axis.text.x = element_blank(), 
         axis.ticks.x = element_blank())

ggsave("Fig6C.tiff",
       plot = Fig6C,
       dpi = 600, 
       height = 3,
       width  = 5, 
       units = "in")  

ggsave("Fig6D.tiff",
       plot = Fig6D,
       dpi = 600, 
       height = 2.3,
       width  = 2, 
       units = "in")  
# Fig7 -----------------
Fig7A <- Seurat::DimPlot(vgm[[2]],reduction = "umap", group.by = "VDJ_available", shuffle = T)+
  guides(color = "none", fill = "none")+
  theme_classic(base_size = 10) +
  ggtitle ("VDJ chain availability")

table(vgm[[2]]$VDJ_available)
pt <- table(vgm[[2]]$VDJ_available, vgm[[2]]$RM_ID, vgm[[2]]$seurat_clusters, vgm[[2]]$DPI)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)
pt$Var2 <- factor (pt$Var2, levels = c("13F", "47G", "69G", "43F", "83H", "84H"))
pt$Var3 <- factor (pt$Var3, levels = c(0:12))

Fig7B <- ggplot(pt, aes(x = Var4, y = Freq, fill = Var1)) +
  theme_bw(base_size = 9) +
  geom_col(position = "fill", width = 0.5) +
  xlab("DPI") +
  ylab("Proportions") +
  theme(legend.title = element_blank(), legend.position = "left")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_discrete (labels = c("Unavailable", "Available"))+    
  scale_x_discrete (labels = c("13F" = "RM_1", 
                               "47G" = "RM_2",
                               "69G" = "RM_3",
                               "43F" = "RM_4",
                               "83H" = "RM_5", 
                               "84H" = "RM_6"))+
  facet_wrap (~Var2)

ggsave("Fig7A.tiff",
       plot = Fig7A,
       dpi = 600, 
       height = 2.5,
       width  = 2.5, 
       units = "in")

ggsave("Fig7B.tiff",
       plot = Fig7B,
       dpi = 600, 
       height = 2.5,
       width  = 4.5, 
       units = "in")


vgm_gene_family_function <- function (dataset, IG_chain_gene_segment){
  vgm_all_total <- dataset %>%
    mutate (substring_gene = substr(!!sym(IG_chain_gene_segment),1,6))%>%
    mutate(substring_gene = str_remove(substring_gene, "[-*].*"))%>%
    group_by(seurat_clusters, substring_gene, DPI, txgrp) %>%
    summarize(n = sum(n))%>%
    group_by(seurat_clusters, DPI, txgrp) %>%
    mutate(total_2 = sum(n), Freq = n/total_2*100)
  
  return (vgm_all_total)
}
vgm_gene_family_VDJ_cgene <- vgm_gene_family_function(vgm_VDJ_cgene, "VDJ_cgene")
plot_3b_gene_family_prop_per_cluster_function <- function (dataset, IG_chain_gene_segment) {
  plot_d_prop <- ggplot(subset(dataset, !seurat_clusters %in% c(9)), aes(x=DPI , y= n))+
    geom_col(position = "fill", width = 0.5, aes(fill = substring_gene))+
    theme_bw(base_size=12)+
    xlab("DPI")+
    ylab("Proportions")+
    labs(fill=IG_chain_gene_segment)+
    facet_wrap(txgrp~seurat_clusters, ncol=12)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  
  return (plot_d_prop)
  
}

Fig7C <- plot_3b_gene_family_prop_per_cluster_function (vgm_gene_family_VDJ_cgene, "IGHC")  +theme_bw(base_size=8)+
  theme(axis.text.x=element_text(angle=90))


Fig7D <- ggplot(all_clone_pid_new, aes (x = DPI, y = Freq, fill= factor(mutational_state, levels = mutational_state_levels)))+
  geom_col(position = "fill", width = 0.5)+
  facet_wrap(txgrp~seurat_clusters, ncol = 12)+
  labs(fill="SHM levels")+
  scale_fill_manual(values = c( "#FC4E07", "#E7B800", "#99CC66", "#00AFBB"))+
  theme_bw(base_size=8)+
  theme(axis.text.x=element_text(angle=90))+
  ylab ("Proportions")

ggsave("Fig7C.tiff",
       plot = Fig6C,
       dpi = 600, 
       height = 3,
       width  = 7, 
       units = "in")

ggsave("Fig7D.tiff",
       plot = Fig6D,
       dpi = 600, 
       height = 3,
       width  = 7, 
       units = "in")


Fig7E <-ggboxplot(subset(all_clone_pid_IGH, seurat_clusters %in% c(2) & chain == "IGH"), x="DPI", y="CDRH3_length", color = "DPI", size=0.5, add = c("boxplot")) +
  geom_jitter(position=position_jitter(0.15), cex=1.5, alpha=0.01) +
  ggtitle("CDRH3 lengths of naive B cell cluster 2") +
  labs(x ="DPI", y = "CDRH3 length")+ 
  facet_grid(seurat_clusters~txgrp)+
  scale_color_manual(values=c("#00AFBB", "#99CC66", "#E7B800", "#FC4E07"))+
  theme_bw(base_size=8)+  
  theme (axis.title.x = element_blank(),
         axis.text.x = element_blank(),
         axis.ticks.x = element_blank(),
         strip.text.y= element_blank())+
  theme(legend.position = "right")

ggsave("Fig7E.tiff",
       plot = Fig6E,
       dpi = 600, 
       height = 1.5,
       width  = 5, 
       units = "in")

# Fig8 ---------
#TCR
Fig8A <- quantContig(combined_TCR, cloneCall="gene+nt", scale = T)+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), text = element_text(size = 8), legend.position = "none")  +
  quantContig(combined_TCR, cloneCall="gene+nt", scale = F)+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), text = element_text(size = 8)) 

ggsave("Fig8A.tiff",
       plot = Fig8A,
       dpi = 600, 
       height = 1.8,
       width  = 8, 
       units = "in")  


Clonal_expansion_donut_func <- function (dataset,chain_var, combined_var, sample_var) {
  
  clone_expansion_donut_2 <- function(sample_name, clone_exp){
    
    df <- clone_exp %>%
      # classify clones
      mutate(clone_id = ifelse(clone_count == 1, 'unique', clone_id)) %>%
      arrange(clone_count)
    
    # summarize into one row per clone (unique or expanded)
    df <- rbind(
      tibble(CDR3_concat = 'unique',
             Sample_Name = sample_name,
             clone_count = sum(df$clone_count == 1),
             clone_id = 'unique'),
      df %>% filter(clone_count > 1))
    
    # assign expansion categories
    df <- df %>%
      mutate(label = factor(case_when(
        clone_id == 'unique' ~ 'Singlet',
        clone_count == 2 ~ '2',
        clone_count == 3 ~ '3',
        clone_count == 4 ~ '4',
        clone_count == 5 ~ '5',
        clone_count > 5 & clone_count <= 9 ~ '6-9',
        clone_count > 9 & clone_count <= 19 ~ '10-19',
        clone_count > 19 ~ '>=20'),
        levels = c('Singlet', '2', '3', '4', '5', '6-9', '10-19', '>=20'))
      )
    
    # summarize by group
    df_grouped <- df %>%
      group_by(label) %>%
      summarise(clone_count = sum(clone_count), .groups = "drop") %>%
      arrange(label) %>%
      mutate(
        fraction = clone_count / sum(clone_count),
        ymax = cumsum(fraction),
        ymin = c(0, head(ymax, n = -1))
      )
    
    # y-axis label positions
    breaks <- df_grouped %>%
      mutate(pos = (ymax + ymin) / 2)
    
    # colors
    cols <- c("lightgrey", wesanderson::wes_palette("FantasticFox1", n = length(levels(df_grouped$label)) - 1, type = "continuous"))
    names(cols) <- levels(df_grouped$label)
    
    # plot
    g <- ggplot(df_grouped, aes(ymax = ymax, ymin = ymin, xmax = 5, xmin = 3.7, fill = label)) +
      geom_rect(color = 1, size = 0.1) +
      annotate("text", x = 2, y = 0, label = sum(df_grouped$clone_count), size = 5) +
      coord_polar(theta = "y") +
      xlim(c(2, 5)) +
      scale_y_continuous(breaks = breaks$pos, labels = breaks$label) +
      scale_fill_manual(values = cols) +
      theme(
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "right",
        panel.background = element_rect(fill = "white"),
        plot.title = element_text(hjust = 0.5)
      ) +
      labs(title = sample_name, fill = "Levels")
    
    return(g)
  }
  
  clone_exp <- dataset[[combined_var]] %>%
    filter(!is.na({{chain_var}}))%>%
    group_by(CTaa, sample) %>%
    summarise(clone_count = n()) %>% 
    ungroup() %>%
    dplyr::rename (CDR3_concat = CTaa, Sample_Name = sample)
  
  clone_exp$clone_id<-as.character(1:nrow(clone_exp))
  
  clone_expansion_donut_2(sample_var, clone_exp)
}


sample_vars <- c("13F-D0", "13F-D11", "13F-D43", "13F-D168", "43F-D0", "43F-D11", "43F-D56", "43F-D176", "47G-D0", "47G-D10", "47G-D42", "47G-D168", "69G-D0", "69G-D10", "69G-D42", "69G-D168", "83H-D0", "83H-D11", "83H-D42", "83H-D167", "84H-D0", "84H-D11", "84H-D42", "84H-D167")

plot_list <- list()

for (i in seq_along(sample_vars)) {
  sample_name <- sample_vars[i]
  combined_index <- i
  
  plot <- Clonal_expansion_donut_func(combined_TCR, "TCR2", combined_index, sample_name)
  plot_list[[sample_name]] <- plot
}

Fig8B_1 <- grid.arrange(grobs = plot_list[c(1:4, 9:16)], ncol =4)
Fig8B_2 <-grid.arrange(grobs = plot_list[c(5:8, 17:24)], ncol =4)

ggsave("Fig8B_1.tiff",
       plot = Fig8B_1,
       dpi = 600, 
       height = 7,
       width  = 13, 
       units = "in")  

ggsave("Fig8B_2.tiff",
       plot = Fig8B_2,
       dpi = 600, 
       height = 7,
       width  = 13, 
       units = "in")  


clonalDiversity_table <- clonalDiversity(combined_TCR, cloneCall = "gene", group.by = "sample", n.boots = 100, exportTable = TRUE)

clonalDiversity_table_dataset <- clonalDiversity_table %>%
  separate(sample, into = c("position", "RM_ID", "DPI"), sep = "-", remove=FALSE)%>%
  mutate(txgrp = ifelse (RM_ID %in% c("13F", "47G", "69G"), "LAMV", "WTMeV"))%>%
  mutate(sex = ifelse (RM_ID %in% c("13F", "43F"), "M", "F"))%>%
  mutate(DPI = ifelse (DPI == "D10", "D10/11", DPI))%>%
  mutate(DPI = ifelse (DPI == "D11", "D10/11", DPI))%>%
  mutate(DPI = ifelse (DPI == "D42", "D42/43/56", DPI))%>%
  mutate(DPI = ifelse (DPI == "D43", "D42/43/56", DPI))%>%
  mutate(DPI = ifelse (DPI == "D56", "D42/43/56", DPI))%>%
  mutate(DPI = ifelse (DPI == "D167", "D167/168/176", DPI))%>%
  mutate(DPI = ifelse (DPI == "D168", "D167/168/176", DPI))%>%
  mutate(DPI = ifelse (DPI == "D176", "D167/168/176", DPI))%>%
  mutate(RM_ID = factor(RM_ID, levels = c("13F", "47G", "69G", "43F", "83H", "84H")))%>%
  mutate(position = factor(position, levels = c("01", "02", "03", "04", "09", "10", "11" ,"12" ,"13" ,"14" ,"15" ,"16", "05", "06", "07", "08", "17", "18", "19", "20", "21" ,"22", "23", "24")))%>%
  mutate(DPI = factor (DPI, levels =c('D0', 'D10/11' ,'D42/43/56', 'D167/168/176')))

Fig8C <- ggplot(clonalDiversity_table_dataset) + aes(x= DPI, y= Shannon, color = DPI)+
  geom_boxplot()+
  geom_jitter(alpha = 0.5, width = 0.15)+
  facet_wrap(~txgrp, ncol= 5)+
  theme_bw(base_size=9)+
  xlab("DPI")+
  ylab("Shannon diversity index")+
  labs(color="DPI")+
  scale_color_manual(values=c("#00AFBB", "#99CC66", "#E7B800", "#FC4E07"))+
  theme (legend.position = "none")+
  theme (axis.title.x = element_blank(), 
         axis.text.x = element_blank(), 
         axis.ticks.x = element_blank())+
  
  ggplot(clonalDiversity_table_dataset) + aes(x= DPI, y= Chao, color = DPI)+
  geom_boxplot()+
  geom_jitter(alpha = 0.5, width = 0.15)+
  facet_wrap(~txgrp, ncol= 5)+
  theme_bw(base_size=9)+
  xlab("DPI")+
  ylab("Chao diversity index")+
  labs(color="DPI")+
  scale_color_manual(values=c("#00AFBB", "#99CC66", "#E7B800", "#FC4E07"))+
  theme (axis.title.x = element_blank(), 
         axis.text.x = element_blank(), 
         axis.ticks.x = element_blank())

ggsave("Fig8C.tiff",
       plot = Fig8C,
       dpi = 600, 
       height = 2,
       width  = 7, 
       units = "in") 

# Fig9 -------------
Fig9A_1 <- Seurat::DimPlot(vgm_t_baseline[[2]],reduction = "umap", group.by = "VDJ_available", shuffle = T)+
  guides(color = "none", fill = "none")+
  theme_classic(base_size = 8) +
  ggtitle ("VDJ chain availability")

print(table(vgm_t_baseline[[2]]$VDJ_available)) #11875 true and 1974 false cells 

pt <- table(vgm_t_baseline[[2]]$VDJ_available, vgm_t_baseline[[2]]$RM_ID, vgm_t_baseline[[2]]$seurat_clusters)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)
pt$Var2 <- factor (pt$Var2, levels = c("13F", "47G", "69G", "43F", "83H", "84H"))

Fig9A_2 <- ggplot(pt, aes(x = Var3, y = Freq*100, fill = Var1)) +
  theme_bw(base_size = 8) +
  geom_col(position = "fill", width = 0.5) +
  xlab("T cell cluster") +
  ylab("Proportion") +
  theme(legend.title = element_blank())+
  scale_fill_discrete (labels = c("Unavailable", "Available"))+
  theme (legend.position = "left")


ggsave("Fig9A_1.tiff",
       plot = Fig9A_1,
       dpi = 600, 
       height = 2,
       width  = 2, 
       units = "in")

ggsave("Fig9A_2.tiff",
       plot = Fig9A_2,
       dpi = 600, 
       height = 2,
       width  = 3.5, 
       units = "in")


vgm_t_baseline_GEX <- vgm_t_baseline [[2]]
vgm_t_baseline_GEX_subset <- subset(vgm_t_baseline_GEX, Nr_of_VDJ_chains == 1 & Nr_of_VJ_chains == 1)

table(vgm_t_baseline_GEX$Nr_of_VDJ_chains)

vgm_t_baseline_GEX_remove_VDJ_cgene_missing_value <- subset(vgm_t_baseline_GEX, subset = VDJ_cgene == "TRBC1" | VDJ_cgene == "TRBC2"| VDJ_cgene == "TRGC1"| VDJ_cgene == "TRGC2")

table(vgm_t_baseline_GEX_remove_VDJ_cgene_missing_value$Nr_of_VDJ_chains)
table(vgm_t_baseline_GEX_remove_VDJ_cgene_missing_value$Nr_of_VJ_chains)
table(vgm_t_baseline_GEX_remove_VDJ_cgene_missing_value$seurat_clusters)
table(vgm_t_baseline_GEX_remove_VDJ_cgene_missing_value$Nr_of_VDJ_chains, vgm_t_baseline_GEX_remove_VDJ_cgene_missing_value$VDJ_cgene)


Idents(vgm_t_baseline_GEX_remove_VDJ_cgene_missing_value) <- "VDJ_cgene"

Fig8B <- DimPlot(vgm_t_baseline_GEX_remove_VDJ_cgene_missing_value,reduction = "umap", group.by = "VDJ_cgene", shuffle = T) +
  theme_classic(base_size = 8) +
  ggtitle ("C genes on \nT cell Seurat object")

ggsave("Fig9B.tiff",
       plot = Fig9B,
       dpi = 600, 
       height = 2,
       width  = 2, 
       units = "in")

Idents(vgm_t_baseline_GEX_remove_VDJ_cgene_missing_value) <- "VDJ_cgene"

find_all_markers_VDJ_cgene <- FindAllMarkers(object = vgm_t_baseline_GEX_remove_VDJ_cgene_missing_value,
                                             only.pos = FALSE,
                                             logfc.threshold = 0.2)

top_genes_lfc_output_function <- function (find_all_markers_dataset, slice_val) {
  top_genes_lfc_output_dataset <- find_all_markers_dataset %>%
    group_by(cluster) %>%
    filter (avg_log2FC > 0.2) %>%
    filter (pct.1 > 0.2) %>%
    arrange(cluster, desc(avg_log2FC), p_val_adj) %>%
    dplyr :: slice(1:slice_val) %>%
    mutate(category = row_number()) %>%
    ungroup()
  return (top_genes_lfc_output_dataset)
}
top_genes_lfc_output <- top_genes_lfc_output_function(find_all_markers_VDJ_cgene, 15)
table(top_genes_lfc_output$cluster)
cluster_var <- as.character(unique(top_genes_lfc_output$cluster))

plots_list <- list()

for (i in cluster_var){
  plot <- ggplot(subset(top_genes_lfc_output, cluster == i), aes(x = reorder (gene, -category), y = avg_log2FC)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    facet_wrap(~ cluster, scales = "free", ncol = 5) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    geom_hline(yintercept = c(0, 0.25), linetype = "dashed")+
    theme(strip.text.x = element_text(size = 10, face = "bold"))+
    labs(x = NULL, y = NULL)  
  
  plots_list [[i]] <- plot
}

grid.arrange(grobs = plots_list, ncol =4, top = textGrob("Top genes per cluster", gp = gpar(fontsize = 16, fontface = "bold")), bottom = textGrob("Average Log2 Fold Change", gp = gpar(fontsize = 14)), left = textGrob("Gene", rot = 90, gp = gpar(fontsize = 14)))

Idents(vgm_t_baseline_GEX_remove_VDJ_cgene_missing_value) <- "VDJ_cgene"

alldata <- ScaleData(vgm_t_baseline_GEX_remove_VDJ_cgene_missing_value, features = as.character(unique(top_genes_lfc_output$gene)), assay = "RNA")

Fig9C <- DotPlot(subset(alldata, downsample = 1000), features = (as.character(unique(top_genes_lfc_output$gene))), group.by = "VDJ_cgene", assay = "RNA")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
  scale_colour_gradient2(low = "blue", mid = "white", high = "red")&
  theme(text  = element_text(size = 8),
        axis.text  = element_text(size = 7),
        axis.text.x  = element_text(size = 7)) +
  theme (legend.position = "bottom")


ggsave("Fig9C.tiff",
       plot = Fig9C,
       dpi = 600, 
       height = 3.5,
       width  = 7, 
       units = "in") 


vgm_t_gene_subtype_function <- function (TR_chain_gene_segment) {
  dataset <- vgm_t[[1]] %>%
    filter (Nr_of_VDJ_chains != 2)%>%
    filter (Nr_of_VJ_chains != 2) %>%
    filter(!!sym(TR_chain_gene_segment) != '') %>%
    filter(seurat_clusters != '') %>%
    group_by(orig.ident, sample_id, RM_ID, seurat_clusters, DPI, !!sym(TR_chain_gene_segment), txgrp) %>%
    summarize(n = n())  %>%
    group_by(seurat_clusters, sample_id, DPI) %>%
    mutate(ttl = sum(n), Freq = n/ttl*100)
  
  dataset_2 <- dataset %>%
    group_by (!!sym(TR_chain_gene_segment))%>%
    summarize(meanvalue = mean (n))%>%
    arrange(desc(meanvalue))
  
  dataset[[TR_chain_gene_segment]] <- factor(dataset[[TR_chain_gene_segment]], levels = dataset_2[[TR_chain_gene_segment]])
  
  return(dataset)
  
}

vgm_t_gene_family_function <- function (dataset, TR_chain_gene_segment){
  vgm_all_total <- dataset %>%
    mutate (substring_gene = substr(!!sym(TR_chain_gene_segment),1,6))%>%
    mutate(substring_gene = str_remove(substring_gene, "[-*].*"))%>%
    group_by(seurat_clusters, substring_gene, DPI, txgrp) %>%
    summarize(n = sum(n))%>%
    group_by(seurat_clusters, DPI, txgrp) %>%
    mutate(total_2 = sum(n), Freq = n/total_2*100)
  
  return (vgm_all_total)
}


vgm_t_VDJ_vgene <- vgm_t_gene_subtype_function("VDJ_vgene")%>%
  filter (sample_id != "s1" | VDJ_vgene != "TRBV6-3")

vgm_t_VJ_vgene <- vgm_t_gene_subtype_function("VJ_vgene")%>%
  filter (sample_id != "s1" | VJ_vgene != "TRAV1-2")


vgm_t_BVDJ_vgene <- vgm_t_gene_subtype_function("VDJ_vgene")%>%
  filter(substr(VDJ_vgene, 1,3 ) == "TRB")
vgm_t_AVJ_vgene <- vgm_t_gene_subtype_function("VJ_vgene")%>%
  filter(substr(VJ_vgene, 1,3 ) == "TRA")
vgm_t_gene_family_VDJ_vgene <- vgm_t_gene_family_function(vgm_t_VDJ_vgene, "VDJ_vgene")

vgm_t_gene_family_VJ_vgene <- vgm_t_gene_family_function(vgm_t_VJ_vgene, "VJ_vgene")

vgm_t_gene_family_BVDJ_vgene <- vgm_t_gene_family_function(vgm_t_BVDJ_vgene, "VDJ_vgene")

vgm_t_gene_family_AVJ_vgene <- vgm_t_gene_family_function(vgm_t_AVJ_vgene, "VJ_vgene")

plot_3c_gene_family_prop_per_cluster_function <- function (dataset, TR_chain_gene_segment) {
  plot_d_prop <- ggplot(subset(dataset, !seurat_clusters %in% c(6,7,8,10,15,17)), aes(x=seurat_clusters , y= n))+
    geom_col(position = "fill", width = 0.5, aes(fill = substring_gene))+
    theme_bw(base_size=7.5)+
    theme(axis.text.x  = element_blank()) +
    xlab("")+
    ylab("Proportions")+
    labs(fill=TR_chain_gene_segment)+
    facet_wrap(txgrp~DPI, ncol=4)
  
  return (plot_d_prop)
  
}

Fig9D <-plot_3c_gene_family_prop_per_cluster_function (subset(vgm_t_gene_family_VDJ_vgene, seurat_clusters %in% c(3)), "TRB/GV")+ ggtitle ("CD8 Cytotoxic T cells (cluster 3)")+
plot_3c_gene_family_prop_per_cluster_function (subset(vgm_t_gene_family_VJ_vgene, seurat_clusters %in% c(3)), "TRA/DV")

ggsave("Fig9D.tiff",
       plot = Fig9D,
       dpi = 600, 
       height = 4.5,
       width  = 9.2, 
       units = "in") 
