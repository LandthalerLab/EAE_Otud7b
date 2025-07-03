#Step 2

options(rgl.useNULL=TRUE)

#Install Voltron from github
devtools::install_github("Artur-man/VoltRon")

library(Seurat)
library(VoltRon)
library(dplyr)
library(ggplot2)


setwd("/path/to/working/directory/")


all_d0 <- readRDS("all_d0_mols.rds")
all_d15 <- readRDS("all_d15_mols.rds")

the_cols = c(neurons = "#CD9600", astrocytes = "#7CAE00", oligodendrocytes = "#00BF75", microglia = "#fc8d86", `T cells` = "#FF689F", `endothelial cells` = "#00BCD6", 
             pericytes = "#35A2FF", `ventral leptomeningeal cells` = "gray78")
the_cols_d0 = c(neurons = "#CD9600", astrocytes = "#7CAE00", oligodendrocytes = "#00BF75", microglia = "#fc8d86", `endothelial cells` = "#00BCD6", 
                pericytes = "#35A2FF", `ventral leptomeningeal cells` = "gray78")


#First work on d0 object


#Clustering
all_d0 <- normalizeData(all_d0, method = "LogNorm")
selected_features <- vrFeatures(all_d0)
all_d0 <- getPCA(all_d0, features = selected_features, dims = 30, overwrite = TRUE)
all_d0 <- getUMAP(all_d0, dims = 1:30, overwrite = TRUE)

#Make clusters 
all_d0 <- getProfileNeighbors(all_d0, dims = 1:16, k = 10, method = "SNN", data.type = "pca", graph.key = "SNN")

#We used a resolution that gave us 18 clusters but this is somewhat arbitrary
resols <- 1.1

all_d0 <- getClusters(all_d0, resolution = resols, label = paste0("Clusters_", resols), graph = "SNN")
vrEmbeddingPlot(all_d0, group.by = paste0("Clusters_", resols), embedding = "umap", pt.size=0.5)
ggsave(paste0("all_d0_clusters_", resols, ".pdf"))

#Calculate and plot cluster markers
all_d0_seu <- VoltRon::as.Seurat(all_d0, cell.assay = "Xenium", type = "image")
all_d0_seu <- NormalizeData(all_d0_seu)
Idents(all_d0_seu) <- paste0("Clusters_", resols)
markers <- FindAllMarkers(all_d0_seu)
topmarkers <- markers %>%
  filter(pct.1 > 0.4, avg_log2FC > 0.5) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)
DotPlot(all_d0_seu, group.by=paste0("Clusters_", resols),
        features = unique(topmarkers$gene, assay='Xenium')) +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1))
ggsave("all_d0_clustermarkers.pdf", height=length(unique(topmarkers$cluster))/2,
       width=length(unique(topmarkers$gene))/4, useDingbats=FALSE)

#Display cell type marker according to the 10x classification from https://www.10xgenomics.com/products/xenium-panels
marker_genes_file <- read.csv("Xenium_mBrain_v1.1_metadata.csv")
marker_genes <- marker_genes_file %>%
  filter(grepl("Astrocyte|Oligodendro|Neuron", Annotation, ignore.case=TRUE, perl=TRUE))
the_features <- as.character(marker_genes$Genes)
names(the_features) <- marker_genes$Annotation
DotPlot(all_d0_seu, group.by=paste0("Clusters_", resols),
        features = the_features) +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1))
ggsave("all_d0_xeniumgenes.pdf", height=length(unique(topmarkers$cluster))/2,
       width=length(unique(marker_genes$Genes))/4, useDingbats=FALSE)

#Get localization of clusters
vrSpatialPlot(all_d0, group.by = paste0("Clusters_", resols), alpha = 1, plot.segments = TRUE, background.color = "black")
ggsave("all_d0_Spatial_clusters.pdf", useDingbats=FALSE, height=6, width=6)


#Match clusters with cell types
#Note that the order of clusters can change between VoltRon instances, i.e. when you process the data again, cluster order might be different.
#We have used these marker genes for the different cell types. Localization can also be taken into account.
#Neurons: Ebf3, Penk, Necab1, Gad1, Gad2
#Pericytes: Acta2, Carmn, Ano1
#Endothelial cells: Cldn5, Adgrl4 but no Acta2
#Oligodendrocytes: Sox10, Opalin, Gjc3
#Astrocytes: Sox9, Rfx4, Wfs1
#Microglia: Trem2, Siglech, Laptm5
#Ventral leptomeningeal cells: Igfbp6, Col1a1, Dcn (also show strong peripheral localization)
all_d0@metadata@cell <- all_d0@metadata@cell %>% mutate(celltype = case_when(!!sym(paste0("Clusters_", resols)) == 1 ~ "neurons",
                                                                             !!sym(paste0("Clusters_", resols)) == 2 ~ "neurons",
                                                                             !!sym(paste0("Clusters_", resols)) == 3 ~ "endothelial cells",
                                                                             !!sym(paste0("Clusters_", resols)) == 4 ~ "oligodendrocytes",
                                                                             !!sym(paste0("Clusters_", resols)) == 5 ~ "astrocytes",
                                                                             !!sym(paste0("Clusters_", resols)) == 6 ~ "neurons",
                                                                             !!sym(paste0("Clusters_", resols)) == 7 ~ "microglia",
                                                                             !!sym(paste0("Clusters_", resols)) == 8 ~ "oligodendrocytes",
                                                                             !!sym(paste0("Clusters_", resols)) == 9 ~ "astrocytes",
                                                                             !!sym(paste0("Clusters_", resols)) == 10 ~ "pericytes",
                                                                             !!sym(paste0("Clusters_", resols)) == 11 ~ "ventral leptomeningeal cells",
                                                                             !!sym(paste0("Clusters_", resols)) == 12 ~ "oligodendrocytes",
                                                                             !!sym(paste0("Clusters_", resols)) == 13 ~ "oligodendrocytes",
                                                                             !!sym(paste0("Clusters_", resols)) == 14 ~ "ventral leptomeningeal cells",
                                                                             !!sym(paste0("Clusters_", resols)) == 15 ~ "astrocytes",
                                                                             !!sym(paste0("Clusters_", resols)) == 16 ~ "astrocytes",
                                                                             !!sym(paste0("Clusters_", resols)) == 17 ~ "astrocytes",
                                                                             !!sym(paste0("Clusters_", resols)) == 18 ~ "neurons",
                                                                             .default = "other"))

saveRDS(all_d0, "./all_d0_mols_celltypeann.rds")

#Visualize cell types
vrSpatialPlot(all_d0, colors = the_cols_d0, group.by = "celltype", alpha = 1, plot.segments = TRUE, background.color = "black")
ggsave("all_d0_Spatial_celltypes.pdf", useDingbats=FALSE, height=6, width=6)

#cell type markers at d0 and d15
all_d0_seu <- VoltRon::as.Seurat(all_d0, cell.assay = "Xenium", type = "image")
all_d0_seu <- NormalizeData(all_d0_seu)
Idents(all_d0_seu) <- "celltype"
markers <- FindAllMarkers(all_d0_seu)
d0_celltype_topmarkers <- markers %>%
  filter(pct.1 > 0.4, avg_log2FC > 0.5) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)
DotPlot(all_d0_seu, group.by="celltype",
        features = unique(d0_celltype_topmarkers$gene, assay='Xenium')) +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1))
ggsave("all_d0_celltypemarkers.pdf", height=length(unique(d0_celltype_topmarkers$cluster))/2,
       width=length(unique(d0_celltype_topmarkers$gene))/4, useDingbats=FALSE)


#Now the day 15 object
all_d15 <- normalizeData(all_d15, method = "LogNorm")
selected_features <- vrFeatures(all_d15)
all_d15 <- getPCA(all_d15, features = selected_features, dims = 30, overwrite = TRUE)
all_d15 <- getUMAP(all_d15, dims = 1:30, overwrite = TRUE)

#Make clusters 
all_d15 <- getProfileNeighbors(all_d15, dims = 1:16, k = 10, method = "SNN", data.type = "pca", graph.key = "SNN")

#We used resolution 1.2 but others could be better.
resols = 1.3
all_d15 <- getClusters(all_d15, resolution = resols, label = paste0("Clusters_", resols), graph = "SNN")
vrEmbeddingPlot(all_d15, group.by = paste0("Clusters_", resols), embedding = "umap", pt.size=0.5)
ggsave(paste0("all_d15_clusters_", resols, ".pdf"))
  

#Calculate and plot cluster markers
all_d15_seu <- VoltRon::as.Seurat(all_d15, cell.assay = "Xenium", type = "image")
all_d15_seu <- NormalizeData(all_d15_seu)
Idents(all_d15_seu) <- paste0("Clusters_", resols)
markers <- FindAllMarkers(all_d15_seu)
topmarkers <- markers %>%
  filter(pct.1 > 0.4, avg_log2FC > 0.5) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)
DotPlot(all_d15_seu, group.by=paste0("Clusters_", resols),
        features = unique(topmarkers$gene, assay='Xenium')) +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1))
ggsave("all_d15_clustermarkers.pdf", height=length(unique(topmarkers$cluster))/2,
       width=length(unique(topmarkers$gene))/4, useDingbats=FALSE)

#Display cell type marker according to the 10x classification from https://www.10xgenomics.com/products/xenium-panels
marker_genes_file <- read.csv("Xenium_mBrain_v1.1_metadata.csv")
marker_genes <- marker_genes_file %>%
  filter(grepl("Astrocyte|Oligodendro|Neuron", Annotation, ignore.case=TRUE, perl=TRUE))
the_features <- as.character(marker_genes$Genes)
names(the_features) <- marker_genes$Annotation
DotPlot(all_d15_seu, group.by=paste0("Clusters_", resols),
        features = the_features) +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1))
ggsave("all_d15_xeniumgenes.pdf", height=length(unique(topmarkers$cluster))/2,
       width=length(unique(marker_genes$Genes))/4, useDingbats=FALSE)

#Display d0 cell type markers
DotPlot(all_d15_seu, group.by=paste0("Clusters_", resols),
        features = unique(d0_celltype_topmarkers$gene, assay='Xenium')) +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1))
ggsave("all_d15_d0_celltype_topmarkers.pdf", height=length(unique(topmarkers$cluster))/2,
       width=length(unique(d0_celltype_topmarkers$gene))/4, useDingbats=FALSE)


#Get localization of clusters
vrSpatialPlot(all_d15, group.by = paste0("Clusters_", resols), alpha = 1, plot.segments = TRUE, background.color = "black")
ggsave("all_d15_Spatial_clusters.pdf", useDingbats=FALSE, height=6, width=6)

#Neurons: Ebf3, Penk, Necab1, Gad1, Gad2
#Pericytes: Acta2, Carmn, Ano1
#Endothelial cells: Cldn5, Adgrl4 but no Acta2
#Oligodendrocytes: Sox10, Opalin, Gjc3
#Astrocytes: Sox9, Rfx4, Wfs1
#Microglia: Trem2, Siglech, Laptm5, Cd53, Cd68
#Ventral leptomeningeal cells: Igfbp6, Col1a1, Dcn (also show strong peripheral localization)
#T cells: Cd3e
all_d15@metadata@cell <- all_d15@metadata@cell %>% mutate(celltype = case_when(!!sym(paste0("Clusters_", resols)) == 1 ~ "ventral leptomeningeal cells",
                                                                               !!sym(paste0("Clusters_", resols)) == 2 ~ "ventral leptomeningeal cells",
                                                                               !!sym(paste0("Clusters_", resols)) == 3 ~ "neurons",
                                                                               !!sym(paste0("Clusters_", resols)) == 4 ~ "neurons",
                                                                               !!sym(paste0("Clusters_", resols)) == 5 ~ "oligodendrocytes",
                                                                               !!sym(paste0("Clusters_", resols)) == 6 ~ "astrocytes",
                                                                               !!sym(paste0("Clusters_", resols)) == 7 ~ "endothelial cells",
                                                                               !!sym(paste0("Clusters_", resols)) == 8 ~ "astrocytes",
                                                                               !!sym(paste0("Clusters_", resols)) == 9 ~ "T cells",
                                                                               !!sym(paste0("Clusters_", resols)) == 10 ~ "microglia",
                                                                               !!sym(paste0("Clusters_", resols)) == 11 ~ "endothelial cells",
                                                                               !!sym(paste0("Clusters_", resols)) == 12 ~ "pericytes",
                                                                               !!sym(paste0("Clusters_", resols)) == 13 ~ "astrocytes",
                                                                               !!sym(paste0("Clusters_", resols)) == 14 ~ "microglia",
                                                                               !!sym(paste0("Clusters_", resols)) == 15 ~ "astrocytes",
                                                                               !!sym(paste0("Clusters_", resols)) == 16 ~ "ventral leptomeningeal cells",
                                                                               !!sym(paste0("Clusters_", resols)) == 17 ~ "microglia",
                                                                               !!sym(paste0("Clusters_", resols)) == 18 ~ "oligodendrocytes",
                                                                               !!sym(paste0("Clusters_", resols)) == 19 ~ "microglia",
                                                                               !!sym(paste0("Clusters_", resols)) == 20 ~ "microglia",
                                                                               !!sym(paste0("Clusters_", resols)) == 21 ~ "neurons",
                                                                               .default = "other"))

saveRDS(all_d15, "./all_d15_mols_celltypeann.rds")


#Some helper plots

#Visualize cell types
vrSpatialPlot(all_d15, colors = the_cols, group.by = "celltype", alpha = 1, plot.segments = TRUE, background.color = "black")
ggsave("all_d15_Spatial_celltypes.pdf", useDingbats=FALSE, height=6, width=6)

#This is to get cluster similarities between d0 and d15
all_d0@metadata@cell$pre_clusters <- paste0("d0_", all_d0@metadata@cell$Clusters_1.1)
all_d15@metadata@cell$pre_clusters <- paste0("d15_", all_d15@metadata@cell$Clusters_1.2)

combined <- merge(all_d0, all_d15)
selected_features <- vrFeatures(combined)
combined <- getPCA(combined, features = selected_features, dims = 30, overwrite = FALSE, type="combPCA")
combined <- getUMAP(combined, dims = 1:30, data.type="combPCA", overwrite = FALSE, umap.key="combUMAP")

combined_seu <- VoltRon::as.Seurat(combined, cell.assay = "Xenium", type = "image")
combined_seu <- NormalizeData(combined_seu)
Idents(combined_seu) <- combined_seu@meta.data$pre_clusters
#combined_seu <- BuildClusterTree(combined_seu, features=rownames(combined_seu))
#PlotClusterTree(combined_seu)
combined_seu <- BuildClusterTree(combined_seu, dims=TRUE, reduction="combUMAP")
PlotClusterTree(combined_seu)



