library(dplyr)
library(tibble)
library(ggplot2)
library(Seurat)


#Otud7b expression in multiple sclerosis
#Get sNuc data from GSE180759 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE180759)
#Plot as below, rearrange in Adobe Illustrator
setwd("/path/to/GSE180759/")  
GSE180759_mtx <- read.table("GSE180759_expression_matrix.csv.gz", sep=",")
GSE180759_meta <- read.table("GSE180759_annotation.txt.gz", sep="\t", header=TRUE)

seu = CreateSeuratObject(counts = GSE180759_mtx, meta.data = GSE180759_meta, min.cells = 50, min.features = 3)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
seu <- ScaleData(seu)
seu <- RunPCA(seu)
ElbowPlot(seu)
seu <- RunUMAP(seu, dims = 1:19)
UMAPPlot(seu)
seu <- FindNeighbors(seu, dims = 1:15)
seu <- FindClusters(seu, resolution = 0.8)
DimPlot(seu, reduction = "umap", label=TRUE)+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  ggtitle("Cells colored by cluster")
Idents(seu) <- seu@meta.data$cell_type

seu@meta.data$celltype_pathology <- paste(seu@meta.data$cell_type, seu@meta.data$pathology, sep="|")

DotPlot(seu, group.by='celltype_pathology', cols='RdYlBu', idents = c("astrocytes"),
        features = "OTUD7B")
ggsave("MS_sNucSeq_Dotplot_2.pdf", width=6.6, height=1.8)
