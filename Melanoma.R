#Packages#

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(SingleR)
library(celldex)
library(RColorBrewer)

set.seed(1234)
raw_data <- Read10X(data.dir = "raw_data")

#Creating the seurat object#

Seurat_obj <- CreateSeuratObject(counts = raw_data, project = "B16_Melanoma", min.cell = 3, min.features = 200)

#Adding Meta data#

Seurat_obj$genetype <-"WT"
Seurat_obj$sample <- "sample1"
Seurat_obj
saveRDS(Seurat_obj,"output/seurat_obj_raw.rds")

#Calculate the percetage of mt genes per cell:
Seurat_obj[["percent.mt"]] <-PercentageFeatureSet(Seurat_obj, pattern = "^mt")

#Visualize key QC metrices
VlnPlot(Seurat_obj, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        pt.size = 0.1, 
        ncol = 3)
Seurat_obj <- subset(Seurat_obj, subset = nFeature_RNA > 200 & 
                       nFeature_RNA < 6000 & 
                       percent.mt< 15)

#Normalize the data using log narmalization
Seurat_obj <- NormalizeData(Seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)

#Identify 2000 most variable genes across all cells:
Seurat_obj<- FindVariableFeatures(Seurat_obj, selection.method = "vst", nfeatures = 2000)

#Visualize top 20 variable features
var_plot <- VariableFeaturePlot(Seurat_obj)
LabelPoints(plot = var_plot, points = head(VariableFeatures(Seurat_obj), 20), repel = TRUE)

#Scale data
Seurat_obj <- ScaleData(Seurat_obj)

#Run PCA
Seurat_obj <- RunPCA(Seurat_obj, npcs = 50)

#Visualize the standard deviation of each pc
ElbowPlot(Seurat_obj, ndims = 50)
pcs <- 15

#Identify the k-nearest neighbors
Seurat_obj <- FindNeighbors(Seurat_obj, dims = 1:pcs)

#Visualization
#graph_based clustering
Seurat_obj <- FindClusters(Seurat_obj, resolution = 0.6)

#UMAP for 2D visualization
Seurat_obj <- RunUMAP(Seurat_obj, dims = 1:pcs) 

#Plot UMAP clustering
DimPlot(Seurat_obj, reduction = "umap", label = TRUE, repel = TRUE) + ggtitle("UMAP: clustered cells")

#Annotation
ref <- celldex::MouseRNAseqData()
annotations <- SingleR(test = GetAssayData(Seurat_obj, layer = "data"), ref = ref, labels = ref$label.main)

Seurat_obj$celltype <- annotations$labels
DimPlot(Seurat_obj, group.by = "celltype", label = TRUE, repel = TRUE) + ggtitle("UMAP: Annotated cell types")

#Find Marker genes for all clusters
markers <- FindAllMarkers(Seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#Extract top 5
top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(n=5, order_by = avg_log2FC)

#save markers
write.csv(top_markers, "cluster_markers.csv")


cluster0_vs_cluster1_DEGs.csv

#differential expression between two cluster
de_genes <- FindMarkers(Seurat_obj, ident.1 = 0,
                        ident.2 = 1, min.pct = 0.25, logfc.threshold = 0.25)
#save
write.csv(de_genes, "cluster0_vs_cluster1_DEGs.csv")

#view top DEgenes
head(de_genes)


de_genes <- read.csv("cluster0_vs_cluster1_DEGs.csv")

#Add a colomn for significance
de_genes$gene <- rownames(de_genes)
de_genes$significant <- ifelse(de_genes$p_val_adj < 0.05 & abs(de_genes$avg_log2FC) > 0.5, "Yes", "No")

#Plot using ggplot2
library(ggplot2)
ggplot(de_genes, aes(x = avg_log2FC, y = -log10(p_val_adj),
color = significant)) + geom_point(alpha = 0.8) + scale_color_manual(values = c("grey", "red")) + 
theme_minimal() + labs(title = "volcano plot: cluster 0 vs cluster 1", 
                       x = "Log2 Fold change", y = "-Log10 Adjusted p_value")


