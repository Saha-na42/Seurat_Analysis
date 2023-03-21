#Loading libraries

install.packages('Seurat')
install.packages('tidyverse')
install.packages("vctrs")
install.packages("ggplot2")
install.packages("hdf5r")

library(Seurat)
library(tidyverse)
library(ggplot2)
library(hdf5r)

#Loading the dataset 

lc<-Read10X_h5(filename = "/Users/sahanabaskar/Desktop/40k_NSCLC_DTC_3p_HT_nextgem_Multiplex_count_raw_feature_bc_matrix.h5")

#Selecting gene expression from the modalities
str(lc)
countslc<-lc$'Gene Expression'

#creating seurat object with non-normalized data

lc.obj<-CreateSeuratObject(counts = countslc, 
                           project = 'NSCLC',
                           min.cells = 3,
                           min.features = 200)
str(lc.obj)


#Quality control
View(lc.obj@meta.data)
#mitochondrial reads
lc.obj[["percent.mt"]]<-PercentageFeatureSet(lc.obj, pattern = '^MT-')
View(lc.obj@meta.data)

VlnPlot(lc.obj, features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol=3)

FeatureScatter(lc.obj, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA') +
  geom_smooth(method = 'lm')

#Filtering
lc.obj<-subset(lc.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Normalization of data
lc.obj<-NormalizeData(lc.obj)

#Highly variable features
lc.obj<-FindVariableFeatures(lc.obj,selection.method = 'vst', nfeatures = 2000)

#top 10 higly variable genes
mostv<-head(VariableFeatures(lc.obj),10)

mostvplot<-VariableFeaturePlot(lc.obj)
LabelPoints(plot = mostvplot,points = mostv, repel = TRUE)

#Scaling
all.genes<-rownames(lc.obj)
lc.obj<-ScaleData(lc.obj, features = all.genes)

#Linear dimensionality reduction
lc.obj<-RunPCA(lc.obj, features = VariableFeatures(object = lc.obj))

print(lc.obj[['pca']], dims = 1:5, nfeatures = 5)
DimHeatmap(lc.obj, dims=1, cells = 500, balanced = TRUE)

#dimensionality of data
ElbowPlot(lc.obj)

#Clustering
lc.obj<-FindNeighbors(lc.obj, dims = 1:15)
lc.obj<-FindClusters(lc.obj, resolution = c(0.1,0.3,0.5,0.7,1))
View(lc.obj@meta.data)

DimPlot(lc.obj, group.by = 'RNA_snn_res.0.1', label = TRUE)

Idents(lc.obj)<-'RNA_snn_res.0.1'
Idents(lc.obj)

#Non-linear dimensionality reduction
lc.obj<-RunUMAP(lc.obj, dims = 1:15)

DimPlot(lc.obj,reduction = 'umap')
