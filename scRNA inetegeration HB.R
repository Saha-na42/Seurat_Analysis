#Loading Libraries
library(tidyverse)
library(Seurat)
library(ggplot2)
library(gridExtra)
library(tidyr)

#Loading data
direc<-list.dirs(path = "/Users/sahanabaskar/Desktop/data",recursive = FALSE,full.names = F)
direc

for(x in direc){
  name<- gsub('_filtered_feature_bc_matrix',"",x)
  
  ct<-ReadMtx(mtx = paste0('/Users/sahanabaskar/Desktop/data/',x,'/matrix.mtx'),
              features = paste0('/Users/sahanabaskar/Desktop/data/',x,'/features.tsv'),
              cells = paste0('/Users/sahanabaskar/Desktop/data/',x,'/barcodes.tsv'))
  
  assign(name,CreateSeuratObject(counts = ct))}

#merging dataset

merged_obj<-merge(HB17_background, y = c(HB17_PDX,HB17_tumor,HB30_PDX,HB30_tumor,
                                         HB53_tumor,HB53_background),
                  add.cell.ids = ls()[3:9],
                  project = 'Hepatoblastoma')

merged_obj

adj.matrix <- NULL

#Quality Control
View(merged_obj@meta.data)

dim(merged_obj@meta.data)
summary(merged_obj@meta.data$nCount_RNA)
summary(merged_obj@meta.data$nFeature_RNA)

#creating a sample column
merged_obj$sample<-rownames(merged_obj@meta.data)

merged_obj@meta.data<-separate(merged_obj@meta.data, col ="sample",
                               into = c('Patient','type','barcode'),
                               sep = '_')
#mito percentage
merged_obj$mitoperc<-PercentageFeatureSet(merged_obj,pattern = '^MT-')

#merged_obj$rb_perc<-PercentageFeatureSet(merged_obj,pattern = "^RP[SL]")

View(merged_obj@meta.data)

VlnPlot(merged_obj, features = c("nFeature_RNA","nCount_RNA","mitoperc"),ncol = 4,pt.size = 0.1) & 
  theme(plot.title = element_text(size=10))

#filtering

filtered_obj<-subset(merged_obj, subset = nCount_RNA > 800 & nFeature_RNA > 500 &
                       mitoperc < 5)
filtered_obj

#Normalization
filtered_obj<-NormalizeData(object = filtered_obj)

#highly variable genes
filtered_obj<-FindVariableFeatures(object = filtered_obj)

#top 10 higly variable genes
hv<-head(VariableFeatures(filtered_obj),10)

hvplot<-VariableFeaturePlot(filtered_obj)
LabelPoints(plot = hvplot,points = hv, repel = TRUE)

#Scaling
filtered_obj<-ScaleData(object = filtered_obj)

#Linear Dimensionality Reduction
filtered_obj<-RunPCA(object = filtered_obj)

print(filtered_obj[['pca']], dims = 1:5, nfeatures = 5)
DimHeatmap(filtered_obj, dims=1, cells = 500, balanced = TRUE)

ElbowPlot(filtered_obj)

#Clustering
filtered_obj<-FindNeighbors(object = filtered_obj, dims = 1:20)
filtered_obj<-FindClusters(filtered_obj, resolution = c(0.1,0.3,0.5,0.7,1))
View(filtered.obj@meta.data)

DimPlot(filtered_obj, group.by = 'RNA_snn_res.0.1', label = TRUE)

Idents(filtered_obj)<-'RNA_snn_res.0.1'
Idents(filtered_obj)

filtered_obj<-RunUMAP(object = filtered_obj,dims =1:20)

#visualization
plt1<-DimPlot(filtered_obj,reduction = "umap", group.by = 'Patient')
plt2<-DimPlot(filtered_obj,reduction = "umap", group.by = 'type',
        cols = c('red','green','blue'))

grid.arrange(plt1,plt2,ncol = 2, nrow = 2)

#integration to correct batch effects
obj.list<-SplitObject(filtered_obj,split.by = 'Patient')
for(i in 1:length(obj.list)){
  obj.list[[i]] <-NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <-FindVariableFeatures
  
}


#integration features
features<-SelectIntegrationFeatures(object,list = obj.list)

#integration anchors(CCA)
anchors_obj<-FindIntegrationAnchors(object.list = obj.list,
                                anchor.features = features)

#integrate data
int<-IntegrateData(anchorset = anchors_obj)

#Scaling,pca & Umap on integrated data
int<-ScaleData(object = int)
int<-RunPCA(object = int)
int<-RunUMAP(object = int,dims = 1:50)

p3<-DimPlot(int,reduction = "umap",group.by ='Patient')
p4<-DimPlot(int,reduction = "umap",group.by ='type')

grid.arrange(p3,p4, ncol =2)