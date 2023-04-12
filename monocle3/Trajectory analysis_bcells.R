#setwd()

#Loading libraries

install.packages(c("remotes", "devtools", "BiocManager"))
BiocManager::install(c("SingleCellExperiment", "DelayedArray", "DelayedMatrixStats"))

devtools::install_version('speedglm', '0.3-4', repos = 'https://packagemanager.rstudio.com/cran/2023-03-31')

devtools::install_github("cole-trapnell-lab/monocle3")

remotes::install_github("satijalab/seurat-wrappers") 

library(monocle3)
library(SeuratWrappers)
library(Seurat)
library(ggplot2)
library(tidyverse)

#Reading data
markers<- read.delim('ABC_Marker.txt',header = TRUE)
metadata<-read.delim('ABC_Meta.txt',header = TRUE)
matrix<-read.delim('ABC_umi_matrix_7551_cells.csv',header = TRUE,sep=',')

#Creating Seuratobj

expr_t<-t(matrix)

exp_obj<-CreateSeuratObject(counts = expr_t)
View(exp_obj@meta.data)

exp_obj@meta.data<-merge(exp_obj@meta.data,metadata, by.x='row.names',by.y='cell_id')
exp_obj@meta.data<- exp_obj@meta.data %>%
  column_to_rownames(var = 'Row.names')

exp_obj$mito_perc<-PercentageFeatureSet(exp_obj,pattern = '^MT-')

#Filtering
exp_obj.filtered<-subset(exp_obj, subset= nCount_RNA > 800 &
         nFeature_RNA > 500 &
         mito_perc < 10)

# Selecting B cells

View(exp_obj.filtered@meta.data)
unique(exp_obj.filtered@meta.data$population)

Idents(exp_obj.filtered)<-unique(exp_obj.filtered$population)
b_cell.obj <- subset(exp_obj.filtered, population == "b")
b_cell.obj
unique(b_cell.obj@meta.data$redefined_cluster)
View(b.cell.obj@meta.data)

#Seurat Pre-processing

b_cell.obj<-NormalizeData(b_cell.obj)
b_cell.obj<-FindVariableFeatures(b_cell.obj)
b_cell.obj<-ScaleData(b_cell.obj)
b_cell.obj<-RunPCA(b_cell.obj)
b_cell.obj <- FindNeighbors(b_cell.obj, dims = 1:30, graph.name = "b_cell_graph")
b_cell.obj <- FindClusters(b_cell.obj, resolution = 0.9, graph.name = "b_cell_graph")
b_cell.obj<-RunUMAP(b_cell.obj, dims = 1:30, n.neighbors = 50)

b_plot1<-DimPlot(b_cell.obj, reduction = 'umap', group.by = 'redefined_cluster',
                 label = TRUE)
b_plot2<-DimPlot(b_cell.obj, reduction = 'umap', group.by = 'seurat_clusters',
                 label = TRUE)
b_plot1|b_plot2

#Converting seurat obj into cell_data_Set object

b.cds<-as.cell_data_set(b_cell.obj)

colData(b.cds)#cell metadata
fData(b.cds)#gene metadata
rownames(fData(b.cds))[1:10]
fData(b.cds)$gene_shortname<-rownames(fData(b.cds))
counts(b.cds)#count matrix

#assigning partitions
recreate.partition<-c(rep(1,length(b.cds@colData@rownames)))
names(recreate.partition)<-b.cds@colData@rownames
recreate.partition<-as.factor(recreate.partition)

b.cds@clusters$UMAP$partitions<-recreate.partition

#assigning cluster info
list_cluster<-b_cell.obj@active.ident
b.cds@clusters$UMAP$clusters<-list_cluster

#assigning umap coordinated to cell embeddings
b.cds@int_colData@listData$reducedDims$UMAP <-b_cell.obj@reductions$umap@cell.embeddings

before.trajectory<-plot_cells(b.cds,
           color_cells_by = 'cluster',
           label_groups_by_cluster = FALSE,
           group_label_size = 5) +
  theme(legend.position = 'right')

cluster.names<-plot_cells(b.cds,
           color_cells_by = 'redefined_cluster',
           label_groups_by_cluster = FALSE,
           group_label_size = 5) +
  scale_color_manual(values = c('red','blue','green','yellow','grey','cyan','orange')) +
  theme(legend.position = 'right')

before.trajectory|cluster.names

#trajectory graph
b.cds<-learn_graph(b.cds, use_partition = FALSE)

plot_cells(b.cds,
           color_cells_by = 'redefined_cluster',
           label_groups_by_cluster = F,
           label_branch_points =  F,
           label_roots = F,
           label_leaves = F,
           group_label_size = 5)

#ordering the cells by pseudotime
b.cds<-order_cells(b.cds, reduction_method = 'UMAP', root_cells = colnames(b.cds[,clusters(b.cds) == 5]))

plot_cells(b.cds,
           color_cells_by = 'pseudotime',
           label_branch_points = F,
           label_roots = F,
           label_leaves = F)
#cells ordered by monocle3 pseudotime

pseudotime(b.cds)
b.cds$monocle3_pseudotime<-pseudotime(b.cds)

data.pseudo<-as.data.frame(colData(b.cds))

ggplot(data.pseudo,aes(monocle3_pseudotime,redefined_cluster,
                       reorder(redefined_cluster,monocle3_pseudotime,median), 
                       fill = redefined_cluster))+
  geom_boxplot()

#identifying genes that change as a function of pseudotime

deg_bcells<-graph_test(b.cds,neighbor_graph = 'principal_graph',cores = 4)

deg_bcells %>%
  arrange(q_value) %>%
  filter(status == 'OK') %>%
  head()

FeaturePlot(b_cell.obj,features = c('STMN1','CD52','HMGN2'))

b_cell.obj$pseudotime<-pseudotime(b.cds)
Idents(b_cell.obj) <-b_cell.obj$redefined_cluster
FeaturePlot(b_cell.obj, features = 'pseudotime',label = TRUE)
