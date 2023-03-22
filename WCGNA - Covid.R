#WCGNA
install.packages("BiocManager")
BiocManager::install("WGCNA")
BiocManager::install("DESeq2")
BiocManager::install("BiocParallel")
BiocManager::install("GEOquery")
install.packages("CorLevelPlot")
install.packages("remotes")
remotes::install_github("kevinblighe/CorLevelPlot")


#loading libraries

library(WGCNA)
library(BiocParallel)
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)
library(dplyr)

#calculation of correlations in the presence of missing data

enableWGCNAThreads() 

#loading data
data<-read.delim('/Users/sahanabaskar/Downloads/GSE152418_p20047_Study1_RawCounts.txt',header = TRUE)

#metadata
geo_id<-"GSE152418"
gse<-getGEO(geo_id, GSEMatrix = TRUE)
phenoData<-pData(phenoData(gse[[1]]))

head(phenoData)
phenoData<-phenoData[,c(1,2,46:50)]

#expression data
data[1:10,1:10]

data<-data %>%
  gather(key = 'samples', value = 'counts', -ENSEMBLID) %>%
           mutate(samples = gsub('\\.','_',samples)) %>%
           inner_join(.,phenoData,by = c('samples'  = 'title')) %>%
           select(1,3,4) %>%
           spread(key = "geo_accession", value = "counts") %>%
           column_to_rownames(var = 'ENSEMBLID')

#Quality control

#outlier detection

gd<-goodSamplesGenes(t(data))
summary(gd)
gd$allOK

table(gd$goodGenes)
table(gd$goodSamples)

data<-data[gd$goodGenes == TRUE,]

#Hhierachial clustering for outliers

htree<-hclust(dist(t(data)), method ="average")
plot(htree)

#pca - for outlier detection

pca<-prcomp(t(data))
pca1<-pca$x

pca.var<-pca$sdev^2
pca.var.perc<-round(pca.var/sum(pca.var)*100,digits = 2)

pca1<- as.data.frame(pca1)

ggplot(pca1, aes(PC1,PC2))+
  geom_point() +
  geom_text(label = rownames(pca1)) +
  labs(x = paste0("PC1:",pca.var.perc[1],'%'),
       y = paste0("PC2:",pca.var.perc[2],'%'))

#excluding samples

samples.excluded<-c('GSM4614993','GSM4614995','GSM4615000')
data.subset <- data[,!(colnames(data) %in% samples.excluded)]       

#Normalization

#rempve outliers from phenodaata

colData<-phenoData %>%
  subset(!row.names(.) %in% samples.excluded)

names(colData)
names(colData) <-gsub(':ch1',"",names(colData))
names(colData) <-gsub('\\s',"_",names(colData))

#making the rownames similar to colnames

all(rownames(colData) %in% colnames(data.subset))
all(rownames(colData) == colnames(data.subset))


cat("Sample names in colData:\n")
cat(rownames(colData), sep = "\n")

cat("\nSample names in data.subset:\n")
cat(colnames(data.subset), sep = "\n")



#counts < 15 , > 75%  of samples (31*0.75 = 23.25)

ds <- DESeqDataSetFromMatrix(countData = data.subset,
                             colData = colData,
                             design = ~ 1) 





dds75 <- dds[rowSums(counts(dds) >= 15) >= 24,]
nrow(dds75) 


# perform variance stabilization
dds_norm <- vst(dds75)


# get normalized counts
norm.counts <- assay(dds_norm) %>% 
  t()

#network construction 
# Choose a set of soft-thresholding powers
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# Call the network topology analysis function
sft <- pickSoftThreshold(norm.counts,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)


sft.data <- sft$fitIndices

# visualization to pick power

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()


a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()


grid.arrange(a1, a2, nrow = 2)


# convert matrix to numeric
norm.counts[] <- sapply(norm.counts, as.numeric)

soft_power <- 18
temp_cor <- cor
cor <- WGCNA::cor


# memory estimate w.r.t blocksize
bwnet <- blockwiseModules(norm.counts,
                          maxBlockSize = 14000,
                          TOMType = "signed",
                          power = soft_power,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3)


cor <- temp_cor


# 5. Module Eigengenes ---------------------------------------------------------
module_eigengenes <- bwnet$MEs


# Print out a preview
head(module_eigengenes)


# get number of genes for each module
table(bwnet$colors)

# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)




# grey module = all genes that doesn't fall into other modules were assigned to the grey module


# module trait associations



# create traits file - binarize categorical variables
traits <- colData %>% 
  mutate(disease_state_bin = ifelse(grepl('COVID', disease_state), 1, 0)) %>% 
  select(8)


# binarize categorical variables

colData$severity <- factor(colData$severity, levels = c("Healthy", "Convalescent", "ICU", "Moderate", "Severe"))

severity.out <- binarizeCategoricalColumns(colData$severity,
                                           includePairwise = FALSE,
                                           includeLevelVsAll = TRUE,
                                           minCount = 1)


traits <- cbind(traits, severity.out)


# Define numbers of genes and samples
nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)


module.trait.corr <- cor(module_eigengenes, traits, use = 'p')
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)



# module-trait association as a heatmap

heatmap.data <- merge(module_eigengenes, traits, by = 'row.names')

head(heatmap.data)

heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = 'Row.names')




CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[18:22],
             y = names(heatmap.data)[1:17],
             col = c("blue1", "skyblue", "white", "pink", "red"))



module.gene.mapping <- as.data.frame(bwnet$colors)
module.gene.mapping %>% 
  filter(`bwnet$colors` == 'turquoise') %>% 
  rownames()



# 6B. Intramodular analysis: Identifying driver genes ---------------



# Calculate the module membership and the associated p-values

module.membership.measure <- cor(module_eigengenes, norm.counts, use = 'p')
module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nSamples)


module.membership.measure.pvals[1:10,1:10]


# Calculate the gene significance and associated p-values

gene.signf.corr <- cor(norm.counts, traits$data.Severe.vs.all, use = 'p')
gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples)


gene.signf.corr.pvals %>% 
  as.data.frame() %>% 
  arrange(V1) %>% 
  head(25)






