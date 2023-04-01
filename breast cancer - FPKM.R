# Manipulation of FPKM normalized gene data

#Loading libraries
library(dplyr)
library(tidyverse)
library(GEOquery)

# Load the dataset
gene_data<-read.csv(file = "/Users/sahanabaskar/Documents/datasets/GSE183947_fpkm.csv")

# Get the metadata
gse<-getGEO(GEO = 'GSE183947',GSEMatrix = TRUE)
gse
metadata<-pData(phenoData(gse[[1]]))
head(metadata)

# Selecting the columns
metadata_s<- metadata %>%
  select(c(1,10,11,17)) %>% rename(tissuesample = characteristics_ch1,
                                   metastasis = characteristics_ch1.1) %>%
  mutate(tissuesample  = gsub("tissue:", "",tissuesample)) %>%
  mutate(metastasis = gsub("metastasis:", "", metastasis))

head(metadata_s)

# Modifying gene_data

gene_data<-gene_data %>% rename( Gene = X) %>%
  gather(., key = "samples", value = "FPKM", -Gene) 

head(gene_data)

#Combining the columns from metadata_s with gene_data

gene_data<- gene_data %>% left_join(., metadata_s,
          by = c("samples" = "description"))

head(gene_data)

# Grouping by BRCA1 & BRCA2

data_N <- gene_data %>% 
          filter(Gene == "BRCA1" | Gene == "BRCA2") %>%
          group_by(Gene,tissuesample)  %>%
             summarise(mean_FPKM = mean(FPKM))

head(data_N)

# Visualization of the data

gene_int<-gene_data %>%
  group_by(Gene) %>%
  summarise(mean_fpkm = mean(FPKM)) 

gene_int

# brca1 & brca2 - bar plot
gene_data %>%
    filter(Gene == "BRCA1" | Gene == "BRCA2") %>%
    ggplot(., aes(x = samples, y = FPKM, fill = Gene)) +
                          geom_col()

gene_data %>%
  filter(Gene == "BRCA2") %>%
  ggplot(., aes(x = samples, y = FPKM, fill = tissuesample)) +
  geom_col()

# density plot to visualise the distribution of FPKM

gene_data %>%
  filter(Gene == "BRCA1" | Gene == "BRCA2" ) %>%
  ggplot(., aes(x = FPKM, fill = tissuesample)) +
  geom_density(alpha = 0.5)

# boxplot - metastasis distribution

gene_data %>%
  filter(Gene == "BRCA1" | Gene == "BRCA2" ) %>%
  ggplot(., aes(x = metastasis, y = FPKM, fill = Gene)) +
  geom_boxplot()

# Scatterplot for BRCA1 & BRCA2
gene_data %>%
  filter(Gene == "BRCA1" | Gene == "BRCA2" ) %>%
  spread(key = Gene, value = FPKM) %>%
  ggplot(., aes(x = BRCA1, y = BRCA2, color = tissuesample)) +
  geom_point() +
  geom_smooth(method = "lm",se = TRUE)

#Heatmap for genes of interest

gene.interest<- c("BRCA1","BRCA2","TP53","MYCN","ALK")

gene_data %>%
  filter(Gene %in% gene.interest) %>%
  ggplot(., aes(x = samples, y = Gene, fill = FPKM)) +
  geom_tile()+
  scale_fill_gradient(low = 'white', high = 'darkgreen')
