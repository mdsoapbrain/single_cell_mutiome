library(Seurat)
library(WGCNA)
library(igraph)
library(devtools)
library(sctransform)
library(dplyr)
library(scDblFinder)
library(SingleCellExperiment)
library(Signac)
library(GenomicRanges)
library(rtracklayer)
library(EnsDb.Rnorvegicus.v105)
library(scales)
library(future) 
library(harmony)
library(ggplot2)


rm(list = ls())
data<-readRDS("mutiome_final.rds")

table(data$phenotype)
table(data$predict.id)

library(dplyr)
library(ggplot2)

library(dplyr)
library(ggplot2)

# 从Seurat对象中提取元数据
metadata <- data@meta.data

# 计算每个predict.id在每个phenotype中的数量
count_data <- metadata %>% 
  group_by(phenotype, predict.id) %>% 
  summarise(count = n())

# 计算每个predict.id在每个phenotype中的比例
count_data <- count_data %>% 
  group_by(phenotype) %>% 
  mutate(proportion = count / sum(count))

# 绘制累积条形图
ggplot(count_data, aes(x = phenotype, y = proportion, fill = predict.id)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(y = "Proportion", x = "Phenotype", fill = "Celltypes") +
  theme_minimal()

write.csv(count_data, 'count_data.csv')






