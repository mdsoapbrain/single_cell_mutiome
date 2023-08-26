setwd("/Volumes/Seagate/multiome_cp")
# single-cell analysis package
library(Seurat)

# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)

# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)

# gene enrichment packages
library(enrichR)
library(GeneOverlap)

# packages for TF motif analysis
library(JASPAR2020)
library(motifmatchr)
library(TFBSTools)
library(EnsDb.Rnorvegicus.v105)
library(GenomicRanges)
# pseudotime
library(monocle3)
library(SeuratWrappers)

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)

# optionally enable multithreading
#enableWGCNAThreads(nThreads = 8)

# load snRNA-seq dataset
data<-readRDS("mutiome_final.rds")
Drd2_data<-subset(x= data, subset = predict.id == 'Drd2-MSN')
Olig_data<-subset(x= data, subset = predict.id == 'Olig')
Drd1_data<-subset(x= data, subset = predict.id == 'Drd1-MSN')
Astrocyte_data<-subset(x= data, subset = predict.id == 'Astrocyte')
Microglia_data<-subset(x= data, subset = predict.id == 'Microglia')
OPC_data<-subset(x= data, subset = predict.id == 'OPC')
Glutamatergic_data<-subset(x= data, subset = predict.id == 'Glutamatergic')
Interneuron_data<-subset(x= data, subset = predict.id == 'Interneuron')
###########Drd case control#######
Drd2_ROT<-subset(x= Drd2_data, subset = phenotype == 'ROT')
Drd2_control<-subset(x= Drd2_data, subset = phenotype == 'control')

########### pre-processing using seurat#########
Drd2_ROT <- NormalizeData(Drd2_ROT)
Drd2_ROT <- FindVariableFeatures(Drd2_ROT)
Drd2_ROT <- ScaleData(Drd2_ROT)
Drd2_ROT <- RunPCA(Drd2_ROT)
Drd2_ROT <- FindNeighbors(Drd2_ROT, dims = 1:30)
Drd2_ROT <- FindClusters(Drd2_ROT, resolution = 0.9)
Drd2_ROT <- RunUMAP(Drd2_ROT, dims = 1:30, n.neighbors = 50)
DimPlot(Drd2_ROT, reduction = 'umap', group.by = 'seurat_clusters', label = T)

library(monocle3)
library(SeuratWrappers)

# convert the seurat object to CDS
cds <- as.cell_data_set(Drd2_ROT)
# assign paritions
reacreate.partition <- c(rep(1,length(cds@colData@rownames)))
names(reacreate.partition) <- cds@colData@rownames
reacreate.partition <- as.factor(reacreate.partition)
cds@clusters$UMAP$partitions <- reacreate.partition 
list_cluster <- Drd2_ROT@active.ident
cds@clusters$UMAP$clusters <- list_cluster
cds@int_colData@listData$reducedDims$UMAP <- Drd2_ROT@reductions$wnn.umap@cell.embeddings
cds <- cluster_cells(cds, reduction_method='UMAP')

# learn graph for pseudotime
cds <- learn_graph(cds)

# plot the pseudotime graph:
p1 <- plot_cells(
  cds = cds,
  color_cells_by = "seurat_clusters",
  show_trajectory_graph = TRUE,
  label_principal_points = TRUE
) 

# plot the UMAP partitions from the clustering algorithm
p2 <- plot_cells(
  cds = cds,
  color_cells_by = "partition",
  show_trajectory_graph = FALSE
)

pdf(('umap_monocle_drd2_rot.pdf'),  width=4, height=6)
p1
dev.off()

# get principal node & order cells
principal_node <- 'Y_7'
cds <- order_cells(cds,root_pr_nodes = principal_node)

# add pseudotime to seurat object:
Drd2_ROT$pseudotime <- pseudotime(cds)

Drd2_ROT$UMAP1 <- Drd2_ROT@reductions$wnn.umap@cell.embeddings[,1]
Drd2_ROT$UMAP2 <- Drd2_ROT@reductions$wnn.umap@cell.embeddings[,2]
library(viridis)
Drd2_ROT$Drd2_ROT <- ifelse(Drd2_ROT$predict.id %in% c('Drd2-MSN'),
                                    Drd2_ROT$pseudotime, NA)
p1 <- Drd2_ROT@meta.data %>%
  ggplot(aes(x=UMAP1, y=UMAP2, color=Drd2_ROT)) +
  ggrastr::rasterise(geom_point(size=1), dpi=500, scale=0.75) +
  coord_equal() +
  scale_color_gradientn(colors=plasma(256), na.value='grey') +
  umap_theme()

pdf(('pseudotime_monocle_drd2_rot.pdf'),  width=4, height=6)
p1
dev.off()









###########

Drd2_data <- NormalizeData(Drd2_data)
Drd2_data <- FindVariableFeatures(Drd2_data)
Drd2_data <- ScaleData(Drd2_data)
Drd2_data <- RunPCA(Drd2_data)
Drd2_data <- FindNeighbors(Drd2_data, dims = 1:30)
Drd2_data <- FindClusters(Drd2_data, resolution = 0.9)
Drd2_data <- RunUMAP(Drd2_data, dims = 1:30, n.neighbors = 50)
DimPlot(Drd2_data, reduction = 'umap', group.by = 'seurat_clusters', label = T)

# convert the seurat object to CDS
cds2 <- as.cell_data_set(Drd2_data)
# assign paritions
reacreate.partition <- c(rep(1,length(cds2@colData@rownames)))
names(reacreate.partition) <- cds2@colData@rownames
reacreate.partition <- as.factor(reacreate.partition)
cds2@clusters$UMAP$partitions <- reacreate.partition 
list_cluster <- Drd2_data@active.ident
cds2@clusters$UMAP$clusters <- list_cluster
cds2@int_colData@listData$reducedDims$UMAP <- Drd2_data@reductions$wnn.umap@cell.embeddings
cds2 <- cluster_cells(cds2, reduction_method='UMAP')
# learn graph for pseudotime
cds2 <- learn_graph(cds2)

# plot the pseudotime graph:
p1 <- plot_cells(
  cds = cds2,
  color_cells_by = "seurat_clusters",
  show_trajectory_graph = TRUE,
  label_principal_points = TRUE
) 

pdf(('umap_monocle_drd2.pdf'),  width=4, height=6)
p1
dev.off()

# get principal node & order cells
principal_node <- 'Y_12'
cds2 <- order_cells(cds2,root_pr_nodes = principal_node)

# add pseudotime to seurat object:
Drd2_data$pseudotime <- pseudotime(cds2)

Drd2_data$UMAP1 <- Drd2_data@reductions$wnn.umap@cell.embeddings[,1]
Drd2_data$UMAP2 <- Drd2_data@reductions$wnn.umap@cell.embeddings[,2]
library(viridis)
Drd2_data$Drd2_ROT <- ifelse(Drd2_data$phenotype %in% c('ROT'),
                              Drd2_data$pseudotime, NA)
p1 <- Drd2_data@meta.data %>%
  ggplot(aes(x=UMAP1, y=UMAP2, color=Drd2_ROT)) +
  ggrastr::rasterise(geom_point(size=1), dpi=500, scale=0.75) +
  coord_equal() +
  scale_color_gradientn(colors=viridis(256), na.value='grey') +
  umap_theme()

pdf(('pseudotime_monocle_drd2_control.pdf'),  width=4, height=6)
p1
dev.off()

Drd2_data$Drd2_control <- ifelse(Drd2_data$phenotype %in% c('control'),
                             Drd2_data$pseudotime, NA)
p2 <- Drd2_data@meta.data %>%
  ggplot(aes(x=UMAP1, y=UMAP2, color=Drd2_control)) +
  ggrastr::rasterise(geom_point(size=1), dpi=500, scale=0.75) +
  coord_equal() +
  scale_color_gradientn(colors=viridis(256), na.value='grey') +
  umap_theme()

pdf(('pseudotime_monocle_drd2.pdf'),  width=12, height=6)
p1+p2
dev.off()











###########

Drd1_data <- NormalizeData(Drd1_data)
Drd1_data <- FindVariableFeatures(Drd1_data)
Drd1_data <- ScaleData(Drd1_data)
Drd1_data <- RunPCA(Drd1_data)
Drd1_data <- FindNeighbors(Drd1_data, dims = 1:30)
Drd1_data <- FindClusters(Drd1_data, resolution = 0.9)
Drd1_data <- RunUMAP(Drd1_data, dims = 1:30, n.neighbors = 50)
DimPlot(Drd1_data, reduction = 'umap', group.by = 'seurat_clusters', label = T)

# convert the seurat object to CDS
cds3 <- as.cell_data_set(Drd1_data)
# assign paritions
reacreate.partition <- c(rep(1,length(cds3@colData@rownames)))
names(reacreate.partition) <- cds3@colData@rownames
reacreate.partition <- as.factor(reacreate.partition)
cds3@clusters$UMAP$partitions <- reacreate.partition 
list_cluster <- Drd1_data@active.ident
cds3@clusters$UMAP$clusters <- list_cluster
cds3@int_colData@listData$reducedDims$UMAP <- Drd1_data@reductions$wnn.umap@cell.embeddings
cds3 <- cluster_cells(cds3, reduction_method='UMAP')
# learn graph for pseudotime
cds3 <- learn_graph(cds3)

# plot the pseudotime graph:
p1 <- plot_cells(
  cds = cds3,
  color_cells_by = "seurat_clusters",
  show_trajectory_graph = TRUE,
  label_principal_points = TRUE
) 

pdf(('umap_monocle_drd1.pdf'),  width=4, height=6)
p1
dev.off()

# get principal node & order cells
principal_node <- 'Y_18'
cds3 <- order_cells(cds3,root_pr_nodes = principal_node)

# add pseudotime to seurat object:
Drd1_data$pseudotime <- pseudotime(cds3)

Drd1_data$UMAP1 <- Drd1_data@reductions$wnn.umap@cell.embeddings[,1]
Drd1_data$UMAP2 <- Drd1_data@reductions$wnn.umap@cell.embeddings[,2]
library(viridis)
Drd1_data$Drd1_ROT <- ifelse(Drd1_data$phenotype %in% c('ROT'),
                             Drd1_data$pseudotime, NA)
p1 <- Drd1_data@meta.data %>%
  ggplot(aes(x=UMAP1, y=UMAP2, color=Drd1_ROT)) +
  ggrastr::rasterise(geom_point(size=1), dpi=500, scale=0.75) +
  coord_equal() +
  scale_color_gradientn(colors=plasma(256), na.value='grey') +
  umap_theme()

pdf(('pseudotime_monocle_drd1_ROT.pdf'),  width=4, height=6)
p1
dev.off()

Drd1_data$Drd1_control <- ifelse(Drd1_data$phenotype %in% c('control'),
                                 Drd1_data$pseudotime, NA)
p2 <- Drd1_data@meta.data %>%
  ggplot(aes(x=UMAP1, y=UMAP2, color=Drd1_control)) +
  ggrastr::rasterise(geom_point(size=1), dpi=500, scale=0.75) +
  coord_equal() +
  scale_color_gradientn(colors=plasma(256), na.value='grey') +
  umap_theme()

pdf(('pseudotime_monocle_drd1_control.pdf'),  width=12, height=6)
p1+p2
dev.off()

###########Drd2#############
DefaultAssay(Drd2_data)<-'RNA'
Drd2_data <- SetupForWGCNA(
  Drd2_data,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "Drd2_data" # the name of the hdWGCNA experiment
)

# construct metacells  in each group
Drd2_data <- MetacellsByGroups(
  seurat_obj = Drd2_data,
  group.by = c("predict.id", "phenotype"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'wnn.umap', # select the dimensionality reduction to perform KNN on
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'predict.id' # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
Drd2_data <- NormalizeMetacells(Drd2_data)
Drd2_data <- SetDatExpr(
  Drd2_data,
  group_name = "Drd2-MSN", # the name of the group of interest in the group.by column
  group.by='predict.id', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  slot = 'data' # using normalized data
)
# Test different soft powers:
Drd2_data <- TestSoftPowers(
  Drd2_data,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(Drd2_data)

# assemble with patchwork
wrap_plots(plot_list, ncol=2)
# construct co-expression network:
Drd2_data <- ConstructNetwork(
  Drd2_data, soft_power=5,
  setDatExpr=FALSE,overwrite_tom = TRUE,
  tom_name = 'Drd2-MSN' # name of the topoligical overlap matrix written to disk
)
PlotDendrogram(Drd2_data, main='Drd2-MSN hdWGCNA Dendrogram')

# need to run ScaleData first or else harmony throws an error:
Drd2_data <- ScaleData(Drd2_data, features=VariableFeatures(Drd2_data))

# compute all MEs in the full single-cell dataset
Drd2_data <- ModuleEigengenes(
  Drd2_data,
  group.by.vars="phenotype"
)
# harmonized module eigengenes:
hMEs <- GetMEs(Drd2_data)
# module eigengenes:
MEs <- GetMEs(Drd2_data, harmonized=FALSE)

# compute eigengene-based connectivity (kME):
Drd2_data <- ModuleConnectivity(
  Drd2_data,
  group.by = 'predict.id', group_name = 'Drd2-MSN'
)
# rename the modules
Drd2_data <- ResetModuleNames(
  Drd2_data,
  new_name = "Drd2-MSN"
)
# plot genes ranked by kME for each module

hub_df <-GetHubGenes(Drd2_data, n_hubs = 10)

pdf("Drd2_hubgebenw.pdf", width=6, height=8)
HubGeneNetworkPlot(
  Drd2_data,
  n_hubs = 3, n_other=5,
  edge_prop = 0.75,
  mods = 'all'
)
dev.off()
# get the module assignment table:
modules <- GetModules(Drd2_data)


#DME analysis comparing two groups
group1 <- Drd2_data@meta.data %>% subset(phenotype == 'control') %>% rownames
group2 <- Drd2_data@meta.data %>% subset(phenotype == 'ROT') %>% rownames

head(group1)

DMEs <- FindDMEs(
  Drd2_data,
  barcodes1 = group1,
  barcodes2 = group2,
  test.use='wilcox',
  wgcna_name='Drd2_data'
)

head(DMEs)
PlotDMEsLollipop(
  Drd2_data, 
  DMEs, 
  wgcna_name='Drd2_data', 
  pvalue = "p_val_adj"
)
library(ggrepel)

#This plot shows the fold-change for each of the modules, 
#and the size of each dot corresponds to the number of genes in that module.
#An “X” is placed over each point that does not reach statistical significance.
PlotDMEsLollipop(
  Drd2_data, 
  DMEs, 
  wgcna_name='Drd2_data', 
  group.by = "metacell_grouping", 
  comparison = c("Drd2-MSN#control", "Drd2-MSN#ROT"),  
  pvalue = "p_val_adj"
) 
PlotDMEsVolcano(
  Drd2_data,
  DMEs,
  wgcna_name = 'Drd2_data'
)

Drd2_modules<-modules[modules$color!='grey',]
write_csv(Drd2_modules, 'Drd1_modules.csv')


# enrichr databases to test
dbs <- c('GO_Biological_Process_2021','GO_Cellular_Component_2021','GO_Molecular_Function_2021')

# perform enrichment tests
Drd2_data <- RunEnrichr(
  Drd2_data,
  dbs=dbs, # character vector of enrichr databases to test
  max_genes = 100 # number of genes per module to test
)
library(dplyr)
# retrieve the output table
Drd2_enrich_df <- GetEnrichrTable(Drd2_data)
write_csv(Drd2_enrich_df, 'Drd2_enrich_df.csv')
# enrichr dotplot

pdf("Drd2_wgcna_GObp.pdf", width=6, height=8)
EnrichrDotPlot(
  Drd2_data,
  mods = "all", # use all modules (this is the default behavior)
  database = "GO_Biological_Process_2021", # this has to be one of the lists we used above!!!
  n_terms=1 # number of terms for each module
)
dev.off()
#########motif#######
#BiocManager::install("BSgenome.Rnorvegicus.UCSC.rn5.masked")
# get the pfm from JASPAR2020 using TFBSTools
pfm_core <- TFBSTools::getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

Drd2_data <- dannymotif(
  Drd2_data,
  species_genome = 'BSgenome.Rnorvegicus.UCSC.rn7',
  pfm = pfm_core,
  EnsDb = EnsDb.Rnorvegicus.v105
)

# TF target genes
Drd2_TF_target_genes <- GetMotifTargets(Drd2_data)
# overlap between modules & TF target genes:
Drd2_data<- OverlapModulesMotifs(Drd2_data)
# look at the overlap data
head(GetMotifOverlap(Drd2_data))

######
df <- GetMotifOverlap(Drd2_data)
df2<-df
df2_no_duplicates <- df2[!duplicated(df2$tf), ]
table(df2_no_duplicates$fdr<0.05)
##Top 3 TF in DAE (HES1,HEY2,HES2 )
cur_df <- df %>% subset(tf == 'HES1')
plot_var <- 'odds_ratio'
p <- cur_df %>%
  ggplot(aes(y=reorder(module, odds_ratio), x=odds_ratio)) +
  geom_bar(stat='identity', fill=cur_df$color) +
  geom_vline(xintercept = 1, linetype='dashed', color='gray') +
  geom_text(aes(label=Significance), color='black', size=3.5, hjust='center') +
  ylab('') +
  xlab("Odds Ratio") +
  ggtitle("HES1 overlap") +
  theme(
    plot.title = element_text(hjust = 0.5)
  )


png('HES1_motif_overlap_or.png', width=3, height=4, units='in', res=400)
p
dev.off()

cur_df <- df %>% subset(tf == 'HEY2')
plot_var <- 'odds_ratio'
p <- cur_df %>%
  ggplot(aes(y=reorder(module, odds_ratio), x=odds_ratio)) +
  geom_bar(stat='identity', fill=cur_df$color) +
  geom_vline(xintercept = 1, linetype='dashed', color='gray') +
  geom_text(aes(label=Significance), color='black', size=3.5, hjust='center') +
  ylab('') +
  xlab("Odds Ratio") +
  ggtitle("HEY2 overlap") +
  theme(
    plot.title = element_text(hjust = 0.5)
  )


png('HEY2_motif_overlap_or.png', width=3, height=4, units='in', res=400)
p
dev.off()

cur_df <- df %>% subset(tf == 'HES2')
plot_var <- 'odds_ratio'
p <- cur_df %>%
  ggplot(aes(y=reorder(module, odds_ratio), x=odds_ratio)) +
  geom_bar(stat='identity', fill=cur_df$color) +
  geom_vline(xintercept = 1, linetype='dashed', color='gray') +
  geom_text(aes(label=Significance), color='black', size=3.5, hjust='center') +
  ylab('') +
  xlab("Odds Ratio") +
  ggtitle("HES2 overlap") +
  theme(
    plot.title = element_text(hjust = 0.5)
  )


png('HES2_motif_overlap_or.png', width=3, height=4, units='in', res=400)
p
dev.off()




#########pseudotime########
# convert the seurat object to CDS
cds <- as.cell_data_set(Drd2_data)
# run the monocle clustering
# assign paritions
reacreate.partition <- c(rep(1,length(cds@colData@rownames)))
names(reacreate.partition) <- cds@colData@rownames
reacreate.partition <- as.factor(reacreate.partition)
cds@clusters$UMAP$partitions <- reacreate.partition 
list_cluster <- Drd2_data@active.ident
cds@clusters$UMAP$clusters <- list_cluster
cds@int_colData@listData$reducedDims$UMAP <- Drd2_data@reductions$wnn.umap@cell.embeddings
cds <- cluster_cells(cds, reduction_method='UMAP')
# learn graph for pseudotime
cds <- learn_graph(cds)
plot_cells(
  cds = cds,
  color_cells_by = "phenotype",
  show_trajectory_graph = TRUE,
  label_principal_points = TRUE
) 













###########Drd1#############
DefaultAssay(Drd1_data)<-'RNA'
Drd1_data <- SetupForWGCNA(
  Drd1_data,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "Drd1_data" # the name of the hdWGCNA experiment
)

# construct metacells  in each group
Drd1_data <- MetacellsByGroups(
  seurat_obj = Drd1_data,
  group.by = c("predict.id", "phenotype"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'wnn.umap', # select the dimensionality reduction to perform KNN on
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'predict.id' # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
Drd1_data <- NormalizeMetacells(Drd1_data)
Drd1_data <- SetDatExpr(
  Drd1_data,
  group_name = "Drd1-MSN", # the name of the group of interest in the group.by column
  group.by='predict.id', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  slot = 'data' # using normalized data
)
# Test different soft powers:
Drd1_data <- TestSoftPowers(
  Drd1_data,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(Drd1_data)

# assemble with patchwork
wrap_plots(plot_list, ncol=2)
# construct co-expression network:
Drd1_data <- ConstructNetwork(
  Drd1_data, soft_power=5,
  setDatExpr=FALSE,overwrite_tom = TRUE,
  tom_name = 'Drd1-MSN' # name of the topoligical overlap matrix written to disk
)
PlotDendrogram(Drd1_data, main='Drd1-MSN hdWGCNA Dendrogram')

# need to run ScaleData first or else harmony throws an error:
Drd1_data <- ScaleData(Drd1_data, features=VariableFeatures(Drd1_data))

# compute all MEs in the full single-cell dataset
Drd1_data <- ModuleEigengenes(
  Drd1_data,
  group.by.vars="phenotype"
)
# harmonized module eigengenes:
hMEs <- GetMEs(Drd1_data)
# module eigengenes:
MEs <- GetMEs(Drd1_data, harmonized=FALSE)

# compute eigengene-based connectivity (kME):
Drd1_data <- ModuleConnectivity(
  Drd1_data,
  group.by = 'predict.id', group_name = 'Drd1-MSN'
)
# rename the modules
Drd1_data <- ResetModuleNames(
  Drd1_data,
  new_name = "Drd1-MSN"
)
# plot genes ranked by kME for each module

hub_df <-GetHubGenes(Drd1_data, n_hubs = 10)

pdf("Drd1_hubgebenw.pdf", width=6, height=8)
HubGeneNetworkPlot(
  Drd1_data,
  n_hubs = 3, n_other=5,
  edge_prop = 0.75,
  mods = 'all'
)
dev.off()
###########Astrocyte#############
DefaultAssay(Astrocyte_data)<-'RNA'
Astrocyte_data <- SetupForWGCNA(
  Astrocyte_data,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "Astrocyte_data" # the name of the hdWGCNA experiment
)

# construct metacells  in each group
Astrocyte_data <- MetacellsByGroups(
  seurat_obj = Astrocyte_data,
  group.by = c("predict.id", "phenotype"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'wnn.umap', # select the dimensionality reduction to perform KNN on
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'predict.id' # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
Astrocyte_data <- NormalizeMetacells(Astrocyte_data)
Astrocyte_data <- SetDatExpr(
  Astrocyte_data,
  group_name = "Astrocyte", # the name of the group of interest in the group.by column
  group.by='predict.id', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  slot = 'data' # using normalized data
)
# Test different soft powers:
Astrocyte_data <- TestSoftPowers(
  Astrocyte_data,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(Astrocyte_data)

# assemble with patchwork
pdf('soft_astrrocyte.pdf', width = 6, height = 8)
wrap_plots(plot_list, ncol=2)
dev.off()
# construct co-expression network:
Astrocyte_data <- ConstructNetwork(
  Astrocyte_data, soft_power=5,
  setDatExpr=FALSE,overwrite_tom = TRUE,
  tom_name = 'Astrocyte' # name of the topoligical overlap matrix written to disk
)
PlotDendrogram(Astrocyte_data, main='Astrocyte hdWGCNA Dendrogram')

# need to run ScaleData first or else harmony throws an error:
Astrocyte_data <- ScaleData(Astrocyte_data, features=VariableFeatures(Astrocyte_data))

# compute all MEs in the full single-cell dataset
Astrocyte_data <- ModuleEigengenes(
  Astrocyte_data,
  group.by.vars="phenotype"
)
# harmonized module eigengenes:
hMEs <- GetMEs(Astrocyte_data)
# module eigengenes:
MEs <- GetMEs(Astrocyte_data, harmonized=FALSE)

# compute eigengene-based connectivity (kME):
Astrocyte_data <- ModuleConnectivity(
  Astrocyte_data,
  group.by = 'predict.id', group_name = 'Astrocyte'
)
# rename the modules
Astrocyte_data <- ResetModuleNames(
  Astrocyte_data,
  new_name = "Astrocyte"
)
# plot genes ranked by kME for each module

hub_df <-GetHubGenes(Astrocyte_data, n_hubs = 10)

pdf("Astrocyte_hubgebenw.pdf", width=6, height=8)
HubGeneNetworkPlot(
  Astrocyte_data,
  n_hubs = 3, n_other=5,
  edge_prop = 0.75,
  mods = 'all'
)
dev.off()

DefaultAssay(Astrocyte_data)<-'RNA'
Astrocyte_data <- SetupForWGCNA(
  Astrocyte_data,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "Astrocyte_data" # the name of the hdWGCNA experiment
)

# construct metacells  in each group
Astrocyte_data <- MetacellsByGroups(
  seurat_obj = Astrocyte_data,
  group.by = c("predict.id", "phenotype"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'wnn.umap', # select the dimensionality reduction to perform KNN on
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'predict.id' # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
Astrocyte_data <- NormalizeMetacells(Astrocyte_data)
Astrocyte_data <- SetDatExpr(
  Astrocyte_data,
  group_name = "Astrocyte", # the name of the group of interest in the group.by column
  group.by='predict.id', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  slot = 'data' # using normalized data
)
# Test different soft powers:
Astrocyte_data <- TestSoftPowers(
  Astrocyte_data,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(Astrocyte_data)

# assemble with patchwork
pdf('soft_astrrocyte.pdf', width = 6, height = 8)
wrap_plots(plot_list, ncol=2)
dev.off()
# construct co-expression network:
Astrocyte_data <- ConstructNetwork(
  Astrocyte_data, soft_power=5,
  setDatExpr=FALSE,overwrite_tom = TRUE,
  tom_name = 'Astrocyte' # name of the topoligical overlap matrix written to disk
)
PlotDendrogram(Astrocyte_data, main='Astrocyte hdWGCNA Dendrogram')

# need to run ScaleData first or else harmony throws an error:
Astrocyte_data <- ScaleData(Astrocyte_data, features=VariableFeatures(Astrocyte_data))

# compute all MEs in the full single-cell dataset
Astrocyte_data <- ModuleEigengenes(
  Astrocyte_data,
  group.by.vars="phenotype"
)
# harmonized module eigengenes:
hMEs <- GetMEs(Astrocyte_data)
# module eigengenes:
MEs <- GetMEs(Astrocyte_data, harmonized=FALSE)

# compute eigengene-based connectivity (kME):
Astrocyte_data <- ModuleConnectivity(
  Astrocyte_data,
  group.by = 'predict.id', group_name = 'Astrocyte'
)
# rename the modules
Astrocyte_data <- ResetModuleNames(
  Astrocyte_data,
  new_name = "Astrocyte"
)
# plot genes ranked by kME for each module

hub_df <-GetHubGenes(Astrocyte_data, n_hubs = 10)

pdf("Astrocyte_hubgebenw.pdf", width=6, height=8)
HubGeneNetworkPlot(
  Astrocyte_data,
  n_hubs = 3, n_other=5,
  edge_prop = 0.75,
  mods = 'all'
)
dev.off()






###########Olig#############
DefaultAssay(Olig_data)<-'RNA'
Olig_data <- SetupForWGCNA(
  Olig_data,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "Olig_data" # the name of the hdWGCNA experiment
)

# construct metacells  in each group
Olig_data <- MetacellsByGroups(
  seurat_obj = Olig_data,
  group.by = c("predict.id", "phenotype"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'wnn.umap', # select the dimensionality reduction to perform KNN on
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'predict.id' # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
Olig_data <- NormalizeMetacells(Olig_data)
Olig_data <- SetDatExpr(
  Olig_data,
  group_name = "Olig", # the name of the group of interest in the group.by column
  group.by='predict.id', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  slot = 'data' # using normalized data
)
# Test different soft powers:
Olig_data <- TestSoftPowers(
  Olig_data,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(Olig_data)

# assemble with patchwork
pdf('soft_olig.pdf', width = 6, height = 8)
wrap_plots(plot_list, ncol=2)
dev.off()
# construct co-expression network:
Olig_data <- ConstructNetwork(
  Olig_data, soft_power=5,
  setDatExpr=FALSE,overwrite_tom = TRUE,
  tom_name = 'olig' # name of the topoligical overlap matrix written to disk
)
PlotDendrogram(Olig_data, main='Olig hdWGCNA Dendrogram')

# need to run ScaleData first or else harmony throws an error:
Olig_data <- ScaleData(Olig_data, features=VariableFeatures(Olig_data))

# compute all MEs in the full single-cell dataset
Olig_data <- ModuleEigengenes(
  Olig_data,
  group.by.vars="phenotype"
)
# harmonized module eigengenes:
hMEs <- GetMEs(Olig_data)
# module eigengenes:
MEs <- GetMEs(Olig_data, harmonized=FALSE)

# compute eigengene-based connectivity (kME):
Olig_data <- ModuleConnectivity(
  Olig_data,
  group.by = 'predict.id', group_name = 'Olig'
)
# rename the modules
Olig_data <- ResetModuleNames(
  Olig_data,
  new_name = "Olig"
)
# plot genes ranked by kME for each module

hub_df <-GetHubGenes(Olig_data)

pdf("Olig_hubgebenw.pdf", width=6, height=8)
HubGeneNetworkPlot(
  Olig_data,
  n_hubs = 3, n_other=5,
  edge_prop = 0.75,
  mods = 'all'
)
dev.off()




###########Micro#############
DefaultAssay(Microglia_data)<-'RNA'
Microglia_data <- SetupForWGCNA(
  Microglia_data,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "Microglia" # the name of the hdWGCNA experiment
)

# construct metacells  in each group
Microglia_data <- MetacellsByGroups(
  seurat_obj = Microglia_data,
  group.by = c("predict.id", "phenotype"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'wnn.umap', # select the dimensionality reduction to perform KNN on
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'predict.id' # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
Microglia_data <- NormalizeMetacells(Microglia_data)
Microglia_data <- SetDatExpr(
  Microglia_data,
  group_name = "Microglia", # the name of the group of interest in the group.by column
  group.by='predict.id', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  slot = 'data' # using normalized data
)
# Test different soft powers:
Microglia_data <- TestSoftPowers(
  Microglia_data,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(Microglia_data)

# assemble with patchwork
pdf('soft_micro.pdf', width = 6, height = 8)
wrap_plots(plot_list, ncol=2)
dev.off()
# construct co-expression network:
Microglia_data <- ConstructNetwork(
  Microglia_data, soft_power=4,
  setDatExpr=FALSE,overwrite_tom = TRUE,
  tom_name = 'Microglia' # name of the topoligical overlap matrix written to disk
)
PlotDendrogram(Microglia_data, main='Micro_hdWGCNA Dendrogram')

# need to run ScaleData first or else harmony throws an error:
Microglia_data <- ScaleData(Microglia_data, features=VariableFeatures(Microglia_data))

# compute all MEs in the full single-cell dataset
Microglia_data <- ModuleEigengenes(
  Microglia_data,
  group.by.vars="phenotype"
)
# harmonized module eigengenes:
hMEs <- GetMEs(Microglia_data)
# module eigengenes:
MEs <- GetMEs(Microglia_data, harmonized=FALSE)

# compute eigengene-based connectivity (kME):
Microglia_data <- ModuleConnectivity(
  Microglia_data,
  group.by = 'predict.id', group_name = 'Microglia'
)
# rename the modules
Microglia_data <- ResetModuleNames(
  Microglia_data,
  new_name = "Microglia"
)
# plot genes ranked by kME for each module

hub_df <-GetHubGenes(Microglia_data)

pdf("Micro_hubgebenw.pdf", width=6, height=8)
HubGeneNetworkPlot(
  Microglia_data,
  n_hubs = 3, n_other=5,
  edge_prop = 0.75,
  mods = 'all'
)
dev.off()



