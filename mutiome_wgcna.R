rm(list = ls())
# single-cell analysis package
library(Seurat)

# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)
library(devtools)
# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)
library(GeneOverlap)
# using the cowplot theme for ggplot
theme_set(theme_cowplot())
set.seed(12345)
data<-readRDS("mutiome_final.rds")
DefaultAssay(data)<-'RNA'
DimPlot(data, group.by='predict.id', label=TRUE)

data <- SetupForWGCNA(
  data,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "tutorial" # the name of the hdWGCNA experiment
)
# construct metacells  in each group
data <- MetacellsByGroups(
  seurat_obj = data,
  group.by = c("predict.id", "phenotype"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'wnn.umap', # select the dimensionality reduction to perform KNN on
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'predict.id' # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
data <- NormalizeMetacells(data)
##############
# get list of cell populations
groups <- unique(GetMetacellObject(data)$predict.id)

# set up gene expression matrix
data <- SetDatExpr(
  data,
  group.by='predict.id',
  group_name = groups,
  use_metacells=TRUE,
  slot = 'data',
  assay = 'RNA'
)
# Test different soft powers:
data <- TestSoftPowers(
  data,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# test soft powers
data <- TestSoftPowers(data)
plot_list <- PlotSoftPowers(data)

# assemble with patchwork
pdf("soft_power_hdWGCNA_Dendrogram.pdf", width=8,height=4)
wrap_plots(plot_list, ncol=2)
dev.off()
##################
# construct metacells  in each group
data <- MetacellsByGroups(
  seurat_obj = data,
  group.by = c("predict.id", "phenotype"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'wnn.umap', # select the dimensionality reduction to perform KNN on
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'predict.id' # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
data <- NormalizeMetacells(data)

data <- SetDatExpr(
  data,
  group_name = c("Drd2-MSN", "Drd1-MSN","Olig","Astrocyte","Microglia"
                 ,"OPC","Glutamatergic","Interneuron"),
  group.by='predict.id',assay = "RNA"
)
# Test different soft powers:
data <- TestSoftPowers(
  data,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(data)

# assemble with patchwork
wrap_plots(plot_list, ncol=2)
############Drd2#################
data <- SetupForWGCNA(
  data,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "drd2" # the name of the hdWGCNA experiment
)

# construct metacells  in each group
data <- MetacellsByGroups(
  seurat_obj = data,
  group.by = c("predict.id"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'wnn.umap', # select the dimensionality reduction to perform KNN on
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'predict.id' # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
data <- NormalizeMetacells(data)



data <- SetDatExpr(
  data,
  group_name = "Drd2-MSN", # the name of the group of interest in the group.by column
  group.by='predict.id', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  slot = 'data' # using normalized data
)
# Test different soft powers:
data <- TestSoftPowers(
  data,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(data)

# assemble with patchwork
wrap_plots(plot_list, ncol=2)

# construct co-expression network:
data <- ConstructNetwork(
  data, soft_power=5,
  setDatExpr=FALSE,
  tom_name = 'Drd2', overwrite_tom = TRUE # name of the topoligical overlap matrix written to disk
)
pdf('Drd2_hdWGCNA_Dendrogram.pdf', width=8, height=4)
PlotDendrogram(data, main='Drd2 hdWGCNA Dendrogram')
dev.off()

# need to run ScaleData first or else harmony throws an error:
data <- ScaleData(data, features=VariableFeatures(data))

# compute all MEs in the full single-cell dataset
data <- ModuleEigengenes(
  data,
  group.by.vars="phenotype", exclude_grey = TRUE
)

# harmonized module eigengenes:
hMEs <- GetMEs(data)
# compute eigengene-based connectivity (kME):
data <- ModuleConnectivity(
  data,
  group.by = 'predict.id', harmonized = FALSE
)

# rename the modules
data <- ResetModuleNames(
  data,
  new_name = "Drd2"
)
# get the module assignment table:
modules <- GetModules(data)
write_csv(modules, 'Drd2_modules.csv')
# show the first 6 columns:
head(modules[,1:6])
# get hub genes
hub_df <- GetHubGenes(data, n_hubs = 10)

head(hub_df)

# compute gene scoring for the top 25 hub genes by kME for each module
# with UCell method
library(UCell)
data <- ModuleExprScore(
  data,
  n_genes = 25,
  method='UCell'
)
# plot module correlagram
library(corrplot)
ModuleCorrelogram(data)

# get hMEs from seurat object
MEs <- GetMEs(data, harmonized=TRUE)
mods <- colnames(MEs); mods <- mods[mods != 'grey']

# add hMEs to Seurat meta-data:
data@meta.data <- cbind(data@meta.data, MEs)

##
# gene enrichment packages
library(enrichR)
library(GeneOverlap)
# enrichr databases to test
dbs <- c('GO_Biological_Process_2021','GO_Cellular_Component_2021','GO_Molecular_Function_2021')

# perform enrichment tests
data <- RunEnrichr(
  data,
  dbs=dbs, # character vector of enrichr databases to test
  max_genes = 100 # number of genes per module to test
)



# enrichr dotplot
pdf('Drd2_EnrichrDotPlot.pdf', width=8, height=9)
EnrichrDotPlot(
  data,
  mods = "all", # use all modules (this is the default behavior)
  database = "GO_Biological_Process_2021", # this has to be one of the lists we used above!!!
  n_terms=1 # number of terms for each module
)
dev.off()
###### test
# enrichr databases to test
dbs <- c('GO_Biological_Process_2021','GO_Cellular_Component_2021','GO_Molecular_Function_2021')

# perform enrichment tests
data <- RunEnrichr(
  data,
  dbs=dbs, # character vector of enrichr databases to test
  max_genes = 100 # number of genes per module to test
)

# retrieve the output table
enrich_df <- GetEnrichrTable(seurat_obj)






##########Drd1#################
data <- SetDatExpr(
  data,
  group_name = "Drd1-MSN", # the name of the group of interest in the group.by column
  group.by='predict.id', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  slot = 'data' # using normalized data
)
# Test different soft powers:
data <- TestSoftPowers(
  data,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(data)

# assemble with patchwork
wrap_plots(plot_list, ncol=2)

# construct co-expression network:
data <- ConstructNetwork(
  data, soft_power=5,
  setDatExpr=FALSE,
  tom_name = 'Drd1',overwrite_tom = TRUE # name of the topoligical overlap matrix written to disk
)
pdf('Drd1_hdWGCNA_Dendrogram.pdf', width=8, height=4)
PlotDendrogram(data, main='Drd1 hdWGCNA Dendrogram')
dev.off()
# need to run ScaleData first or else harmony throws an error:
data <- ScaleData(data, features=VariableFeatures(data))

# compute all MEs in the full single-cell dataset
data <- ModuleEigengenes(
  data,
  group.by.vars="phenotype"
)

# harmonized module eigengenes:
hMEs <- GetMEs(data)
# compute eigengene-based connectivity (kME):
data <- ModuleConnectivity(
  data,
  group.by = 'predict.id', group_name = 'Drd1'
)

# rename the modules
#data <- ResetModuleNames(
#  data,
#  new_name = "Microglia"
#)
# get the module assignment table:
modules <- GetModules(data)
write_csv(modules, 'Drd1_modules.csv')
# show the first 6 columns:
head(modules[,1:6])
# get hub genes
hub_df <- GetHubGenes(data, n_hubs = 10)

head(hub_df)

# compute gene scoring for the top 25 hub genes by kME for each module
# with UCell method
library(UCell)
data <- ModuleExprScore(
  data,
  n_genes = 25,
  method='UCell'
)
# plot module correlagram
library(corrplot)
ModuleCorrelogram(data)

# get hMEs from seurat object
MEs <- GetMEs(data, harmonized=TRUE)
mods <- colnames(MEs); mods <- mods[mods != 'grey']

# add hMEs to Seurat meta-data:
data@meta.data <- cbind(data@meta.data, MEs)

##
# gene enrichment packages
library(enrichR)
library(GeneOverlap)
# enrichr databases to test
dbs <- c('GO_Biological_Process_2021','GO_Cellular_Component_2021','GO_Molecular_Function_2021')

# perform enrichment tests
data <- RunEnrichr(
  data,
  dbs=dbs, # character vector of enrichr databases to test
  max_genes = 100 # number of genes per module to test
)



# enrichr dotplot
pdf('Drd1_EnrichrDotPlot.pdf', width=8, height=9)
EnrichrDotPlot(
  data,
  mods = "all", # use all modules (this is the default behavior)
  database = "GO_Biological_Process_2021", # this has to be one of the lists we used above!!!
  n_terms=1 # number of terms for each module
)
dev.off()
#########trajectory############
library(monocle3)
library(SeuratWrappers)
data<-readRDS("mutiome_final.rds")
# convert the seurat object to CDS
#choose celltype 
#Drd2-MSN Drd1-MSN Olig Astrocyte Microglia OPC Glutamatergic Interneuron
data<-subset(x=data, subset = predict.id =='Microglia')



data$combo_group <- paste0(data$predict.id, "_", data$phenotype)
cds <- as.cell_data_set(data)
# assign paritions
reacreate.partition <- c(rep(1,length(cds@colData@rownames)))
names(reacreate.partition) <- cds@colData@rownames
reacreate.partition <- as.factor(reacreate.partition)
cds@clusters$UMAP$partitions <- reacreate.partition                    
# Assign the cluster info 

list_cluster <- data@active.ident
cds@clusters$UMAP$clusters <- list_cluster
cds@int_colData@listData$reducedDims$UMAP <- data@reductions$wnn.umap@cell.embeddings

cluster.before.trajectory <- plot_cells(cds,
                                        color_cells_by = 'cluster',
                                        label_groups_by_cluster = FALSE,
                                        group_label_size = 5) +
  theme(legend.position = "right")

cluster.names <- plot_cells(cds,
                            color_cells_by = "predict.id",
                            label_groups_by_cluster = FALSE,
                            group_label_size = 5) +
  scale_color_manual(values = c('blue','red', 'green', 'maroon', 'yellow', 'grey', 'cyan')) +
  theme(legend.position = "right")

cluster.before.trajectory | cluster.names

cds <- learn_graph(cds)
# plot the pseudotime graph:
p1 <- plot_cells(
  cds = cds,
  color_cells_by = "predict.id",
  show_trajectory_graph = TRUE,
  label_principal_points = TRUE
) 
# plot the UMAP partitions from the clustering algorithm
p2 <-  plot_cells(
  cds = cds,
  color_cells_by = "phenotype",
  show_trajectory_graph = TRUE,
  label_principal_points = TRUE
) 

p3 <-  plot_cells(
  cds = cds,
  color_cells_by = "combo_group",
  show_trajectory_graph = TRUE,
  label_principal_points = TRUE
) 

png(('Astrocyteumap_monocle.png'),  width=8, height=4, res=500, units='in')
p1 + p3
dev.off()


# get principal node & order cells
principal_node <- 'Y_13'
cds <- order_cells(cds,root_pr_nodes = principal_node)

# add pseudotime to seurat object:
data$pseudotime <- pseudotime(cds)

# separate pseudotime trajectories by the different mature cells
data$drd1_pseudotime <- ifelse(data$combo_group %in% c("Drd1-MSN_control", "Drd1-MSN_ROT"), data$pseudotime, NA)
data$drd2_pseudotime <- ifelse(data$combo_group %in% c("Drd2-MSN_control", "Drd2-MSN_ROT"), data$pseudotime, NA)
data$glu_pseudotime <- ifelse(data$combo_group %in% c("Drd1-MSN_control", "Drd1-MSN_ROT", "Drd2-MSN_control", "Drd2-MSN_ROT",'Glutamatergic_control','Glutamatergic_ROT' ), data$pseudotime, NA)
data$olig_pseudotime <- ifelse(data$combo_group %in% c("Astrocyte_control", "Astrocyte_ROT"), data$pseudotime, NA)
#seurat_obj$mono_pseudotime <- ifelse(seurat_obj$celltype %in% c("HSC", "HMP", 'Mono'), seurat_obj$pseudotime, NA)
#seurat_obj$dc_pseudotime <- ifelse(seurat_obj$celltype %in% c("HSC", "HMP", 'DCPre', 'cDC', 'pDC'), seurat_obj$pseudotime, NA)
#seurat_obj$clp_pseudotime <- ifelse(seurat_obj$celltype %in% c("HSC", "HMP", 'CLP'), seurat_obj$pseudotime, NA)

data$UMAP1 <- data@reductions$wnn.umap@cell.embeddings[,1]
data$UMAP2 <- data@reductions$wnn.umap@cell.embeddings[,2]
library(viridis)

p1 <- data@meta.data %>%
  ggplot(aes(x=UMAP1, y=UMAP2, color=drd1_pseudotime)) +
  ggrastr::rasterise(geom_point(size=1), dpi=500, scale=0.75) +
  coord_equal() +
  scale_color_gradientn(colors=plasma(256), na.value='grey') +
  umap_theme()

p2 <- data@meta.data %>%
  ggplot(aes(x=UMAP1, y=UMAP2, color=drd2_pseudotime)) +
  ggrastr::rasterise(geom_point(size=1), dpi=500, scale=0.75) +
  coord_equal() +
  scale_color_gradientn(colors=viridis(256), na.value='grey') +
  umap_theme()
p3 <- data@meta.data %>%
  ggplot(aes(x=UMAP1, y=UMAP2, color=glu_pseudotime)) +
  ggrastr::rasterise(geom_point(size=1), dpi=500, scale=0.75) +
  coord_equal() +
  scale_color_gradientn(colors=viridis(256), na.value='grey') +
  umap_theme()
p4 <- data@meta.data %>%
  ggplot(aes(x=UMAP1, y=UMAP2, color=olig_pseudotime)) +
  ggrastr::rasterise(geom_point(size=1), dpi=500, scale=0.75) +
  coord_equal() +
  scale_color_gradientn(colors=viridis(256), na.value='grey')+
  umap_theme()
# assemble with patchwork
(p1 | p2) / (p3 + p4) + plot_layout(ncol=1, guides='collect')

#############################
# construct co-expression network:
data <- ConstructNetwork(
  data,
  soft_power=5,
  setDatExpr=FALSE,
  detectCutHeight=0.995,
  mergeCutHeight=0.05,
  TOM_name = 'variable',
  minModuleSize=30,
  overwrite_tom = TRUE
)
pdf("hdWGCNA_Dendrogram.pdf", width=8,height=4)
PlotDendrogram(data, main='hdWGCNA Dendrogram')
dev.off()

data <- ScaleData(data, features=rownames(data))
# harmony correction by Sample
data <- ModuleEigengenes(
  data,
  group.by.vars = 'phenotype',
  verbose=TRUE
)       


#install.packages("qlcMatrix")
# compute module connectivity:
modules_orig <- GetModules(data)
data <- ModuleConnectivity(data)
modules <- GetModules(data)


# get the MEisos from the seurat object and add it to the metadata
MEiso <- GetMEs(data)
meta <- data@meta.data
data@meta.data <- cbind(meta, MEiso)

# get a list of features to plot
modules <- GetModules(data)
mods <- levels(modules$module)
mods <- mods[mods!='grey']

# make dotplot
p <- DotPlot(
  data,
  group.by='predict.id',
  features = rev(mods)
) + RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue') + xlab('') + ylab('') +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  )

# show plot

pdf("dotplot_hdWGCNA.pdf", width=8,height=4)
p
dev.off()
# restore the original metadata
data@meta.data <- meta

#########
data <- RunModuleUMAP(
  data,
  n_hubs =5,
  n_neighbors=15,
  min_dist=0.2,
  spread=1
  #supervised=TRUE,
  #target_weight=0.5
)


# get the hub gene UMAP table from the seurat object
umap_df <- GetModuleUMAP(
  data
)

# plot with ggplot
p <- ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
  geom_point(
    color=umap_df$color,
    size=umap_df$kME*2
  ) +
  umap_theme()

pdf('variable_test_hubgene_umap_ggplot_uns.pdf', width=8, height=4)
p
dev.off()

png('variable_coex_umap.png', width=6, height=6, units='in', res=500)
ModuleUMAPPlot(
  data,
  edge.alpha=0.5,
  sample_edges=TRUE,
  keep_grey_edges=FALSE,
  edge_prop=0.075, # taking the top 20% strongest edges in each module
  #   label_genes = umap_df$isoform_name[umap_df$gene_name %in% plot_genes],
  #label_genes = c(''),
  label_hubs=2 # how many hub genes to plot per module?
)
dev.off()

saveRDS(data, file=paste0(data_dir, 'hdWGCNA_wip.rds'))

#######
# hubgene network
pdf('HubGeneNetworkPlot.pdf', width=8, height=4)
HubGeneNetworkPlot(
  data,
  n_hubs = 3, n_other=5,
  edge_prop = 0.75,
  mods = 'all'
)
dev.off()
##########Olig#################
# construct co-expression network:
data <- ConstructNetwork(
  data, soft_power=5,
  setDatExpr=FALSE,
  tom_name = 'Olig',overwrite_tom = TRUE # name of the topoligical overlap matrix written to disk
)

PlotDendrogram(data, main='Olig hdWGCNA Dendrogram')

# need to run ScaleData first or else harmony throws an error:
data <- ScaleData(data, features=VariableFeatures(data))

# compute all MEs in the full single-cell dataset
data <- ModuleEigengenes(
  data,
  group.by.vars="phenotype"
)

# harmonized module eigengenes:
hMEs <- GetMEs(data)
# compute eigengene-based connectivity (kME):
data <- ModuleConnectivity(
  data,
  group.by = 'predict.id', group_name = 'Olig'
)

# rename the modules
data <- ResetModuleNames(
  data,
  new_name = "Olig"
)
# get the module assignment table:
modules <- GetModules(data)

# show the first 6 columns:
head(modules[,1:6])
# get hub genes
hub_df <- GetHubGenes(data, n_hubs = 10)

head(hub_df)

# compute gene scoring for the top 25 hub genes by kME for each module
# with UCell method
library(UCell)
data <- ModuleExprScore(
  data,
  n_genes = 25,
  method='UCell'
)
# plot module correlagram
library(corrplot)
ModuleCorrelogram(data)

# get hMEs from seurat object
MEs <- GetMEs(data, harmonized=TRUE)
mods <- colnames(MEs); mods <- mods[mods != 'grey']

# add hMEs to Seurat meta-data:
data@meta.data <- cbind(data@meta.data, MEs)

##
# gene enrichment packages
library(enrichR)
library(GeneOverlap)
# enrichr databases to test
dbs <- c('GO_Biological_Process_2021','GO_Cellular_Component_2021','GO_Molecular_Function_2021')

# perform enrichment tests
data <- RunEnrichr(
  data,
  dbs=dbs, # character vector of enrichr databases to test
  max_genes = 100 # number of genes per module to test
)

# make GO term plots:
EnrichrBarPlot(
  data,
  outdir = "enrichr_plots_Olig", # name of output directory
  n_terms = 10, # number of enriched terms to show (sometimes more show if there are ties!!!)
  plot_size = c(5,7), # width, height of the output .pdfs
  logscale=TRUE # do you want to show the enrichment as a log scale?
)


# retrieve the output table
enrich_df <- GetEnrichrTable(data)
# enrichr dotplot
EnrichrDotPlot(
  data,
  mods = "all", # use all modules (this is the default behavior)
  database = "GO_Biological_Process_2021", # this has to be one of the lists we used above!!!
  n_terms=1 # number of terms for each module
)




##########Astrocyte#######
data <- SetDatExpr(
  data,
  group_name = "Astrocyte", # the name of the group of interest in the group.by column
  group.by='predict.id', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  slot = 'data' # using normalized data
)
# Test different soft powers:
data <- TestSoftPowers(
  data,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(data)

# assemble with patchwork
wrap_plots(plot_list, ncol=2)

# construct co-expression network:
data <- ConstructNetwork(
  data, soft_power=4,
  setDatExpr=FALSE,
  tom_name = 'Astrocyte',overwrite_tom = TRUE # name of the topoligical overlap matrix written to disk
)
pdf('Astrocyte_hdWGCNA_Dendrogram.pdf', width=8, height=4)
PlotDendrogram(data, main='Astrocyte hdWGCNA Dendrogram')
dev.off()
# need to run ScaleData first or else harmony throws an error:
data <- ScaleData(data, features=VariableFeatures(data))

# compute all MEs in the full single-cell dataset
data <- ModuleEigengenes(
  data,
  group.by.vars="phenotype"
)

# harmonized module eigengenes:
hMEs <- GetMEs(data)
# compute eigengene-based connectivity (kME):
data <- ModuleConnectivity(
  data,
  group.by = 'predict.id', group_name = 'Astrocyte'
)

# rename the modules
data <- ResetModuleNames(
  data,
  new_name = "Astrocyte"
)
# get the module assignment table:
modules <- GetModules(data)
write_csv(modules, 'Astrocyte_modules.csv')
# show the first 6 columns:
head(modules[,1:6])
# get hub genes
hub_df <- GetHubGenes(data, n_hubs = 10)

head(hub_df)

# compute gene scoring for the top 25 hub genes by kME for each module
# with UCell method
library(UCell)
data <- ModuleExprScore(
  data,
  n_genes = 25,
  method='UCell'
)
# plot module correlagram
library(corrplot)
ModuleCorrelogram(data)

# get hMEs from seurat object
MEs <- GetMEs(data, harmonized=TRUE)
mods <- colnames(MEs); mods <- mods[mods != 'grey']

# add hMEs to Seurat meta-data:
data@meta.data <- cbind(data@meta.data, MEs)

##
# gene enrichment packages
library(enrichR)
library(GeneOverlap)
# enrichr databases to test
dbs <- c('GO_Biological_Process_2021','GO_Cellular_Component_2021','GO_Molecular_Function_2021')

# perform enrichment tests
data <- RunEnrichr(
  data,
  dbs=dbs, # character vector of enrichr databases to test
  max_genes = 100 # number of genes per module to test
)



# enrichr dotplot
pdf('Astrocyte_EnrichrDotPlot.pdf', width=8, height=4)
EnrichrDotPlot(
  data,
  mods = "all", # use all modules (this is the default behavior)
  database = "GO_Biological_Process_2021", # this has to be one of the lists we used above!!!
  n_terms=1 # number of terms for each module
)
dev.off()
##########Micriglia########
data <- SetDatExpr(
  data,
  group_name = "Microglia", # the name of the group of interest in the group.by column
  group.by='predict.id', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  slot = 'data' # using normalized data
)
# Test different soft powers:
data <- TestSoftPowers(
  data,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(data)

# assemble with patchwork
wrap_plots(plot_list, ncol=2)

# construct co-expression network:
data <- ConstructNetwork(
  data, soft_power=4,
  setDatExpr=FALSE,
  tom_name = 'Microglia',overwrite_tom = TRUE # name of the topoligical overlap matrix written to disk
)
pdf('Microglia_hdWGCNA_Dendrogram.pdf', width=8, height=4)
PlotDendrogram(data, main='Microglia hdWGCNA Dendrogram')
dev.off()
# need to run ScaleData first or else harmony throws an error:
data <- ScaleData(data, features=VariableFeatures(data))

# compute all MEs in the full single-cell dataset
data <- ModuleEigengenes(
  data,
  group.by.vars="phenotype"
)

# harmonized module eigengenes:
hMEs <- GetMEs(data)
# compute eigengene-based connectivity (kME):
data <- ModuleConnectivity(
  data,
  group.by = 'predict.id', group_name = 'Microglia'
)

# rename the modules
data <- ResetModuleNames(
  data,
  new_name = "Microglia"
)
# get the module assignment table:
modules <- GetModules(data)
write_csv(modules, 'Microglia_modules.csv')
# show the first 6 columns:
head(modules[,1:6])
# get hub genes
hub_df <- GetHubGenes(data, n_hubs = 10)

head(hub_df)

# compute gene scoring for the top 25 hub genes by kME for each module
# with UCell method
library(UCell)
data <- ModuleExprScore(
  data,
  n_genes = 25,
  method='UCell'
)
# plot module correlagram
library(corrplot)
ModuleCorrelogram(data)

# get hMEs from seurat object
MEs <- GetMEs(data, harmonized=TRUE)
mods <- colnames(MEs); mods <- mods[mods != 'grey']

# add hMEs to Seurat meta-data:
data@meta.data <- cbind(data@meta.data, MEs)

##
# gene enrichment packages
library(enrichR)
library(GeneOverlap)
# enrichr databases to test
dbs <- c('GO_Biological_Process_2021','GO_Cellular_Component_2021','GO_Molecular_Function_2021')

# perform enrichment tests
data <- RunEnrichr(
  data,
  dbs=dbs, # character vector of enrichr databases to test
  max_genes = 100 # number of genes per module to test
)



# enrichr dotplot
pdf('Microglia_EnrichrDotPlot.pdf', width=8, height=4)
EnrichrDotPlot(
  data,
  mods = "all", # use all modules (this is the default behavior)
  database = "GO_Biological_Process_2021", # this has to be one of the lists we used above!!!
  n_terms=1 # number of terms for each module
)
dev.off()
##########OPC#######
data <- SetDatExpr(
  data,
  group_name = "OPC", # the name of the group of interest in the group.by column
  group.by='predict.id', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  slot = 'data' # using normalized data
)
# Test different soft powers:
data <- TestSoftPowers(
  data,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(data)

# assemble with patchwork
wrap_plots(plot_list, ncol=2)

# construct co-expression network:
data <- ConstructNetwork(
  data, soft_power=4,
  setDatExpr=FALSE,
  tom_name = 'OPC',overwrite_tom = TRUE # name of the topoligical overlap matrix written to disk
)
pdf('OPC_hdWGCNA_Dendrogram.pdf', width=8, height=4)
PlotDendrogram(data, main='OPC hdWGCNA Dendrogram')
dev.off()
# need to run ScaleData first or else harmony throws an error:
data <- ScaleData(data, features=VariableFeatures(data))

# compute all MEs in the full single-cell dataset
data <- ModuleEigengenes(
  data,
  group.by.vars="phenotype"
)

# harmonized module eigengenes:
hMEs <- GetMEs(data)
# compute eigengene-based connectivity (kME):
data <- ModuleConnectivity(
  data,
  group.by = 'predict.id', group_name = 'OPC'
)

# rename the modules
#data <- ResetModuleNames(
#  data,
#  new_name = "Microglia"
#)
# get the module assignment table:
modules <- GetModules(data)
write_csv(modules, 'OPC_modules.csv')
# show the first 6 columns:
head(modules[,1:6])
# get hub genes
hub_df <- GetHubGenes(data, n_hubs = 10)

head(hub_df)

# compute gene scoring for the top 25 hub genes by kME for each module
# with UCell method
library(UCell)
data <- ModuleExprScore(
  data,
  n_genes = 25,
  method='UCell'
)
# plot module correlagram
library(corrplot)
ModuleCorrelogram(data)

# get hMEs from seurat object
MEs <- GetMEs(data, harmonized=TRUE)
mods <- colnames(MEs); mods <- mods[mods != 'grey']

# add hMEs to Seurat meta-data:
data@meta.data <- cbind(data@meta.data, MEs)

##
# gene enrichment packages
library(enrichR)
library(GeneOverlap)
# enrichr databases to test
dbs <- c('GO_Biological_Process_2021','GO_Cellular_Component_2021','GO_Molecular_Function_2021')

# perform enrichment tests
data <- RunEnrichr(
  data,
  dbs=dbs, # character vector of enrichr databases to test
  max_genes = 100 # number of genes per module to test
)



# enrichr dotplot
pdf('OPC_EnrichrDotPlot.pdf', width=8, height=9)
EnrichrDotPlot(
  data,
  mods = "all", # use all modules (this is the default behavior)
  database = "GO_Biological_Process_2021", # this has to be one of the lists we used above!!!
  n_terms=1 # number of terms for each module
)
dev.off()
###########
# make GO term plots:
# enrichr databases to test
dbs <- c('GO_Biological_Process_2021','GO_Cellular_Component_2021','GO_Molecular_Function_2021')

# perform enrichment tests
data <- RunEnrichr(
  data,
  dbs=dbs, # character vector of enrichr databases to test
  max_genes = 100 # number of genes per module to test
)

# retrieve the output table
enrich_df <- GetEnrichrTable(data)
EnrichrBarPlot(
  data,
  outdir = "enrichr_plots", # name of output directory
  n_terms = 10, # number of enriched terms to show (sometimes more show if there are ties!!!)
  plot_size = c(5,7), # width, height of the output .pdfs
  logscale=TRUE # do you want to show the enrichment as a log scale?
)
pdf('EnrichrBarPlot.pdf', width=8, height=4)
# enrichr dotplot
EnrichrDotPlot(
  data,
  mods = "all", # use all modules (this is the default behavior)
  database = "GO_Biological_Process_2021", # this has to be one of the lists we used above!!!
  n_terms=1 # number of terms for each module
)
dev.off()

##Marker gene overlap analysis
# compute cell-type marker genes with Seurat:
Idents(data) <- data$predict.id
markers <- Seurat::FindAllMarkers(
  data,
  only.pos = TRUE,
  logfc.threshold=1
)

# compute marker gene overlaps
overlap_df <- OverlapModulesDEGs(
  data,
  deg_df = markers,
  fc_cutoff = 1 # log fold change cutoff for overlap analysis
)

#Visualize overlaps
# overlap barplot, produces a plot for each cell type
plot_list <- OverlapBarPlot(overlap_df)

# stitch plots with patchwork
pdf('OverlapBarPlot.pdf', width=8, height=4)
wrap_plots(plot_list, ncol=3)
dev.off()

library(monocle3)
library(SeuratWrappers)

# convert the seurat object to CDS
cds <- as.cell_data_set(seurat_obj)

# run the monocle clustering
cds <- cluster_cells(cds, reduction_method='UMAP')

# learn graph for pseudotime
cds <- learn_graph(cds)

# plot the pseudotime graph:
p1 <- plot_cells(
  cds = cds,
  color_cells_by = "celltype",
  show_trajectory_graph = TRUE,
  label_principal_points = TRUE
) 

# plot the UMAP partitions from the clustering algorithm
p2 <- plot_cells(
  cds = cds,
  color_cells_by = "partition",
  show_trajectory_graph = FALSE
)

png(paste0(fig_dir, 'umap_monocle.png'),  width=8, height=4, res=500, units='in')
p1 + p2
dev.off()
##########
#Monocle3 pseudotime analysis
BiocManager::install(version = "3.17")
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'HDF5Array',
                       'terra', 'ggrastr'))
install.packages("devtools")
devtools::install_github('cole-trapnell-lab/monocle3')
library(monocle3)
library(SeuratWrappers)

head(modules)



