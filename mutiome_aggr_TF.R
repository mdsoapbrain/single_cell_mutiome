# packages for TF motif analysis
library(JASPAR2020)
library(motifmatchr)
library(TFBSTools)
library(EnsDb.Rnorvegicus.v105)
library(GenomicRanges)





dannymotif <- function (seurat_obj, species_genome, pfm, EnsDb, wgcna_name = NULL) 
{
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  motif_df <- data.frame(motif_name = purrr::map(1:length(pfm), 
                                                 function(i) {
                                                   pfm[[i]]@name
                                                 }) %>% unlist, motif_ID = purrr::map(1:length(pfm), function(i) {
                                                   pfm[[i]]@ID
                                                 }) %>% unlist)
  gene.promoters <- ensembldb::promoters(EnsDb, filter = ~gene_biotype == 
                                           "protein_coding") %>% subset(seqnames %in% c(1:100))
  gene.coords <- ensembldb::genes(EnsDb, filter = ~gene_biotype == 
                                    "protein_coding") %>% subset(seqnames %in% c(1:100))
  gene.promoters$symbol <- gene.coords$symbol[match(gene.promoters$gene_id, 
                                                    names(gene.coords))]
  gene.promoters <- keepSeqlevels(gene.promoters, value = levels(droplevels(seqnames(gene.promoters))))
  old_levels <- levels(seqnames(gene.promoters))
  new_levels <- ifelse(old_levels %in% c("X", "Y"), old_levels, paste0("chr", old_levels))
  gene.promoters <- renameSeqlevels(gene.promoters, old_levels)
  genome(seqinfo(gene.promoters)) <- species_genome
  my_promoters <- GRanges(seqnames = paste0("chr", droplevels(seqnames(gene.promoters))), 
                          IRanges(start = start(gene.promoters), end = end(gene.promoters)), 
                          symbol = gene.promoters$symbol, genome = species_genome)
  standard_chromosomes <- paste0("chr", c(1:20))
  my_promoters <- keepSeqlevels(my_promoters, standard_chromosomes, pruning.mode = "coarse")
  
  print("Matching motifs...")
  motif_ix <- motifmatchr::matchMotifs(pfm, my_promoters, genome = species_genome)
  tf_match <- motifmatchr::motifMatches(motif_ix)
  rownames(tf_match) <- my_promoters$symbol
  colnames(tf_match) <- motif_df$motif_name
  gene_list <- rownames(seurat_obj)
  gene_list <- gene_list[gene_list %in% rownames(tf_match)]
  tf_match <- tf_match[gene_list, ]
  print("Getting TF target genes...")
  tfs <- motif_df$motif_name
  tf_targets <- list()
  n_targets <- list()
  for (cur_tf in tfs) {
    tf_targets[[cur_tf]] <- names(tf_match[, cur_tf][tf_match[, 
                                                              cur_tf]])
    n_targets[[cur_tf]] <- length(tf_targets[[cur_tf]])
  }
  n_targets <- unlist(n_targets)
  motif_df$n_targets <- n_targets
  seurat_obj <- SetMotifMatrix(seurat_obj, tf_match)
  seurat_obj <- SetMotifs(seurat_obj, motif_df)
  seurat_obj <- SetMotifTargets(seurat_obj, tf_targets)
  seurat_obj <- SetPFMList(seurat_obj, pfm)
  seurat_obj
}

#BiocManager::install("BSgenome.Rnorvegicus.UCSC.rn5.masked")
# get the pfm from JASPAR2020 using TFBSTools
pfm_core <- TFBSTools::getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

data <- dannymotif(
  data,
  species_genome = 'BSgenome.Rnorvegicus.UCSC.rn5.masked',
  pfm = pfm_core,
  EnsDb = EnsDb.Rnorvegicus.v105
)

dim(GetMotifMatrix(data))
# TF target genes
target_genes <- GetMotifTargets(data)
###########overlap###########
data <- SetupForWGCNA(
  data,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "tutorial" # the name of the hdWGCNA experiment
)
# construct metacells 
data <- MetacellsByGroups(
  seurat_obj = data,
  group.by = "predict.id",
  k = 50,
  target_metacells=250,
  ident.group = 'predict.id',
  min_cells=0,
  max_shared=5,
)
data <- NormalizeMetacells(data)
# get list of cell populations
groups <- unique(GetMetacellObject(data)$predict.id)

# set up gene expression matrix
##
data <- SetDatExpr(
  data,
  group.by='predict.id',
  group_name = 'Drd2-MSN',# try Drd2 first, group_name = groups
  use_metacells=TRUE,
  slot = 'data',
  assay = 'RNA'
)
# test soft power parameter
data <- TestSoftPowers(data)

# construct the co-expression network
data <- ConstructNetwork(
  data, 
  tom_name='TF', 
  overwrite_tom=TRUE
)

# compute module eigengenes & connectivity
data <- ModuleEigengenes(data)
data <- ModuleConnectivity(data)

# get the module assignment table:
modules <- GetModules(data)

# overlap between modules & TF target genes:
data<- OverlapModulesMotifs(data)

# look at the overlap data
head(GetMotifOverlap(data))
# plot the top TFs overlapping with
MotifOverlapBarPlot(
  data,
  #motif_font = 'xkcd_regular',
  outdir = 'MotifOverlaps_drd2/',
  plot_size=c(5,6)
)

data <- MotifTargetScore(
  data,
  method='Seurat'
)

data_drd2 <- GetMotifOverlap(data)

write_csv(data,'drds_TF_module_overlap.csv')
cur_df <- data %>% subset(tf == 'HES1')

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


png(paste0('', 'HES1_motif_overlap_or.png'), width=3, height=4, units='in', res=400)
p
dev.off()

ModuleTFNetwork(
  data,
  edge.alpha=0.75,
  cor_thresh = 0.25,
  tf_name = "Hes1",
  tf_gene_name = "Apoe",
  tf_x = 0,
  tf_y = 0
)



#############################
HES1_df<-as.data.frame(target_genes[["HES1"]])
modules_gene <- as.data.frame(modules[modules$color == 'turquoise', c("gene_name")])
colnames(modules_gene)<-'gene'

gene_in_list<-HES1_df[HES1_df$`target_genes[["HES1"]]`%in%modules_gene$`gene`,]

gene_in_list

