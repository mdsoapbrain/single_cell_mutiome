## This script performs single-cell RNA-seq analysis using Seurat and hdWGCNA.
## It includes two main functions: 'determine_soft_power' and 'continue_analysis'.
## The script is designed to handle multiple cell types and performs various analyses
## such as differential module expression, motif analysis, and gene ontology enrichment.


# Library Imports
library(Seurat)
library(WGCNA)
library(hdWGCNA)
library(ggplot2)
library(patchwork)
library(dplyr)
library(TFBSTools)

# Dummy Variable Initialization for demonstration
# Replace these with your actual data
Drd2_data <- NULL
Drd1_data <- NULL
Astrocyte_data <- NULL
Olig_data <- NULL
Microglia_data <- NULL

# Function Definitions
determine_soft_power <- function(seurat_obj, group_name) {
  tryCatch({
  # Setup
  seurat_obj <- SetupForWGCNA(
    seurat_obj,
    gene_select = "fraction",
    fraction = 0.05,
    wgcna_name = group_name
  )
  
  # Construct metacells
  seurat_obj <- MetacellsByGroups(
    seurat_obj = seurat_obj,
    group.by = c("predict.id", "phenotype"),
    reduction = 'wnn.umap',
    k = 25,
    max_shared = 10,
    ident.group = 'predict.id'
  )
  
  # Normalize metacell expression matrix
  seurat_obj <- NormalizeMetacells(seurat_obj)
  seurat_obj <- SetDatExpr(
    seurat_obj,
    group_name = group_name,
    group.by = 'predict.id',
    assay = 'RNA',
    slot = 'data'
  )
  
  # Test different soft powers
  seurat_obj <- TestSoftPowers(seurat_obj, networkType = 'signed')
  
  # Plot the results
  plot_list <- PlotSoftPowers(seurat_obj)
  print(plot_list)
  
  return(seurat_obj)
  }, error = function(e) {
    print(paste("Error in determine_soft_power: ", e))
  })
}

continue_analysis <- function(seurat_obj, soft_power, group_name, phenotype1, phenotype2, output_prefix) {
  tryCatch({
  # Construct co-expression network
  seurat_obj <- ConstructNetwork(
    seurat_obj, 
    soft_power = soft_power,
    setDatExpr = FALSE,
    tom_name = group_name,overwrite_tom = TRUE
  )
  
  # Plot dendrogram
  pdf_file <- paste0(group_name, "_dendrogram.pdf")
  pdf(pdf_file)
  PlotDendrogram(seurat_obj, main = paste0(group_name, ' hdWGCNA Dendrogram'))
  dev.off()
  
  # Scale data
  seurat_obj <- ScaleData(seurat_obj, features=VariableFeatures(seurat_obj))
  
  # Compute all MEs in the full single-cell dataset
  seurat_obj <- ModuleEigengenes(seurat_obj, group.by.vars="phenotype")
  
  # Compute eigengene-based connectivity (kME)
  seurat_obj <- ModuleConnectivity(seurat_obj, group.by = 'predict.id', group_name = group_name)
  
  # Rename the modules
  seurat_obj <- ResetModuleNames(seurat_obj, new_name = group_name)
  
  # DME analysis comparing two groups
  group1 <- seurat_obj@meta.data %>% subset(phenotype == 'control') %>% rownames
  group2 <- seurat_obj@meta.data %>% subset(phenotype == 'ROT') %>% rownames
  
  DMEs <- FindDMEs(
    seurat_obj,
    barcodes1 = group1,
    barcodes2 = group2,
    test.use = 'wilcox',
    wgcna_name = group_name
  )
  head(DMEs)
  
  pdf(paste0(group_name, "_DME_Lollipop.pdf"), width=6, height=8)
  PlotDMEsLollipop(
    seurat_obj, 
    DMEs, 
    wgcna_name=group_name, 
    pvalue = "p_val_adj"
  )
  dev.off()
  pdf(paste0(group_name, "_DME_Volcano.pdf"), width=6, height=8)
  PlotDMEsVolcano_Simplified(
    seurat_obj,
    DMEs,
    wgcna_name = group_name
  )
  dev.off()
  #Enrich
  seurat_obj <- RunEnrichr(
    seurat_obj,
    dbs=dbs, # character vector of enrichr databases to test
    max_genes = 100 # number of genes per module to test
  )
  # enrichr dotplot
  
  pdf(paste0(group_name, "_Go_Enr.pdf"), width=6, height=8)
  EnrichrDotPlot(
    seurat_obj,
    mods = "all", # use all modules (this is the default behavior)
    database = "GO_Biological_Process_2021", # this has to be one of the lists we used above!!!
    n_terms=1 # number of terms for each module
  )
  dev.off()
  
  
  
  # Output
  modules <- GetModules(seurat_obj)
  module_file <- paste0(output_prefix, "_modules.csv")
  write.csv(modules[modules$color!='grey',], module_file, row.names = FALSE)
  
  enrich_df <- GetEnrichrTable(seurat_obj)
  enrich_file <- paste0(output_prefix, "_enrich_df.csv")
  write.csv(enrich_df, enrich_file, row.names = FALSE)
  
  # Motif analysis
  pfm_core <- TFBSTools::getMatrixSet(
    x = JASPAR2020,
    opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
  )
  
  seurat_obj <- dannymotif(
    seurat_obj,
    species_genome = 'BSgenome.Rnorvegicus.UCSC.rn7',
    pfm = pfm_core,
    EnsDb = EnsDb.Rnorvegicus.v105
  )
  
  seurat_obj <- OverlapModulesMotifs(seurat_obj)
  
  return(seurat_obj)
  }, error = function(e) {
    print(paste("Error in continue_analysis: ", e))
  })
}

# Main Code
tryCatch({
  Drd2_data <- determine_soft_power(Drd2_data, 'Drd2-MSN')
  Drd1_data <- determine_soft_power(Drd1_data, 'Drd1-MSN')
  Astrocyte_data <- determine_soft_power(Astrocyte_data, 'Astrocyte')
  Olig_data <- determine_soft_power(Olig_data, 'Olig')
  Microglia_data <- determine_soft_power(Microglia_data, 'Microglia')
  
  selected_soft_power <- 4
  
  Drd2_data <- continue_analysis(Drd2_data, 4, "Drd2-MSN", "control", "ROT", "Drd2")
  Drd1_data <- continue_analysis(Drd1_data, 4, "Drd1-MSN", "control", "ROT", "Drd1")
  Astrocyte_data <- continue_analysis(Astrocyte_data, 4, "Astrocyte", "control", "ROT", "Astrocyte")
  Olig_data <- continue_analysis(Olig_data, 4, "Olig", "control", "ROT", "Olig")
  Microglia_data <- continue_analysis(Microglia_data, 4, "Microglia", "control", "ROT", "Microglia")
  
  # DME analysis comparing two groups
  group1 <- Drd1_data@meta.data %>% subset(phenotype == 'control') %>% rownames
  group2 <- Drd1_data@meta.data %>% subset(phenotype == 'ROT') %>% rownames
  
  DMEs <- FindDMEs(
    Drd1_data,
    barcodes1 = group1,
    barcodes2 = group2,
    test.use = 'wilcox',
    wgcna_name = "Drd1-MSN"
  )
  PlotDMEsLollipop(
    Drd1_data, 
    DMEs, 
    wgcna_name='Drd1-MSN', 
    pvalue = "p_val_adj"
  )
  PlotDMEsVolcano(
    Drd1_data,
    DMEs,
    wgcna_name = 'Drd1-MSN'
  )
  
  
  pdf("Drd2_wgcna_GObp.pdf", width=6, height=8)
  EnrichrDotPlot(
    Drd2_data,
    mods = "all", # use all modules (this is the default behavior)
    database = "GO_Biological_Process_2021", # this has to be one of the lists we used above!!!
    n_terms=1 # number of terms for each module
  )
  dev.off()
  
  ##To3 3 TF in Drd1 HINFP, ZFP57, CENPB
  df <- GetMotifOverlap(Drd1_data)
  
  cur_df <- df %>% subset(tf == 'HINFP')
  plot_var <- 'odds_ratio'
  p <- cur_df %>%
    ggplot(aes(y=reorder(module, odds_ratio), x=odds_ratio)) +
    geom_bar(stat='identity', fill=cur_df$color) +
    geom_vline(xintercept = 1, linetype='dashed', color='gray') +
    geom_text(aes(label=Significance), color='black', size=3.5, hjust='center') +
    ylab('') +
    xlab("Odds Ratio") +
    ggtitle("HINFP overlap") +
    theme(
      plot.title = element_text(hjust = 0.5)
    )
  png('HINFP_motif_overlap_or.png', width=3, height=4, units='in', res=400)
  p
  dev.off()
  
  cur_df <- df %>% subset(tf == 'ZFP57')
  plot_var <- 'odds_ratio'
  p <- cur_df %>%
    ggplot(aes(y=reorder(module, odds_ratio), x=odds_ratio)) +
    geom_bar(stat='identity', fill=cur_df$color) +
    geom_vline(xintercept = 1, linetype='dashed', color='gray') +
    geom_text(aes(label=Significance), color='black', size=3.5, hjust='center') +
    ylab('') +
    xlab("Odds Ratio") +
    ggtitle("ZFP57 overlap") +
    theme(
      plot.title = element_text(hjust = 0.5)
    )
  png('ZFP57_motif_overlap_or.png', width=3, height=4, units='in', res=400)
  p
  dev.off()
  
  cur_df <- df %>% subset(tf == 'CENPB')
  plot_var <- 'odds_ratio'
  p <- cur_df %>%
    ggplot(aes(y=reorder(module, odds_ratio), x=odds_ratio)) +
    geom_bar(stat='identity', fill=cur_df$color) +
    geom_vline(xintercept = 1, linetype='dashed', color='gray') +
    geom_text(aes(label=Significance), color='black', size=3.5, hjust='center') +
    ylab('') +
    xlab("Odds Ratio") +
    ggtitle("CENPB overlap") +
    theme(
      plot.title = element_text(hjust = 0.5)
    )
  png('CENPB_motif_overlap_or.png', width=3, height=4, units='in', res=400)
  p
  dev.off()
  
  
  
  
  
  
  
  
  
  ###Astrocyte######
  
  # DME analysis comparing two groups
  group1 <- Astrocyte_data@meta.data %>% subset(phenotype == 'control') %>% rownames
  group2 <- Astrocyte_data@meta.data %>% subset(phenotype == 'ROT') %>% rownames
  modules <- GetModules(Astrocyte_data, "Astrocyte") %>% subset(module !="grey")
  
  DMEs <- FindDMEs(
    Astrocyte_data,
    barcodes1 = group1,
    barcodes2 = group2,
    test.use = 'wilcox',
    wgcna_name = "Astrocyte"
  )
  PlotDMEsLollipop(
    Astrocyte_data, 
    DMEs, 
    wgcna_name='Astrocyte', 
    pvalue = "p_val_adj"
  )
  PlotDMEsVolcano_Simplified(
    Astrocyte_data,
    DMEs,
    wgcna_name = 'Astrocyte'
  )
  
  pdf("Astrocyte_wgcna_GObp.pdf", width=6, height=8)
  EnrichrDotPlot(
    Astrocyte_data,
    mods = "all", # use all modules (this is the default behavior)
    database = "GO_Biological_Process_2021", # this has to be one of the lists we used above!!!
    n_terms=1 # number of terms for each module
  )
  dev.off()
  
  
  
  ##To3 3 TF in astrocyte ISX, Shox2, LHX9
  df <- GetMotifOverlap(Astrocyte_data)
  
  cur_df <- df %>% subset(tf == 'ISX')
  plot_var <- 'odds_ratio'
  p <- cur_df %>%
    ggplot(aes(y=reorder(module, odds_ratio), x=odds_ratio)) +
    geom_bar(stat='identity', fill=cur_df$color) +
    geom_vline(xintercept = 1, linetype='dashed', color='gray') +
    geom_text(aes(label=Significance), color='black', size=3.5, hjust='center') +
    ylab('') +
    xlab("Odds Ratio") +
    ggtitle("ISX overlap") +
    theme(
      plot.title = element_text(hjust = 0.5)
    )
  png('ISX_motif_overlap_or.png', width=3, height=4, units='in', res=400)
  p
  dev.off()
  
  cur_df <- df %>% subset(tf == 'Shox2')
  plot_var <- 'odds_ratio'
  p <- cur_df %>%
    ggplot(aes(y=reorder(module, odds_ratio), x=odds_ratio)) +
    geom_bar(stat='identity', fill=cur_df$color) +
    geom_vline(xintercept = 1, linetype='dashed', color='gray') +
    geom_text(aes(label=Significance), color='black', size=3.5, hjust='center') +
    ylab('') +
    xlab("Odds Ratio") +
    ggtitle("Shox2 overlap") +
    theme(
      plot.title = element_text(hjust = 0.5)
    )
  png('Shox2_motif_overlap_or.png', width=3, height=4, units='in', res=400)
  p
  dev.off()
  
  cur_df <- df %>% subset(tf == 'LHX9')
  plot_var <- 'odds_ratio'
  p <- cur_df %>%
    ggplot(aes(y=reorder(module, odds_ratio), x=odds_ratio)) +
    geom_bar(stat='identity', fill=cur_df$color) +
    geom_vline(xintercept = 1, linetype='dashed', color='gray') +
    geom_text(aes(label=Significance), color='black', size=3.5, hjust='center') +
    ylab('') +
    xlab("Odds Ratio") +
    ggtitle("LHX9 overlap") +
    theme(
      plot.title = element_text(hjust = 0.5)
    )
  png('LHX9_motif_overlap_or.png', width=3, height=4, units='in', res=400)
  p
  dev.off()
  ########Olig########
  # DME analysis comparing two groups
  group1 <- Olig_data@meta.data %>% subset(phenotype == 'control') %>% rownames
  group2 <- Olig_data@meta.data %>% subset(phenotype == 'ROT') %>% rownames
  modules <- GetModules(Olig_data, "Olig") %>% subset(module !="grey")
  
  DMEs <- FindDMEs(
    Olig_data,
    barcodes1 = group1,
    barcodes2 = group2,
    test.use = 'wilcox',
    wgcna_name = "Olig"
  )
  PlotDMEsLollipop(
    Olig_data, 
    DMEs, 
    wgcna_name='Olig', 
    pvalue = "p_val_adj"
  )
  PlotDMEsVolcano_Simplified(
    Olig_data,
    DMEs,
    wgcna_name = 'Olig'
  )
  
  pdf("Olig_wgcna_GObp.pdf", width=6, height=8)
  EnrichrDotPlot(
    Olig_data,
    mods = "all", # use all modules (this is the default behavior)
    database = "GO_Biological_Process_2021", # this has to be one of the lists we used above!!!
    n_terms=1 # number of terms for each module
  )
  dev.off()
  
  
  
  ##To3 3 TF in Olig E2F2, ZBED1, E2F4
  df <- GetMotifOverlap(Olig_data)
  
  cur_df <- df %>% subset(tf == 'E2F2')
  plot_var <- 'odds_ratio'
  p <- cur_df %>%
    ggplot(aes(y=reorder(module, odds_ratio), x=odds_ratio)) +
    geom_bar(stat='identity', fill=cur_df$color) +
    geom_vline(xintercept = 1, linetype='dashed', color='gray') +
    geom_text(aes(label=Significance), color='black', size=3.5, hjust='center') +
    ylab('') +
    xlab("Odds Ratio") +
    ggtitle("E2F2 overlap") +
    theme(
      plot.title = element_text(hjust = 0.5)
    )
  png('E2F2_motif_overlap_or.png', width=3, height=4, units='in', res=400)
  p
  dev.off()
  
  cur_df <- df %>% subset(tf == 'ZBED1')
  plot_var <- 'odds_ratio'
  p <- cur_df %>%
    ggplot(aes(y=reorder(module, odds_ratio), x=odds_ratio)) +
    geom_bar(stat='identity', fill=cur_df$color) +
    geom_vline(xintercept = 1, linetype='dashed', color='gray') +
    geom_text(aes(label=Significance), color='black', size=3.5, hjust='center') +
    ylab('') +
    xlab("Odds Ratio") +
    ggtitle("ZBED1 overlap") +
    theme(
      plot.title = element_text(hjust = 0.5)
    )
  png('ZBED1_motif_overlap_or.png', width=3, height=4, units='in', res=400)
  p
  dev.off()
  
  cur_df <- df %>% subset(tf == 'E2F4')
  plot_var <- 'odds_ratio'
  p <- cur_df %>%
    ggplot(aes(y=reorder(module, odds_ratio), x=odds_ratio)) +
    geom_bar(stat='identity', fill=cur_df$color) +
    geom_vline(xintercept = 1, linetype='dashed', color='gray') +
    geom_text(aes(label=Significance), color='black', size=3.5, hjust='center') +
    ylab('') +
    xlab("Odds Ratio") +
    ggtitle("E2F4 overlap") +
    theme(
      plot.title = element_text(hjust = 0.5)
    )
  png('E2F4_motif_overlap_or.png', width=3, height=4, units='in', res=400)
  p
  dev.off()
  ########Microglia#########
  # DME analysis comparing two groups
  group1 <- Microglia_data@meta.data %>% subset(phenotype == 'control') %>% rownames
  group2 <- Microglia_data@meta.data %>% subset(phenotype == 'ROT') %>% rownames
  modules <- GetModules(Microglia_data, "Microglia") %>% subset(module !="grey")
  
  DMEs <- FindDMEs(
    Microglia_data,
    barcodes1 = group1,
    barcodes2 = group2,
    test.use = 'wilcox',
    wgcna_name = "Microglia"
  )
  PlotDMEsLollipop(
    Microglia_data, 
    DMEs, 
    wgcna_name='Microglia', 
    pvalue = "p_val_adj"
  )
  PlotDMEsVolcano_Simplified(
    Microglia_data,
    DMEs,
    wgcna_name = 'Microglia'
  )
  
  pdf("Microglia_wgcna_GObp.pdf", width=6, height=8)
  EnrichrDotPlot(
    Microglia_data,
    mods = "all", # use all modules (this is the default behavior)
    database = "GO_Biological_Process_2021", # this has to be one of the lists we used above!!!
    n_terms=1 # number of terms for each module
  )
  dev.off()
  
  
  ##To3 3 TF in Microglia E2F3, E2F1, E2F2
  df <- GetMotifOverlap(Microglia_data)
  
  cur_df <- df %>% subset(tf == 'E2F3')
  plot_var <- 'odds_ratio'
  p <- cur_df %>%
    ggplot(aes(y=reorder(module, odds_ratio), x=odds_ratio)) +
    geom_bar(stat='identity', fill=cur_df$color) +
    geom_vline(xintercept = 1, linetype='dashed', color='gray') +
    geom_text(aes(label=Significance), color='black', size=3.5, hjust='center') +
    ylab('') +
    xlab("Odds Ratio") +
    ggtitle("E2F3 overlap") +
    theme(
      plot.title = element_text(hjust = 0.5)
    )
  png('E2F3_motif_overlap_or.png', width=3, height=4, units='in', res=400)
  p
  dev.off()
  
  cur_df <- df %>% subset(tf == 'E2F1')
  plot_var <- 'odds_ratio'
  p <- cur_df %>%
    ggplot(aes(y=reorder(module, odds_ratio), x=odds_ratio)) +
    geom_bar(stat='identity', fill=cur_df$color) +
    geom_vline(xintercept = 1, linetype='dashed', color='gray') +
    geom_text(aes(label=Significance), color='black', size=3.5, hjust='center') +
    ylab('') +
    xlab("Odds Ratio") +
    ggtitle("E2F1 overlap") +
    theme(
      plot.title = element_text(hjust = 0.5)
    )
  png('E2F1_motif_overlap_or.png', width=3, height=4, units='in', res=400)
  p
  dev.off()
  
  cur_df <- df %>% subset(tf == 'E2F2')
  plot_var <- 'odds_ratio'
  p <- cur_df %>%
    ggplot(aes(y=reorder(module, odds_ratio), x=odds_ratio)) +
    geom_bar(stat='identity', fill=cur_df$color) +
    geom_vline(xintercept = 1, linetype='dashed', color='gray') +
    geom_text(aes(label=Significance), color='black', size=3.5, hjust='center') +
    ylab('') +
    xlab("Odds Ratio") +
    ggtitle("E2F2 overlap") +
    theme(
      plot.title = element_text(hjust = 0.5)
    )
  png('E2F2_motif_overlap_or.png', width=3, height=4, units='in', res=400)
  p
  dev.off()
  
}, error = function(e) {
  print(paste("Error in main code: ", e))
})
