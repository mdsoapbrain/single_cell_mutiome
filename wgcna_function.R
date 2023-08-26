EnrichrDotPlot<-function (seurat_obj, database, mods = "all", n_terms = 3, break_ties = TRUE, 
     logscale = TRUE, wgcna_name = NULL, ...) 
 {
     if (is.null(wgcna_name)) {
         wgcna_name <- seurat_obj@misc$active_wgcna
     }
     modules <- GetModules(seurat_obj, wgcna_name)
     if (mods == "all") {
         mods <- levels(modules$module)
         mods <- mods[mods != "grey"]
     }
     enrichr_df <- GetEnrichrTable(seurat_obj, wgcna_name)
     mod_colors <- dplyr::select(modules, c(module, color)) %>% distinct
     enrichr_df$color <- mod_colors[match(enrichr_df$module, mod_colors$module), 
         "color"]
     wrapText <- function(x, len) {
         sapply(x, function(y) paste(strwrap(y, len), collapse = "\n"), 
             USE.NAMES = FALSE)
     }
     plot_df <- enrichr_df %>% subset(db == database & module %in% 
         mods) %>% group_by(module) %>% top_n(n_terms, wt = Combined.Score)
     if (break_ties) {
         plot_df <- do.call(rbind, lapply(plot_df %>% group_by(module) %>% 
             group_split, function(x) {
             x[sample(n_terms), ]
         }))
     }
     plot_df$Term <- wrapText(plot_df$Term, 45)
     plot_df$module <- factor(as.character(plot_df$module), levels = levels(modules$module))
     plot_df <- arrange(plot_df, module)
     plot_df$Term <- factor(as.character(plot_df$Term), levels = unique(as.character(plot_df$Term)))
     if (logscale) {
         plot_df$Combined.Score <- log(plot_df$Combined.Score)
         lab <- "Enrichment\nlog(combined score)"
         x <- 0.2
     }
     else {
         lab <- "Enrichment\n(combined score)"
         x <- 5
     }
     p <- plot_df %>% ggplot(aes(x = module, y = Term)) + geom_point(aes(size = Combined.Score), 
         color = plot_df$color) + RotatedAxis() + ylab("") + xlab("") + 
         labs(size = lab) + scale_y_discrete(limits = rev) + ggtitle(database) + 
         theme(plot.title = element_text(hjust = 0.5), axis.line.x = element_blank(), 
             axis.line.y = element_blank(), panel.border = element_rect(colour = "black", 
                 fill = NA, size = 1))
     p
 }



dannymotif <- function (seurat_obj, species_genome, pfm, EnsDb, wgcna_name = 'Drd2_data') 
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



PlotDMEsVolcano_Simplified <- function(seurat_obj, DMEs, plot_labels = TRUE, 
                                       mod_point_size = 4, label_size = 4, 
                                       show_cutoff = TRUE, wgcna_name = NULL) {
  
  modules <- GetModules(seurat_obj, wgcna_name) %>% subset(module != "grey") %>% 
    mutate(module = droplevels(module))
  
  module_colors <- modules %>% dplyr::select(c(module, color)) %>% distinct
  
  mod_colors <- module_colors$color
  names(mod_colors) <- as.character(module_colors$module)
  
  DMEs$anno <- ifelse(DMEs$p_val_adj < 0.05, DMEs$module, "")
  
  p <- DMEs %>% ggplot(aes(x = avg_log2FC, y = -log10(p_val_adj), 
                           fill = module, color = module))
  
  if (show_cutoff) {
    p <- p + geom_vline(xintercept = 0, linetype = "dashed", 
                        color = "grey75", alpha = 0.8) + 
      geom_rect(data = DMEs[1, ], 
                aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -log10(0.05)), 
                fill = "grey75", alpha = 0.8, color = NA)
  }
  
  p <- p + geom_point(size = mod_point_size, pch = 21, color = "black")
  
  if (plot_labels) {
    p <- p + geom_text_repel(aes(label = anno), color = "black", 
                             min.segment.length = 0, max.overlaps = Inf, size = label_size)
  }
  
  p <- p + scale_fill_manual(values = mod_colors) + 
    scale_color_manual(values = mod_colors) + 
    xlab(bquote("Average log"[2] ~ "(Fold Change)")) + 
    ylab(bquote("-log"[10] ~ "(Adj. P-value)")) + 
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1), 
          panel.grid.major = element_blank(), axis.line = element_blank(), 
          plot.title = element_text(hjust = 0.5), legend.position = "bottom") + 
    NoLegend()
  
  return(p + NoLegend())
}
