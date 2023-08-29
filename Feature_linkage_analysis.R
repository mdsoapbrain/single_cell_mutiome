rm(list = ls())
setwd("/Volumes/Seagate/multiome_cp")
#.libPaths("/Volumes/Seagate/multiome_cp/lib/R/library/")
set.seed(1234)
#install.packages("GenomeInfoDb")
library(EnsDb.Rnorvegicus.v105)
library(harmony)# easy aborted...
library(Signac)
library(Seurat)
library(BSgenome.Rnorvegicus.UCSC.rn7)
library(SeuratDisk)
library(dplyr)
library(patchwork)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(SeuratObject)
library(cowplot)
library(reticulate)
library(GenomeInfoDb)

BiocManager::install("rtracklayer", force = TRUE)
library(rtracklayer)
part1<-import("r3_feature_linkage.bedpe", format="bedpe")
part1 <- read.table("ROT_feature_linkage.bedpe", header = FALSE, sep = "\t")

a10<-as.data.frame(part1)

a10<-a10[which(a10$V13!="peak-peak"),]
{
  a10$seqnames <- a10$V1
  start <- as.integer(ifelse(a10$V13 == "peak-gene", (a10$V3 + a10$V2)/2, (a10$V6 + a10$V5)/2))
  end <- ifelse(a10$V13 == "gene-peak", a10$V2, a10$V5)
  a10$start <- ifelse(start < end, start, end)
  a10$end <- ifelse(start < end, end, start)
  a10$peak.seqnames <- paste0("chr", ifelse(a10$V13 == "peak-gene", a10$V1, a10$V4))
  a10$peak.start <- ifelse(a10$V13 == "peak-gene", a10$V2, a10$V5)
  a10$peak.end <- ifelse(a10$V13 == "peak-gene", a10$V3, a10$V6)
  a <- matrix(unlist(strsplit(a10$V7, "><")), ncol = 2, byrow = T)
  gene <- ifelse(a10$V13 == "peak-gene", noquote(gsub(">", "", a[,2])), noquote(gsub("<", "", a[,1])))
  peak <- ifelse(a10$V13 == "gene-peak", noquote(gsub(">", "", a[,2])), noquote(gsub("<", "", a[,1])))
  a10$gene <- gene
  a10$peak <- peak
  a10$gene.strand <- "*"
  a10$peak.strand <- "*"
  a10$gene.seqnames <- paste0("chr", ifelse(a10$V13 == "gene-peak", a10$V1, a10$V4))
  a10$gene.start <- ifelse(a10$V13 == "gene-peak", a10$V2, a10$V5)
  a10$gene.end <- ifelse(a10$V13 == "gene-peak", a10$V3, a10$V6)
  a10$Peak_pos <- paste0(a10$peak.seqnames, "-", a10$peak.start, "-", a10$peak.end)
}# # #  format a10 ( ad)

a10<-a10[which(abs(a10$V8)>0.2),]

annot_peaks<-read.csv("CTpeaks_annotated.csv")
atac_peak_df <- read.table("ROT_atac_peak_annotation.tsv", header = TRUE, sep = "\t")

annot_peaks$X<-NULL
annot_peaks$peaks <-paste0("chr", annot_peaks$seqnames, "-", annot_peaks$start, "-", annot_peaks$end)
atac_peak_df$peaks <-paste0("chr", atac_peak_df$chrom, "-", atac_peak_df$start, "-", atac_peak_df$end)

df <- merge(annot_peaks, atac_peak_df, by = "peaks", all.x = TRUE)

df<-df[!is.na(df$distance),]
a10_ct<-merge(a10, df, by.x="Peak_pos", by.y="peaks")
a10_ct$index <- 1:nrow(a10_ct)
a10_ct$link<-paste0(a10_ct$gene.x,"-", a10_ct$index)
a10_qval<-data.frame(Qval=a10_ct$NA., link=a10_ct$link)
##########################
##########################
part2 <- read.table("Control_feature_linkage.bedpe", header = FALSE, sep = "\t")

c10<-as.data.frame(part2)

c10<-c10[which(c10$V13!="peak-peak"),]
{
  c10$seqnames <- c10$V1
  start <- as.integer(ifelse(c10$V13 == "peak-gene", (c10$V3 + c10$V2)/2, (c10$V6 + c10$V5)/2))
  end <- ifelse(c10$V13 == "gene-peak", c10$V2, c10$V5)
  c10$start <- ifelse(start < end, start, end)
  c10$end <- ifelse(start < end, end, start)
  c10$peak.seqnames <- paste0("chr", ifelse(c10$V13 == "peak-gene", c10$V1, c10$V4))
  c10$peak.start <- ifelse(c10$V13 == "peak-gene", c10$V2, c10$V5)
  c10$peak.end <- ifelse(c10$V13 == "peak-gene", c10$V3, c10$V6)
  a <- matrix(unlist(strsplit(c10$V7, "><")), ncol = 2, byrow = T)
  gene <- ifelse(c10$V13 == "peak-gene", noquote(gsub(">", "", a[,2])), noquote(gsub("<", "", a[,1])))
  peak <- ifelse(c10$V13 == "gene-peak", noquote(gsub(">", "", a[,2])), noquote(gsub("<", "", a[,1])))
  c10$gene <- gene
  c10$peak <- peak
  c10$gene.strand <- "*"
  c10$peak.strand <- "*"
  c10$gene.seqnames <- paste0("chr", ifelse(c10$V13 == "gene-peak", c10$V1, c10$V4))
  c10$gene.start <- ifelse(c10$V13 == "gene-peak", c10$V2, c10$V5)
  c10$gene.end <- ifelse(c10$V13 == "gene-peak", c10$V3, c10$V6)
  c10$Peak_pos <- paste0(c10$peak.seqnames, "-", c10$peak.start, "-", c10$peak.end)
}

# # #  format c10 ( ad)

c10<-c10[which(abs(c10$V8)>0.2),]

#annot_peaks<-read.csv("CTpeaks_annotated.csv")
c10_atac_peak_df <- read.table("Control_atac_peak_annotation.tsv", header = TRUE, sep = "\t")

#annot_peaks$X<-NULL
#annot_peaks$peaks <-paste0("chr", annot_peaks$seqnames, "-", annot_peaks$start, "-", annot_peaks$end)
c10_atac_peak_df$peaks <-paste0("chr", c10_atac_peak_df$chrom, "-", c10_atac_peak_df$start, "-", c10_atac_peak_df$end)

c10_df <- merge(annot_peaks, c10_atac_peak_df, by = "peaks", all.x = TRUE)

c10_df<-c10_df[!is.na(c10_df$distance),]
c10_df$chrom<-NULL
c10_df$start.y<-NULL
c10_df$end.y<-NULL
c10_df$distance<-NULL
c10_ct<-merge(c10, c10_df, by.x="Peak_pos", by.y="peaks")
c10_ct$index <- 1:nrow(c10_ct)
c10_ct$link<-paste0(c10_ct$gene.x,"-", c10_ct$index)
c10_qval<-data.frame(Qval=c10_ct$NA., link=c10_ct$link)







################################
##Common
all<-merge(a10_ct, c10_ct, by=c("gene.x","peak.seqnames","peak.start","peak.end"), all=T)
all$group<-ifelse(is.na(all$V8.x)==T, "Ctrl",
                  ifelse(is.na(all$V8.y)==T, "ROT","common"))

all$CT<-ifelse(is.na(all$V8.x)==T, all$peak_called_in.y, all$peak_called_in.x)
################################ make bed file for LDSC#############
interest_df <- all[,c('Peak_pos.x','Peak_pos.y','group','CT')]

write_bed <- function(df, cell_type, group, filename) {
  df <- df[df$group == group, ]
  
  # For ROT and Ctrl, peak positions can be found in either the Peak_pos.x or Peak_pos.y columns
  if(group == "ROT") {
    df$Peak_pos <- ifelse(is.na(df$Peak_pos.x), df$Peak_pos.y, df$Peak_pos.x)
  } else if(group == "Ctrl") {
    df$Peak_pos <- ifelse(is.na(df$Peak_pos.y), df$Peak_pos.x, df$Peak_pos.y)
  } else {
    df$Peak_pos <- df$Peak_pos.x
  }
  
  # Split the Peak_pos column values into three parts and extract chr, start, and end
  df$chr <- sapply(strsplit(as.character(df$Peak_pos), "-"), "[[", 1)
  df$start <- sapply(strsplit(as.character(df$Peak_pos), "-"), "[[", 2)
  df$end <- sapply(strsplit(as.character(df$Peak_pos), "-"), "[[", 3)
  
  # Check if the CT column contains the specified cell types
  df_filtered <- df[grepl(cell_type, df$CT), ]
  
  write.table(df_filtered[, c("chr", "start", "end")],
              paste0("./LDSC_aggr_bed/", filename, ".bed"),
              quote = F, row.names = F, col.names = F, sep = "\t")
}

cell_types <- c("Astrocyte", "Microglia", "OPC", "Olig", "Drd1-MSN","Drd2-MSN", "Glutamatergic", "Interneuron")
groups <- c("common", "Ctrl", "ROT")

for(cell_type in cell_types) {
  for(group in groups) {
    write_bed(interest_df, cell_type, group, paste0(cell_type, "_", group))
  }
}

write.csv(all, 'all.csv')




#####################




t<-table(all$CT, all$group)
t<-t[order(-t[,2]),]
write.table(t, "Entire_ct_link_tablefiltered.txt", sep="\t")

all$Olig      <-grepl("Olig",all$CT)
all$OPC       <-grepl("OPC",all$CT)
all$GABA_neuron      <-grepl("GABA-neuron",      all$CT)
all$Microglia      <-grepl("Microglia",      all$CT)
all$Glut_neuron      <-grepl("Glut-neuron",      all$CT)
all$Astrocyte             <-grepl("Astrocyte",            all$CT)
all$DA_neuron        <-grepl("DA-neuron",       all$CT)


tab<-rbind(table(all$Olig ,  all$group   ),
           table(all$OPC  ,  all$group   ),
           table(all$GABA_neuron ,  all$group   ),
           table(all$Microglia ,  all$group   ),
           table(all$Glut_neuron ,  all$group   ),
           table(all$Astrocyte, all$group),
           table(all$DA_neuron       , all$group  ))
tab2<-tab[which(rownames(tab)==TRUE),]
rownames(tab2)<-c("Olig","OPC","GABA-neuron",
                  "Microglia","Glut-neuron",
                  "Astrocyte","DA-neuron")

write.table(tab2, "Melt_ct_link_table.txt", sep="\t", quote=F)

# First, transform the data into a "long" format suitable for ggplot
library(reshape2)
tab2_long <- melt(tab2)

library(ggplot2)
ggplot(tab2_long, aes(x=Var2, y=value, fill=Var1)) + 
  geom_bar(stat='identity', position='stack') + 
  labs(x='Group', y='Count', fill='Cell Type') +
  theme_minimal() +
  scale_fill_brewer(palette='Set1') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




###################
install.packages("UpSetR")
library(UpSetR)
upset(fromExpression(tab2), nintersects=12, nsets=6,
      sets.bar.color=c(
  "cornflowerblue","seagreen3","coral3","darkgoldenrod1","mediumorchid3","firebrick"),
  main.bar.color="blue", 
  text.scale=1.3,  mainbar.y.label="Ctrl Links",
  mainbar.y.max=60000, show.numbers=F, keep.order=T, order.by=c("freq","degree"))
###################
###################
######Fig 3A
#c2.df<-as.data.frame(c2)
write_csv(all, 'all.csv')
all<-read.csv('all.csv')
c2<-all
c2$V8<-ifelse(c2$group=="ROT", c2$V8.x, #score=V8
                    ifelse(c2$group=="Ctrl", c2$V8.y, (c2$V8.x+ c2$V8.y)/2))
c2<-c2[order(-c2$V8),]
c2<-c2[!duplicated(c2$gene.x),]
c2$peak<-paste0(c2$peak.seqnames,"-",c2$peak.start,"-", c2$peak.end)
c2<-c2[!duplicated(c2$peak),]
c2$dir<-ifelse(c2$V8<0, "down","up")
c2<-c2[which(c2$dir=="up"),]

data<-readRDS("mutiome_final.rds")

acc<-as.data.frame(AverageExpression(data, group.by=c("predict.id","phenotype"),
                                     assay="CTpeaks",
                                     features=c2$peak))
gex<-as.data.frame(AverageExpression(data, group.by=c("predict.id","phenotype"),
                                     assay="RNA",
                                     features=c2$gene.x))


acc.scale<-t(apply(acc,1, scale))
gex.scale<-t(apply(gex,1,scale))
colnames(acc.scale)<-colnames(acc)
colnames(gex.scale)<-colnames(gex)

acc.scale<-acc.scale[complete.cases(gex.scale),]
gex.scale<-gex.scale[complete.cases(gex.scale),]
c2<-c2[which(c2$gene.x %in% rownames(gex.scale)),]
write.csv(c2,"c2.csv" )
# 将acc.scale保存为文本文件
output_file <- "acc_scale.txt"
write.table(acc.scale, file = output_file, sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
output_file <- "gex_scale.txt"
write.table(gex.scale, file = output_file, sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)


BiocManager::install("ComplexHeatmap")
install.packages("colorRamp2")
library(ComplexHeatmap)
library(colorRamp2)
ha2<-HeatmapAnnotation(celltype=rep(c("Olig","Astrocyte","Microglia","DA-neuron","OPC","GABA-neuron", "Glut-neuron"), each=2),
                       col= list(celltype=c("Olig"="darkgoldenrod1","Astrocyte"="cornflowerblue","Microglia"="seagreen3",
                                            "DA-neuron"="mediumorchid3","OPC"="coral3","GABA-neuron"="firebrick","Glut-neuron"="blue" )), show_legend=F)
ha<-HeatmapAnnotation(Phenotype=rep(c("PD","Ctrl"), 7),
                      col=list(Phenotype=c("PD"="red","Ctrl"="blue")),
                      show_legend=FALSE)
ha<-c( ha, ha2)

ra <- rowAnnotation(Group = c2$group,  `Score` = c2$score, 
                    col = list(Group = c("common"="grey50", "PD"="red", "Ctrl"="blue"), 
                               Dir = c("up"="yellow", "down"="purple"), 
                               `Score` = colorRamp2(c(0,0.25,0.75,1),c("white","blue","green","red"))))

ht1 = Heatmap(acc.scale, cluster_rows = T, show_row_dend = F,
              cluster_columns = F, col = colorRamp2(c(quantile(acc.scale, probs = 0.05), 0, 2), c("#290230", "#CC4678FF", "#F0F921FF" )),
              name = "Accessibility Z",  show_column_names = F, show_row_names = F, column_title = NULL, top_annotation = ha,
              row_dend_reorder = T, use_raster = T, raster_quality = 5,
              row_split = c2$group)

print(length(c2$group))
print(nrow(acc.scale))

install.packages("viridis")
install.packages("magick")
library(magick)
library(viridis)
ht2=Heatmap(gex.scale, cluster_rows=F,show_row_dend=F,
            cluster_columns=F,col=colorRamp2(c(quantile(gex.scale, probs=0.05, na.rm=T),0,quantile(gex.scale, probs=0.95, na.rm=T)),viridis(3)),
            name="Expression  Z",  show_column_names=F, show_row_names=F,
            column_title=NULL,row_names_gp = gpar(fontsize = 3), top_annotation=ha,
            right_annotation=ra,  row_title_gp=gpar(fontsize=0),
            use_raster=T, raster_quality=5)
capabilities("cairo")
output_file <- "/Volumes/Seagate/multiome_cp/Acc_GEX_link_filtered.png"
library(ggplot2)
library(gridExtra)
library(cowplot)
pdf("Acc_GEX_link_filtered.pdf")
ht1+ht2
dev.off()

draw(ht1)





library(ggplot2)
library(reshape2)
library(viridis)


acc.melt <- melt(acc.scale)
gex.melt <- melt(gex.scale)

names(acc.melt) <- c("gene", "sample", "value")
names(gex.melt) <- c("gene", "sample", "value")


acc.melt$group <- c2$group[match(acc.melt$gene, c2$gene.x)]
gex.melt$group <- c2$group[match(gex.melt$gene, c2$gene.x)]

acc.melt$celltype <- substr(acc.melt$sample, 1, regexpr("_", acc.melt$sample) - 1)
gex.melt$celltype <- substr(gex.melt$sample, 1, regexpr("_", gex.melt$sample) - 1)

acc.melt$Phenotype <- ifelse(grepl("PD", acc.melt$sample), "PD", "Ctrl")
gex.melt$Phenotype <- ifelse(grepl("PD", gex.melt$sample), "PD", "Ctrl")


ggplot(acc.melt, aes(x = sample, y = gene, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "#290230", mid = "#CC4678FF", high = "#F0F921FF", midpoint = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(size = 6)) +
  xlab("") +
  ylab("") +
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1.5,
                               title.position = "top", title.hjust = 0.5)) +
  facet_grid(group ~ ., scales = "free", space = "free")

ggplot(gex.melt, aes(x = sample, y = gene, fill = value)) +
  geom_tile() +
  scale_fill_viridis(option = "viridis", begin = 0.05, end = 0.95) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(size = 6)) +
  xlab("") +
  ylab("") +
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1.5,
                               title.position = "top", title.hjust = 0.5)) +
  facet_grid(group ~ ., scales = "free", space = "free")

c2.2<-all


a10.gr<-c2.2[which(c2.2$group !="Ctrl"),]
c10.gr<-c2.2[which(c2.2$group !="ROT"),]
a10.gr$CT <- gsub("_ROT|_control", "", a10.gr$CT)
c10.gr$CT <- gsub("_ROT|_control", "", c10.gr$CT)
groups<-strsplit(a10.gr$CT, split=",")
t<-lapply(groups, function(i) i[order(i)])
t<-lapply(t, paste, collapse="&")
t<-as.character(t)
a10.gr$groups<-t


groups<-strsplit(c10.gr$CT, split=",")
t<-lapply(groups, function(i) i[order(i)])
t<-lapply(t, paste, collapse="&")
t<-as.character(t)
c10.gr$groups<-t



tab<-table(a10.gr$groups)
tab<-head(tab[order(-tab)], n=15)
tab2<-table(c10.gr$groups)
tab2<-head(tab2[order(-tab2)], n=15)

tab<-tab[names(tab) %in% names(tab2)]
tab2<-tab2[names(tab2) %in% names(tab)]

tab2<-tab2[names(tab)] ## do same order

######
install.packages("UpSetR")
library(UpSetR)


u<-upset(fromExpression(tab), nintersects=12, nsets=8,   sets.bar.color=c(
  "cornflowerblue","mediumorchid3","seagreen3","coral3","darkgoldenrod1"), main.bar.color="red3", text.scale=1.3  , mainbar.y.label="ROT Links",
  mainbar.y.max=10000, show.numbers=F, keep.order=T, order.by=c("freq","degree"))

u<-upset(fromExpression(tab), keep.order=T, order.by=c("freq","degree"), sets.bar.color=c(
  "cornflowerblue","mediumorchid3","seagreen3","coral3","darkgoldenrod1"),
  main.bar.color="red3", text.scale=1.3  , mainbar.y.label="ROT Links", show.numbers=F)


pdf("ROT_upset.pdf", width=8, height=5)
u
dev.off()


u2<-upset(fromExpression(tab2), nintersects=12, nsets=7,  sets.bar.color=c(
  "cornflowerblue","seagreen3","coral3","darkgoldenrod1","mediumorchid3","firebrick", "pink"),
  main.bar.color="blue", text.scale=1.3,  mainbar.y.label="Ctrl Links",
  mainbar.y.max=10000, show.numbers=F, keep.order=T, order.by=c("freq","degree"))
u2<-upset(fromExpression(tab2), keep.order=T, order.by=c("freq","degree"), sets.bar.color=c(
  "cornflowerblue","mediumorchid3","seagreen3","coral3","darkgoldenrod1"),
  main.bar.color="blue", text.scale=1.3  , mainbar.y.label="Ctrl Links")


pdf("Ctrl_upset.pdf", width=8, height=5)
u2
dev.off()



