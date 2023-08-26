#.libPaths("/Volumes/Seagate/multiome_cp/lib/R/library/")
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
data<-Read10X(data.dir = "/Volumes/Seagate/multiome_cp/sample_aggr/filtered_feature_bc_matrix/")
rna_counts <- data$`Gene Expression`
atac_counts <- data$Peaks
data <- CreateSeuratObject(counts = rna_counts,project="all")
numCells<-nrow(data@meta.data)

###################################################################################
#     Add atac assay
#   and get TSS, nucleosome signal, blacklist fraction
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Rnorvegicus.v105)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "rn7"

frag.file <- paste0("/Volumes/Seagate/multiome_cp/sample_aggr/atac_fragments.tsv.gz")
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^Mt-")
chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = 'rn7',
  fragments = frag.file,
  min.cells = 10,
  annotation = annotations
)
data[["ATAC"]] <- chrom_assay
DefaultAssay(data)<-"ATAC"
data <- NucleosomeSignal(data)
data <- TSSEnrichment(data)
data@meta.data
####新增phenotype欄位#####
# 引用這些包
library(dplyr)
library(stringr)
orig.ident <-rownames(data@meta.data)
# 使用stringr包來提取名字最後的數字
last.digit <- as.numeric(str_extract(orig.ident, "-\\d+$"))
# 根據這個數字添加phenotype欄位
data$phenotype <- ifelse(last.digit %in% c(-1, -2), "control", "ROT")
#######################
###################################################################################
#                     Filtering, norm, scale, dim reduc
VlnPlot(data, features =  c("nCount_ATAC", "nFeature_RNA","nCount_RNA","percent.mt"
                          ,'nucleosome_signal'),
        ncol = 3, pt.size = 0)
data<-subset(data, subset = nFeature_RNA < 10000 &   
               nFeature_RNA > 200 &
               nucleosome_signal < 2 &
               percent.mt < 5)
pdf("vinplot.pdf", width=14, height=8)
VlnPlot(data, features =  c("nCount_ATAC", "nFeature_RNA","nCount_RNA","percent.mt"
                            ,'nucleosome_signal'),
        ncol = 3, pt.size = 0)
dev.off()
data<-SCTransform(data, vars.to.regress = c("percent.mt"),verbose=F) %>% RunPCA(ndims=30) %>% FindNeighbors(dims = 1:30) %>% 
  RunUMAP(dims = 1:30, reduction.name="umap.rna", reduction.key = "rnaUMAP_") 
RunUMAP(data, dims = 1:30, reduction.name="umap.rna", reduction.key = "rnaUMAP_")
DefaultAssay(data)<-"ATAC"
data<-RunTFIDF(data) %>% FindTopFeatures(min.cutoff='q0')%>% RunSVD()%>% 
  RunUMAP(reduction='lsi',dims=2:50,reduction.name="umap.atac",reduction.key="atacUMAP_") 
step<-c(step,"filter")
numCells<-c(numCells, nrow(data@meta.data))
pdf("Init_Umap.pdf", width=14, height=8)
p1<-DimPlot(data, reduction="umap.rna")
p2<-DimPlot(data, reduction="umap.atac")
p1+p2
dev.off()

#######doublet filter
DefaultAssay(data)<-"SCT"
sce<-as.SingleCellExperiment(data)
dbl<-computeDoubletDensity(sce)
data$isDbl<-ifelse(dbl<3.5,1,0)
data$dbl<-dbl
pdf("Doublet.pdf", width=14, height=8)
p1<-DimPlot( data,group.by="isDbl",reduction="umap.rna")
p2<-FeaturePlot(data, features="dbl",reduction="umap.rna", max.cutoff=10)
p1+p2
dev.off()
data<-subset(data, subset= isDbl==1)
########
data <- SCTransform(data,assay="SCT",vars.to.regress = c("percent.mt"),verbose=F) %>% RunPCA(ndims=30) %>% 
  FindNeighbors(dims = 1:30) %>% RunUMAP(dims = 1:30, reduction.name="umap.rna",reduction.key = "rnaUMAP_")
data<-FindMultiModalNeighbors(data, reduction.list=list("pca","lsi"), dims.list=list(1:30,2:50))
data<-RunUMAP(data,nn.name="weighted.nn", reduction.name="wnn.umap",reduction.key="wnnUMAP_")
data<-FindClusters(data, graph.name="wsnn", algorithm=3)
###################################################################################
pdf("rna_atac.umap.pdf", width=14, height=6)
p1=DimPlot(data, reduction = "umap.rna") + ggtitle("RNA")+theme(legend.position="none")
p2=DimPlot(data, reduction = "umap.atac") + ggtitle("ATAC")+theme(legend.position="none")
p3=DimPlot(data, reduction="wnn.umap", label = TRUE)+ ggtitle("WNN")
p1+p2+p3
dev.off()
Idents(data)<-data$phenotype
DimPlot(data, reduction="wnn.umap", label = TRUE)
###################################################################################
######################Find Marker##############################
BiocManager::install("MAST")
DefaultAssay(data) <- "RNA"
Idents(data)<-data$seurat_clusters
all.markers<-FindAllMarkers(data, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5, test.use = "MAST")
dim(all.markers)
table(all.markers$cluster)
top3_markers <- as.data.frame(all.markers %>% group_by(cluster)%>%top_n(n=3, wt = avg_log2FC))
top3_markers
##############
markers_cluster0 <- FindConservedMarkers(data,
                                         ident.1 = 0,
                                         grouping.var = 'phenotype')
head(markers_cluster0)

markers_cluster1 <- FindConservedMarkers(data,
                                         ident.1 = 1,
                                         grouping.var = 'phenotype')
head(markers_cluster1)

markers_cluster4 <- FindConservedMarkers(data,
                                         ident.1 = 4,
                                         grouping.var = 'phenotype')
head(markers_cluster4)

markers_cluster6 <- FindConservedMarkers(data,
                                         ident.1 = 6,
                                         grouping.var = 'phenotype')
head(markers_cluster6)

markers_cluster17 <- FindConservedMarkers(data,
                                         ident.1 = 17,
                                         grouping.var = 'phenotype')
head(markers_cluster17)

markers_cluster14 <- FindConservedMarkers(data,
                                         ident.1 = 14,
                                         grouping.var = 'phenotype')
head(markers_cluster14)

markers_cluster15 <- FindConservedMarkers(data,
                                          ident.1 = 15,
                                          grouping.var = 'phenotype')
head(markers_cluster15)

markers_cluster3 <- FindConservedMarkers(data,
                                          ident.1 = 3,
                                          grouping.var = 'phenotype')
head(markers_cluster3)


####### naming the cluster######
Idents(data)<-data$seurat_clusters
data <- RenameIdents(data, c("0" = "Drd2-MSN", 
                             "1" = "Drd1-MSN", 
                             "2" = "Olig",
                             "3" = "Olig",
                             "4" = "Drd1-MSN",
                             "5" = "Astrocyte",
                             "6" = "Drd1-MSN",
                             "7" = "Olig",
                             "8" = "Drd1-MSN",
                             "9" = "Drd2-MSN",
                             "10" = "Drd1-MSN",
                             "11" = "Microglia", 
                             "12" = "OPC",
                             "13" = "Olig",
                             "14" = "Glutamatergic",
                             "15" = "Interneuron",
                             "16" = "Microglia",
                             "17" = "Astrocyte",
                             "18" = "Microglia",
                             "19" = "Astrocyte",
                             "20" = "Astrocyte",
                             "24" = "OPC",
                             "25" = "Astrocyte"), 
                     new.names = "seurat_clusters")
data2<-subset(x=data, idents = c("Drd1-MSN","Drd2-MSN","Olig",
                                 "OPC","Microglia","Astrocyte", 'Glutamatergic',
                                 "Interneuron"))
data2$predict.id<-Idents(object = data2)
pdf('Seurat_cluster.pdf',width=14, height=6)
P1<-DimPlot(data, repel = TRUE, reduction = "wnn.umap", group.by = "seurat_clusters")
P2<-DimPlot(data2, repel = TRUE, reduction = "wnn.umap", group.by = "ident")
P1+P2
dev.off()

######################Peak-calling##############################
DefaultAssay(data2)<-"ATAC"
CTpeaks<-CallPeaks(data2, macs2.path="/opt/anaconda3/envs/env01/bin/macs2", group.by="predict.id")

chr<-c()
for (i in 1:22){
  chr<-c(chr, paste0("",i))
}
chr<-c(chr,"chrX")
chr<-c(chr,"chrY")
CTpeaks2<-CTpeaks[seqnames(CTpeaks) %in% chr]


df <- as.data.frame(CTpeaks)
df<-df[which(df$seqnames %in% chr),]
#bed file for reanalyze [sort -k1,1 -k2,2n CT_peak_set.bed]
write.table(df[,1:3], file="mutiome_CT_peak_set.bed", quote=F, sep="\t", row.names=F, col.names=F)
CTpeaks2<-makeGRangesFromDataFrame(df)
macs2_counts<-FeatureMatrix(fragments=Fragments(data2), features=CTpeaks2,cells=colnames(data2))

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Rnorvegicus.v105)
seqlevelsStyle(annotation) <- "UCSC"
data2[["CTpeaks"]]<-CreateChromatinAssay(counts=macs2_counts,  fragments=Fragments(data2))
#get peak list labeled by celltype
write.csv(df,"~/mutiome_CTpeaks_annotated.csv")



saveRDS(data2, "mutiome_final.rds")

# get AD and Ctrl barcodes
df$status<-data2$Status
ROT<-df[which(df$status=="ROT"),1:2]
Ctrl<-df[which(df$status=="Control"),1:2]
write.csv(ROT, "/barcodes_ROT.csv", row.names=F, quote=F)
write.csv(Ctrl, "~/barcodes_Controll.csv", row.names=F, quote=F)
###########
data<-data2
data<-readRDS("mutiome_final.rds")
all.genes <- rownames(data)
data <- ScaleData(data, features = all.genes)



Idents(data)<-'predict.id'
DefaultAssay(data) <- "CTpeaks"

# first compute the GC content for each peak
data<-RegionStats(data,genome=genome_rn7)

# link peaks to genes
data <- LinkPeaks(
  object = data,
  peak.assay = "CTpeaks",
  expression.assay = "SCT",
  genes.use = c("Aqp4")
)
data[['CTpeaks']]

idents.plot <- c("Drd2-MSN", "Drd1-MSN", "Olig",
                 "Astrocyte", "Microglia", "OPC", "Glutamatergic",'Interneuron')

p1 <- CoveragePlot(
  object = data,
  region = "Drd2",
  
  expression.assay = "SCT",
  idents = idents.plot,
  extend.upstream = 500,
  extend.downstream = 10000
)
p2 <- CoveragePlot(
  object = data,
  region = "Drd1",
  expression.assay = "SCT",
  idents = idents.plot,
  extend.upstream = 500,
  extend.downstream = 10000
)
p3 <- CoveragePlot(
  object = data,
  region = "Ppp1r1b",
  expression.assay = "SCT",
  idents = idents.plot,
  extend.upstream = 500,
  extend.downstream = 10000
)
p4 <- CoveragePlot(
  object = data,
  region = "Rgs9",
  expression.assay = "SCT",
  idents = idents.plot,
  extend.upstream = 500,
  extend.downstream = 10000
)
p5 <- CoveragePlot(
  object = data,
  region = "Mbp",
  expression.assay = "SCT",
  idents = idents.plot,
  extend.upstream = 500,
  extend.downstream = 10000
)
p6 <- CoveragePlot(
  object = data,
  region = "Mog",
  expression.assay = "SCT",
  idents = idents.plot,
  extend.upstream = 500,
  extend.downstream = 10000
)
p7 <- CoveragePlot(
  object = data,
  region = "S100b",
  expression.assay = "SCT",
  idents = idents.plot,
  extend.upstream = 500,
  extend.downstream = 10000
)
p8 <- CoveragePlot(
  object = data,
  region = "Aqp4",
  expression.assay = "SCT",
  idents = idents.plot,
  extend.upstream = 500,
  extend.downstream = 10000
)
p9 <- CoveragePlot(
  object = data,
  region = "Cx3cr1",
  expression.assay = "SCT",
  idents = idents.plot,
  extend.upstream = 500,
  extend.downstream = 10000
)
p10 <- CoveragePlot(
  object = data,
  region = "Grin2a",
  expression.assay = "SCT",
  idents = idents.plot,
  extend.upstream = 500,
  extend.downstream = 10000
)
pdf('coverageplot_1.pdf', width=20, height=9)
patchwork::wrap_plots(p1, p2,p4,p5, p6, p7, p8, p9, ncol = 4)
dev.off()
marker_gene<-c('Drd2','Drd1','Ppp1r1b', 'Rgs9','Mbp','Mog', 'S100b', 'Aqp4', 'Cx3cr1', 'Grin2a')






#########################
DefaultAssay(data)<-"RNA"
data<-NormalizeData(data, assay="RNA") %>% ScaleData(features=rownames(data))

avgExp_PC<-as.data.frame(AverageExpression(data, features=rownames(data), group.by=c("predict.id","phenotype"), slot="data", assay="RNA"))


colnames(avgExp_PC)<-gsub("RNA.","",colnames(avgExp_PC))

df<-data.frame(celltype=sapply(strsplit(colnames(avgExp_PC), "_",fixed=T), `[`,1), 
               id=sapply(strsplit(colnames(avgExp_PC), "_",fixed=T), `[`,2))
df$tmp<-substr(df$celltype,start=1,stop=4)

labels<-rownames(data)
labels$tmp<-substr(rownames(data),start=1,stop=4)

# Normalize
mat_scaled<-t(apply(avgExp_PC,1,scale))
mat_scaled<-t(scale(t(avgExp_PC)))

library(ggplot2)
library(RColorBrewer)
library(viridis)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col=sample(col_vector, 15)
names(col)<-unique(df$id)
ha<-rowAnnotation(celltype=labels$tmp,col= list(celltype=c("Astr"="darkgoldenrod1","Micr"="cornflowerblue","OPC"="seagreen3","Olig"="mediumorchid3","Drd1"="coral3","Drd2"="firebrick","Glut"="orange","Inte"="green" )), show_legend=F)
ID_colors <- c("Level1" = "red", "Level2" = "green", "Level3" = "blue")
ha2 <- HeatmapAnnotation(
  celltype = df$tmp,
  ID = df$id,
  col = list(
    celltype = c(
      "Astr" = "darkgoldenrod1",
      "Micr" = "cornflowerblue",
      "OPC" = "seagreen3",
      "Olig" = "mediumorchid3",
      "Drd1" = "coral3",
      "Drd2" = "firebrick",
      "Glut" = "orange",
      "Inte" = "green"
    ),
    ID = ID_colors
  ),
  show_legend = F
)
ht=Heatmap(mat_scaled,cluster_rows=F,cluster_columns=F,col=colorRamp2(c(-1,0,2),viridis(3)), row_split=labels$tmp, top_annotation=ha2, name="Z", left_annotation=ha, show_column_names=F, show_row_names=F, column_title=NULL, row_title_gp=gpar(fontsize=14))

#lab<-which(rownames(mat_scaled) %in% c("AQP4","NRGN","GAD2","P2RY12","MOBP","VCAN"))

pdf("scREAD_heatmap_byID.pdf", width=5,height=8)
ht
dev.off()
####DO DeG ROT/Control##########
DefaultAssay(data)<-"RNA"
all.genes<-rownames(data)
data<-NormalizeData(data) 

Drd1_neuron<-subset(data, subset=predict.id=="Drd1-MSN")
Drd2_neuron<-subset(data, subset=predict.id=="Drd2-MSN")
Olig<-subset(data, subset=predict.id=="Olig")
Astrocyte<-subset(data, subset=predict.id=="Astrocyte")
Microglia<-subset(data, subset=predict.id=="Microglia")
OPC<-subset(data, subset=predict.id=="OPC")
Glutamatergic<-subset(data, subset=predict.id=="Glutamatergic")
Interneuron<-subset(data, subset=predict.id=="Interneuron")


subs<-list(Drd1_neuron,Drd2_neuron,Olig,Astrocyte,Microglia,OPC, Glutamatergic,Interneuron )

library(MAST)
Idents(data)<-data$phenotype
i=1
Idents(subs[[i]])<-subs[[i]]$phenotype
subs[[i]]<-NormalizeData(subs[[i]])
ALL<-FindMarkers(subs[[i]], ident.1="ROT",ident.2="control",min.pct=0.25,test.use="MAST", assay="RNA") 
ALL$celltype<-subs[[i]]$predict.id[1]
ALL$gene<-rownames(ALL)
for (i in 2:length(subs)){
  Idents(subs[[i]])<-subs[[i]]$phenotype
  subs[[i]]<-NormalizeData(subs[[i]])%>% ScaleData(features=row.names(data))
  markers<-FindMarkers(subs[[i]], ident.1="ROT",ident.2="control",min.pct=0.25, test.use="MAST",
                       assay="RNA") 
  markers$celltype<-subs[[i]]$predict.id[1]
  markers$gene<-rownames(markers)
  ALL<-rbind(ALL,markers)
}


ALL$cat<-ifelse(ALL$avg_log2FC >0, "up","down")

significant<-ALL[which(abs(ALL$avg_log2FC)>0.25),]
significant<-significant[which(significant$p_val_adj<0.01),]

write.csv(significant, "mutiome_DEGs_MAST_ROTCtrl.csv")
significant<-read.csv("mutiome_DEGs_MAST_ROTCtrl.csv")

degs<-significant
degs<-degs[order(degs$celltype),]
degs_t<-degs[which(degs$cat=="up"),]
du<-degs_t[,c(8,7)]
gene_ct<-table(du$gene, du$celltype)

mat<-matrix(ncol=8,nrow=8)
colnames(mat)<-unique(degs$celltype)
rownames(mat)<-unique(degs$celltype)
for (j in 1:8){
  for (i in 1:8){
    if (i>=j){
      mat[j,i]=nrow(gene_ct[which(gene_ct[,i]==1 & gene_ct[,j]==1),])
    }
    else{
      mat[j,i]<-0
    }
  }
}




degs_t<-degs[which(degs$cat=="down"),]
du<-degs_t[,c(8,7)]
gene_ct<-table(du$gene, du$celltype)

mat2<-matrix(ncol=8,nrow=8)
colnames(mat2)<-unique(degs$celltype)
rownames(mat2)<-unique(degs$celltype)
for (j in 1:8){
  for (i in 1:8){
    if (i<=j){
      mat2[j,i]=nrow(gene_ct[which(gene_ct[,i]==1 & gene_ct[,j]==1),])
    }
    else{
      mat2[j,i]<-0
    }
  }
}
mat2[is.na(mat2)] <- 0

colnames(mat)<-substr(colnames(mat),start=1,stop=4)
colnames(mat2)<-colnames(mat)
rownames(mat)<-colnames(mat)
rownames(mat2)<-colnames(mat)

ha<-rowAnnotation(celltype=colnames(mat)
                  , col= list(celltype=c("Astr"="darkgoldenrod1","Drd1"="cornflowerblue","Drd2"="seagreen3","Glut"="mediumorchid3","Inte"="coral3","Micr"="firebrick", "Olig"="orange", "OPC"="green")), show_legend=F,annotation_label="")
ha2<-HeatmapAnnotation(celltype=substr(colnames(mat),start=1,stop=4)
                       , col= list(celltype=c("Astr"="darkgoldenrod1","Drd1"="cornflowerblue","Drd2"="seagreen3","Glut"="mediumorchid3","Inte"="coral3","Micr"="firebrick", "Olig"="orange", "OPC"="green")), show_legend=F,annotation_label="")

ht=Heatmap(mat, cluster_rows=F,cluster_columns=F,col=colorRamp2(c(0,20,22,26),c("white","red","red","grey80")),  top_annotation=ha2, name="Up", left_annotation=ha, show_column_names=T, show_row_names=T, column_title=NULL, row_title_gp=gpar(fontsize=14),cell_fun = function(j, i, x, y, w, h, fill) {
  if(mat[i, j] >0) {
    if(i==j){
      grid.text(mat[i,j], x, y, gp=gpar(fontsize=10, col="red"),  vjust=1)
    }
    else grid.text(mat[i,j], x, y, gp=gpar(fontsize=10))
  }
})
ht2=Heatmap(mat2, cluster_rows=F,cluster_columns=F,col=colorRamp2(c(0,34,35,36),c("white","dodgerblue","dodgerblue","grey80")),  top_annotation=ha2, name="Down", left_annotation=ha, show_column_names=T, show_row_names=T, column_title=NULL, row_title_gp=gpar(fontsize=14),cell_fun = function(j, i, x, y, w, h, fill) {
  if(mat2[i, j] >0) {
    if(i==j){
      grid.text(mat2[i,j], x, y, gp=gpar(fontsize=10, col="blue"),  vjust=1)
    }
    else{
      grid.text(mat2[i,j], x, y, gp=gpar(fontsize=10))
    }
  } })

pdf("Shared_degs_heatmap_up.pdf", width=8,height=4)
ht+ht2
dev.off()

###################################
###DO DAE ROT/Control##########
DefaultAssay(data)<-"ATAC"
all.genes<-rownames(data)
data<-NormalizeData(data) 

Drd1_neuron<-subset(data, subset=predict.id=="Drd1-MSN")
Drd2_neuron<-subset(data, subset=predict.id=="Drd2-MSN")
Olig<-subset(data, subset=predict.id=="Olig")
Astrocyte<-subset(data, subset=predict.id=="Astrocyte")
Microglia<-subset(data, subset=predict.id=="Microglia")
OPC<-subset(data, subset=predict.id=="OPC")
Glutamatergic<-subset(data, subset=predict.id=="Glutamatergic")
Interneuron<-subset(data, subset=predict.id=="Interneuron")


subs<-list(Drd1_neuron,Drd2_neuron,Olig,Astrocyte,Microglia,OPC, Glutamatergic,Interneuron )
library(MAST)
Idents(data)<-data$phenotype
i=1
Idents(subs[[i]])<-subs[[i]]$phenotype
subs[[i]]<-NormalizeData(subs[[i]])
ALL<-FindMarkers(subs[[i]], ident.1="control",ident.2="ROT",min.pct=0.25, assay="ATAC") 
ALL$celltype<-subs[[i]]$predict.id[1]
ALL$gene<-rownames(ALL)
for (i in 2:length(subs)){
  Idents(subs[[i]])<-subs[[i]]$phenotype
  subs[[i]]<-NormalizeData(subs[[i]])%>% ScaleData(features=row.names(data))
  markers<-FindMarkers(subs[[i]], ident.1="control",ident.2="ROT",min.pct=0.25,test.use = "MAST",
                       assay="ATAC") 
  markers$celltype<-subs[[i]]$predict.id[1]
  markers$gene<-rownames(markers)
  ALL<-rbind(ALL,markers)
  gc()
}
#subs[[1]]
#Idents(subs[[1]])<-subs[[1]]$phenotype
#Idents(subs[[1]])
#subs[[1]]<-NormalizeData(subs[[1]])%>% ScaleData(features=row.names(data))

ALL$cat<-ifelse(ALL$avg_log2FC >0, "up","down")

significant<-ALL[which(abs(ALL$avg_log2FC)>0.25),]
significant<-significant[which(significant$p_val_adj<0.05),]

write.csv(significant, "mutiome_DAEs_MAST_ROTCtrl.csv")
significant<-read.csv("mutiome_DAEs_MAST_ROTCtrl.csv")

degs<-significant
degs<-degs[order(degs$celltype),]
degs_t<-degs[which(degs$cat=="up"),]
du<-degs_t[,c(6,7)]
gene_ct<-table(du$gene, du$celltype)
#degs
tab<-table(degs$celltype, degs$cat)

mat<-matrix(nrow=length(unique(degs$gene)),ncol=5)
rownames(mat)<-unique(degs$gene)
colnames(mat)<-unique(degs$celltype)

for (i in 1:nrow(mat)){
  for (j in 1:ncol(mat)){
    gene_tmp<-degs[which(degs$gene==rownames(mat)[i]),]
    gene_tmp<-gene_tmp[which(gene_tmp$celltype==colnames(mat)[j]),]
    mat[i,j]<-ifelse(nrow(gene_tmp)>0, gene_tmp$avg_log2FC,0)
  }
}

colnames(mat)<-substr(colnames(mat),start=1,stop=4)
colnames(mat)[colnames(mat) == "Inte"] <- "Micro"

library(ComplexHeatmap)
library(colorRamp2)
# # all labels
noDup<-degs[!duplicated(degs$gene),]
mat<-mat[order(noDup$celltype, noDup$avg_log2FC),]
mat<-mat[,c(5,1,2,3, 8)]
ha2<-HeatmapAnnotation(celltype=colnames(mat), col= list(celltype=c("Astr"="darkgoldenrod1","Olig"="mediumorchid3","Drd1"="coral3","Drd2"="firebrick","Micro"="green" )), show_legend=F)
bar1<-HeatmapAnnotation(up=anno_barplot(tab[,2], gp=gpar(fill="red"), axis_param=list(labels=c("","","","")), ylim=c(0,280)),
                        down=anno_barplot((tab[,1] * -1), gp=gpar(fill="blue"), axis_param=list(labels=c("","","")), ylim=c(-280,0),
                                          show_annotation_name=F,
                                          show_legend=F))
ha<-c( bar1, ha2)
# # 
# 
# # 
noDup<-noDup[order(noDup$celltype, noDup$avg_log2FC),]
#lab<-which(noDup$gene %in% c(ad_related$Gene,"PLCG2","MAPT","ARL17B","SAMD4A","PTPRG","MDGA2","GPR158") & abs(noDup$avg_log2FC)>0.5)
ht=Heatmap(mat, cluster_rows=F,show_row_dend=F,cluster_columns=F,col=colorRamp2(c(-1,0,1),c("blue","white","orange")), name="MAST-adj log2FC",  show_column_names=T, show_row_names=F, column_title=NULL,row_names_gp = gpar(fontsize = 3), top_annotation=ha)
ht=Heatmap(mat, cluster_rows=F,show_row_dend=F,cluster_columns=F,col=colorRamp2(c(-1,0,1),c("blue","white","orange")), name="MAST-adj log2FC",  show_column_names=T, show_row_names=F, column_title=NULL,row_names_gp = gpar(fontsize = 3))


pdf("ROT_DAGs_MAST_heatmap_foldCs.pdf", width=6, height=8)
ht
dev.off()


#########Gene set enrichment###########
library(enrichR)
degs<-significant
setEnrichrSite("Enrichr")

dbs <- c("GO_Molecular_Function_2021", "GO_Cellular_Component_2021", "GO_Biological_Process_2021", "KEGG_2021_Human","Azimuth_Cell_Types_2021")


ast_down<-enrichr(degs[which(degs$celltype=="Astrocyte" & degs$cat=="down"),]$gene, dbs)
ast_up<-enrichr(degs[which(degs$celltype=="Astrocyte" & degs$cat=="up"),]$gene, dbs)


mic_down<-enrichr(degs[which(degs$celltype=="Microglia" & degs$cat=="down"),]$gene, dbs)
mic_up<-enrichr(degs[which(degs$celltype=="Microglia" & degs$cat=="up"),]$gene, dbs)


OPC_down<-enrichr(degs[which(degs$celltype=="OPC" & degs$cat=="down"),]$gene, dbs)
OPC_up<-enrichr(degs[which(degs$celltype=="OPC" & degs$cat=="up"),]$gene, dbs)


Olig_down<-enrichr(degs[which(degs$celltype=="Olig" & degs$cat=="down"),]$gene, dbs)
Olig_up<-enrichr(degs[which(degs$celltype=="Olig" & degs$cat=="up"),]$gene, dbs)


Drd1_down<-enrichr(degs[which(degs$celltype=="Drd1-MSN" & degs$cat=="down"),]$gene, dbs)
Drd1_up<-enrichr(degs[which(degs$celltype=="Drd1-MSN" & degs$cat=="up"),]$gene, dbs)


Drd2_down<-enrichr(degs[which(degs$celltype=="Drd2-MSN" & degs$cat=="down"),]$gene, dbs)
Drd2_up<-enrichr(degs[which(degs$celltype=="Drd2-MSN" & degs$cat=="up"),]$gene, dbs)

Glutamatergic_down<-enrichr(degs[which(degs$celltype=="Glutamatergic" & degs$cat=="down"),]$gene, dbs)
Glutamatergic_up<-enrichr(degs[which(degs$celltype=="Glutamatergic" & degs$cat=="up"),]$gene, dbs)

Interneuron_down<-enrichr(degs[which(degs$celltype=="Interneuron" & degs$cat=="down"),]$gene, dbs)
Interneuron_up<-enrichr(degs[which(degs$celltype=="Interneuron" & degs$cat=="up"),]$gene, dbs)


#### down


a<-ast_down[["GO_Biological_Process_2021"]]
b<-mic_down[["GO_Biological_Process_2021"]]
c<-OPC_down[["GO_Biological_Process_2021"]]
d<-Olig_down[["GO_Biological_Process_2021"]]
e<-Drd1_down[["GO_Biological_Process_2021"]]
f<-Drd2_down[["GO_Biological_Process_2021"]]
g<-Glutamatergic_down[["GO_Biological_Process_2021"]]
h<-Interneuron_down[["GO_Biological_Process_2021"]]


a$celltype<-"Astrocyte"
b$celltype<-"Microglia"
c$celltype<-"OPC"
d$celltype<-"Olig"
e$celltype<-"Drd1-MSN"
f$celltype<-"Drd2-MSN"
g$celltype<-"Glutamatergic"
h$celltype<-"Interneuron"


df<-rbind(a,b,c,d,e,f, g,h)
df2<-df[which(df$Adjusted.P.value<0.05),]
df2<-df2[order(df2$Adjusted.P.value),]
top<-df2[c(1:10),1]
df2_BP_down<-df2

df2<-df[which(df$Adjusted.P.value<0.01),]

# write.csv(df2, "~/scMultiomics_AD/enrichr/top_BP_downAD_across_ct_MAST.csv")
write.csv(df2, "top_BP_downROT_across_ct_MAST.csv")

t<-table(df2_BP_down$Term)
topDown<-head(df2$Term, n=10)
topDown<-df2 %>% group_by(celltype) %>% top_n(n=2, wt=Combined.Score)
topDown<-topDown$Term
tD<-df[which(df$Term %in% topDown),]

mat<-matrix(nrow=6, ncol=8)
colnames(mat)<-unique(df$celltype)
for (j in 1:6){
  for (i in 1:8){
    termTmp<-unique(topDown)[j]
    ctTmp<-colnames(mat)[i]
    sub<-df[which(df$Term==termTmp & df$celltype==ctTmp),]
    if (nrow(sub)>0){
      mat[j,i]<-sub$Odds.Ratio
    }
    else{
      mat[j,i]<-0
    }
  }  
}

rownames(mat)<-sapply(strsplit(unique(topDown), " (", fixed=T), `[`, 1)
colnames(mat)<-substr(colnames(mat),start=1,stop=4)

sig_mat<-matrix(nrow=6, ncol=8)
colnames(sig_mat)<-unique(df$celltype)
for (j in 1:6){
  for (i in 1:8){
    termTmp<-unique(topDown)[j]
    ctTmp<-colnames(sig_mat)[i]
    sub<-df[which(df$Term==termTmp & df$celltype==ctTmp),]
    if (nrow(sub)>0){
      sig_mat[j,i]<-sub$Adjusted.P.value}
    else{
      sig_mat[j,i]<-1}  }}


#mat<-mat[,c(1,3,4,2,5,6)]
#sig_mat<-sig_mat[,c(1,3,4,2,5,6)]
#BiocManager::install("ComplexHeatmap")
#install.packages("colorRamp2")
library(colorRamp2)
library(ComplexHeatmap)
ha2<-HeatmapAnnotation(celltype=colnames(mat)
                       , col= list(celltype=c("Astr"="darkgoldenrod1","Micr"="cornflowerblue","OPC"="seagreen3","Olig"="mediumorchid3","Drd1"="coral3","Drd2"="firebrick","Glut"="orange","Inte"="green" )), show_legend=F,annotation_label="")

ht=Heatmap(mat, cluster_rows=F,cluster_columns=F,col=colorRamp2(c(0,35,100),c("grey95","deepskyblue","dodgerblue4")),  top_annotation=ha2, name="Odds.Ratio",  show_column_names=T, show_row_names=T, column_title=NULL,row_names_side="left", row_title_gp=gpar(fontsize=14),row_names_max_width = max_text_width(
  rownames(mat), 
  gp = gpar(fontsize = 12)
), cell_fun = function(j, i, x, y, w, h, fill) {
  if(sig_mat[i, j] <0.05) {
    grid.text("*", x, y, gp=gpar(fontsize=20, col="white"), vjust="center")
  } })

pdf("Top_downDEG_GObp.pdf", width=13, height=5)
ht
dev.off()

#            BP UP
a<-ast_up[["GO_Biological_Process_2021"]]
b<-mic_up[["GO_Biological_Process_2021"]]
c<-OPC_up[["GO_Biological_Process_2021"]]
d<-Olig_up[["GO_Biological_Process_2021"]]
e<-Drd1_up[["GO_Biological_Process_2021"]]
f<-Drd2_up[["GO_Biological_Process_2021"]]
g<-Glutamatergic_up[["GO_Biological_Process_2021"]]
h<-Interneuron_up[["GO_Biological_Process_2021"]]

a$celltype<-"Astrocyte"
b$celltype<-"Microglia"
c$celltype<-"OPC"
d$celltype<-"Olig"
e$celltype<-"Drd1-MSN"
f$celltype<-"Drd2-MSN"
g$celltype<-"Glutamatergic"
h$celltype<-"Interneuron"

df<-rbind(a,b,c,d,e,f,g,h)
df_all<-df
df2<-df[which(df$Adjusted.P.value<0.05),]
df2<-df2[order(df2$Adjusted.P.value),]
top<-df2[c(1:10),1]


write.csv(df2, "top_BP_upAD_across_ct_MAST.csv")
############
# up
t<-table(df2$Term)
topup<-head(df2$Term, n=10)
topup<-df2 %>% group_by(celltype) %>% top_n(n=5, wt=Combined.Score)
topup<-topup[order(topup$celltype),]
topup<-topup$Term
tD<-df[which(df$Term %in% topup),]

mat<-matrix(nrow=length(unique(topup)), ncol=8)
colnames(mat)<-unique(df$celltype)
for (j in 1:length(unique(topup))){
  for (i in 1:8){
    termTmp<-unique(topup)[j]
    ctTmp<-colnames(mat)[i]
    sub<-df[which(df$Term==termTmp & df$celltype==ctTmp),]
    if (nrow(sub)>0){
      mat[j,i]<- sub$Odds.Ratio
    }
    else{
      mat[j,i]<-0
    }
  }  
}

rownames(mat)<-sapply(strsplit(unique(topup), " (", fixed=T), `[`, 1)
colnames(mat)<-substr(colnames(mat),start=1,stop=4)

sig_mat<-matrix(nrow=length(unique(topup)), ncol=8)
colnames(sig_mat)<-unique(df$celltype)
for (j in 1:length(unique(topup))){
  for (i in 1:8){
    termTmp<-unique(topup)[j]
    ctTmp<-colnames(sig_mat)[i]
    sub<-df[which(df$Term==termTmp & df$celltype==ctTmp),]
    if (nrow(sub)>0){
      sig_mat[j,i]<-sub$Adjusted.P.value}
    else{
      sig_mat[j,i]<-1}  }}

#mat<-mat[,c(1,3,4,2,5,6)]
#sig_mat<-sig_mat[,c(1,3,4,2,5,6)]

ha2<-HeatmapAnnotation(celltype=colnames(mat)
                       , col= list(celltype=c("Astr"="darkgoldenrod1","Micr"="cornflowerblue","OPC"="seagreen3","Olig"="mediumorchid3","Drd1"="coral3","Drd2"="firebrick","Glut"="orange","Inte"="green" )), show_legend=F,annotation_label="")

ht=Heatmap(mat, cluster_rows=F,cluster_columns=F,col=colorRamp2(c(0,35,100),c("grey95","red","red4")),  top_annotation=ha2, name="Odds.Ratio",  show_column_names=T, show_row_names=T, column_title=NULL,row_names_side="left", row_title_gp=gpar(fontsize=14),row_names_max_width = max_text_width(
  rownames(mat), 
  gp = gpar(fontsize = 12)
), cell_fun = function(j, i, x, y, w, h, fill) {
  if(sig_mat[i, j] <0.05) {
    grid.text("*", x, y, gp=gpar(fontsize=20, col="white"), vjust="center")
  } })

pdf("Top_upDEG_GObp.pdf", width=10, height=5)
ht
dev.off()

##################
#degs
tab<-table(degs$celltype, degs$cat)

mat<-matrix(nrow=length(unique(degs$gene)),ncol=8)
rownames(mat)<-unique(degs$gene)
colnames(mat)<-unique(degs$celltype)

for (i in 1:nrow(mat)){
  for (j in 1:ncol(mat)){
    gene_tmp<-degs[which(degs$gene==rownames(mat)[i]),]
    gene_tmp<-gene_tmp[which(gene_tmp$celltype==colnames(mat)[j]),]
    mat[i,j]<-ifelse(nrow(gene_tmp)>0, gene_tmp$avg_log2FC,0)
  }
}

colnames(mat)<-substr(colnames(mat),start=1,stop=4)


# # all labels
noDup<-degs[!duplicated(degs$gene),]
mat<-mat[order(noDup$celltype, noDup$avg_log2FC),]
mat<-mat[,c(5,6,4,1,2,3,7, 8)]
ha2<-HeatmapAnnotation(celltype=colnames(mat), col= list(celltype=c("Astr"="darkgoldenrod1","Micr"="cornflowerblue","OPC"="seagreen3","Olig"="mediumorchid3","Drd1"="coral3","Drd2"="firebrick","Glut"="orange","Inte"="green" )), show_legend=F)
bar1<-HeatmapAnnotation(up=anno_barplot(tab[,2], gp=gpar(fill="red"), axis_param=list(labels=c("","","","")), ylim=c(0,280)),
                        down=anno_barplot((tab[,1] * -1), gp=gpar(fill="blue"), axis_param=list(labels=c("","","")), ylim=c(-280,0),
                                          show_annotation_name=F,
                                          show_legend=F))
ha<-c( bar1, ha2)
# # 

# 
# # 
noDup<-noDup[order(noDup$celltype, noDup$avg_log2FC),]
#lab<-which(noDup$gene %in% c(ad_related$Gene,"PLCG2","MAPT","ARL17B","SAMD4A","PTPRG","MDGA2","GPR158") & abs(noDup$avg_log2FC)>0.5)
ht=Heatmap(mat, cluster_rows=F,show_row_dend=F,cluster_columns=F,col=colorRamp2(c(-1,0,1),c("blue","white","red")), name="MAST-adj log2FC",  show_column_names=T, show_row_names=F, column_title=NULL,row_names_gp = gpar(fontsize = 3), top_annotation=ha)


pdf("ROT_DEGs_MAST_heatmap_foldCs.pdf", width=6, height=8)
ht
dev.off()

#############
degs<-degs[order(degs$celltype),]
degs_t<-degs[which(degs$cat=="up"),]
du<-degs_t[,c(7,6)]
gene_ct<-table(du$gene, du$celltype)

mat<-matrix(ncol=6,nrow=6)
colnames(mat)<-unique(degs$celltype)
rownames(mat)<-unique(degs$celltype)
for (j in 1:6){
  for (i in 1:6){
    if (i>=j){
      mat[j,i]=nrow(gene_ct[which(gene_ct[,i]==1 & gene_ct[,j]==1),])
    }
    else{
      mat[j,i]<-0
    }
  }
}




degs_t<-degs[which(degs$cat=="down"),]
du<-degs_t[,c(7,6)]
gene_ct<-table(du$gene, du$celltype)

mat2<-matrix(ncol=6,nrow=6)
colnames(mat2)<-unique(degs$celltype)
rownames(mat2)<-unique(degs$celltype)
for (j in 1:6){
  for (i in 1:6){
    if (i<=j){
      mat2[j,i] = sum(gene_ct[,i]==1 & gene_ct[,j]==1)
      
    }
    else{
      mat2[j,i]<-0
    }
  }
}


colnames(mat)<-substr(colnames(mat),start=1,stop=3)
colnames(mat2)<-colnames(mat)
rownames(mat)<-colnames(mat)
rownames(mat2)<-colnames(mat)

ha<-rowAnnotation(celltype=colnames(mat)
                  , col= list(celltype=c("DA-neuron"="darkgoldenrod1","GABA-neuron"="cornflowerblue","Olig"="seagreen3","Astrocyte"="mediumorchid3","Microglia"="coral3","OPC"="firebrick")), show_legend=F,annotation_label="")
ha2<-HeatmapAnnotation(celltype=substr(colnames(mat),start=1,stop=3)
                       , col= list(celltype=c("DA-neuron"="darkgoldenrod1","GABA-neuron"="cornflowerblue","Olig"="seagreen3","Astrocyte"="mediumorchid3","Microglia"="coral3","OPC"="firebrick")), show_legend=F,annotation_label="")

ht=Heatmap(mat, cluster_rows=F,cluster_columns=F,col=colorRamp2(c(0,20,25,26),c("white","red","red","grey80")),  top_annotation=ha2, name="Up", left_annotation=ha, show_column_names=T, show_row_names=T, column_title=NULL, row_title_gp=gpar(fontsize=14),cell_fun = function(j, i, x, y, w, h, fill) {
  if(mat[i, j] >0) {
    if(i==j){
      grid.text(mat[i,j], x, y, gp=gpar(fontsize=10, col="red"),  vjust=-1.5)
    }
    else grid.text(mat[i,j], x, y, gp=gpar(fontsize=10))
  }
})
ht2=Heatmap(mat2, cluster_rows=F,cluster_columns=F,col=colorRamp2(c(0,34,35,36),c("white","dodgerblue","dodgerblue","grey80")),  top_annotation=ha2, name="Down", left_annotation=ha, show_column_names=T, show_row_names=T, column_title=NULL, row_title_gp=gpar(fontsize=14),cell_fun = function(j, i, x, y, w, h, fill) {
  if(mat2[i, j] >0) {
    if(i==j){
      grid.text(mat2[i,j], x, y, gp=gpar(fontsize=10, col="blue"),  vjust=1.5)
    }
    else{
      grid.text(mat2[i,j], x, y, gp=gpar(fontsize=10))
    }
  } })

pdf("Shared_degs_heatmap_up.pdf", width=8,height=4)
ht+ht2
dev.off()
############################
library(tidyverse)
data<-readRDS("mutiome_final.rds")
DefaultAssay(data)<-"CTpeaks"
##
peaks <- CallPeaks(
  object = data,
  group.by = "predict.id",
  Idents = unique(data$combo_group),
  combine.peaks = TRUE,
  macs2.path = "/opt/anaconda3/envs/env01/bin/macs2",
  outdir = "/Volumes/Seagate/multiome_cp"
)
##
data$combo_group <- paste0(data$predict.id, "_", data$phenotype)
peaks <- CallPeaks(
  object = data,
  group.by = "combo_group",
  Idents = unique(data$combo_group),
  combine.peaks = TRUE,
  macs2.path = "/opt/anaconda3/envs/env01/bin/macs2",
  outdir = "/Volumes/Seagate/multiome_cp"
)

peaks_df<-as.data.frame(peaks)
#peaks_df$seqnames<- paste0('chr', '', peaks_df$seqnames)
write.csv(peaks_df, 'CTpeaks_annotated.csv')
dfc<-read_csv('CTpeaks_annotated.csv')
table(dfc$peak_called_in)
peaks_df_expanded <- dfc %>%
  separate_rows(peak_called_in, sep = ",") %>%
  separate(peak_called_in, into = c("cluster", "phenotype"), sep = "_")
## ROT
ROT_df <- peaks_df_expanded %>%
  dplyr::filter(phenotype == "ROT") %>%
  dplyr::select(seqnames, start, end, width, cluster)
ROT_df_unique <- ROT_df[!duplicated(ROT_df[, c("seqnames", "start", "end")]), ]

##Control
Control_df <- peaks_df_expanded %>%
  dplyr::filter(phenotype == "control") %>%
  dplyr::select(seqnames, start, end, width, cluster)
# remove dup
Control_df_unique <- Control_df[!duplicated(Control_df[, c("seqnames", "start", "end")])]

Common_df <- peaks_df_expanded %>%
  dplyr::filter(paste(seqnames, start, end, sep = "_") %in% paste(ROT_df$seqnames, ROT_df$start, ROT_df$end, sep = "_") &
                  paste(seqnames, start, end, sep = "_") %in% paste(Control_df$seqnames, Control_df$start, Control_df$end, sep = "_")) %>%
  dplyr::select(seqnames, start, end, width, cluster)
Common_df_unique <- Common_df[!duplicated(Common_df[, c("seqnames", "start", "end")]), ]

table(Common_df)
################ prepare for LDSC################
# 創建一個新的列來組合細胞類型和表型
peaks_df_expanded$cell_phenotype <- paste(peaks_df_expanded$cluster, peaks_df_expanded$phenotype, sep = "_")

# 獲取所有唯一的細胞類型和表型組合
cell_phenotypes <- unique(peaks_df_expanded$cell_phenotype)

# 創建一個空的數據框來存儲每個細胞類型和表型組合的區域數
counts <- data.frame(cell_phenotype = character(), count = integer())
# 對每種細胞類型和表型組合遍歷
for (cp in cell_phenotypes) {
  # 從數據框中選擇特定的細胞類型和表型組合
  df_sub <- peaks_df_expanded[peaks_df_expanded$cell_phenotype == cp, ]
  
  # 選擇你需要的列
  df_bed <- df_sub[, c("seqnames", "start", "end", "cell_phenotype")]
  
  # 將你的數據框寫入一個 .bed 文件
  write.table(df_bed, file = paste0(cp, ".bed"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  # 計算區域數並將結果添加到counts數據框中
  counts <- rbind(counts, data.frame(cell_phenotype = cp, count = nrow(df_sub)))
}

# 打印counts數據框以查看每個細胞類型和表型組合的區域數
print(counts)
###############
#接著做有common peak的三種peak files
df<-peaks_df_expanded
# 獲取所有唯一的細胞類型
celltypes <- unique(df$cluster)

# 對每種細胞類型遍歷
for (ct in celltypes) {
  # 從數據框中選擇特定的細胞類型
  df_sub_rot <- df[df$cluster == ct & df$phenotype == "ROT", c("seqnames", "start", "end")]
  df_sub_control <- df[df$cluster == ct & df$phenotype == "control", c("seqnames", "start", "end")]
  
  # 在seqnames前加上"chr"
  df_sub_rot$seqnames <- paste0("chr", df_sub_rot$seqnames)
  df_sub_control$seqnames <- paste0("chr", df_sub_control$seqnames)
  
  # 找出在ROT和Control中都有表達的區域
  common <- intersect(df_sub_rot, df_sub_control)
  
  # 找出只在ROT或Control中存在的區域
  rot_unique <- setdiff(df_sub_rot, df_sub_control)
  control_unique <- setdiff(df_sub_control, df_sub_rot)
  
  # 將你的數據框寫入一個 .bed 文件
  write.table(common, file = paste0(ct, "_common.bed"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(rot_unique, file = paste0(ct, "_ROT_unique.bed"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(control_unique, file = paste0(ct, "_Control_unique.bed"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  # 計算並打印每種類型的區域數量
  print(paste0(ct, " common: ", nrow(common)))
  print(paste0(ct, " ROT unique: ", nrow(rot_unique)))
  print(paste0(ct, " Control unique: ", nrow(control_unique)))
}
#######lift over from rn7 to hg38#####

library(rtracklayer)

# 獲取所有BED文件的列表
bed_files <- list.files('./bed_file_rn7', pattern = "*.bed", full.names = TRUE)
output_path <- "./lifted_bed_files_hg/"
# 指定chain文件的路徑
chain_path <- "rn7ToHg38.over.chain"
chain_path2 <- "hg38ToHg19.over.chain"
chain_data <- rtracklayer::import.chain(chain_path)
chain_data2 <- rtracklayer::import.chain(chain_path2)

# 對每個BED文件進行轉換
for (bed_file in bed_files) {
  # 讀取BED文件
  bed_gr <- rtracklayer::import(bed_file)
  
  # 使用liftOver進行轉換
  lifted_bed <- liftOver(bed_gr, chain_data)
  lifted_bed2 <- liftOver(lifted_bed, chain_data2)
  # 過濾掉空的GRanges對象
  filtered_bed <- lifted_bed2[sapply(lifted_bed2, function(x) length(x) > 0)]
  # 將GRangesList轉換為一個GRanges物件
  flattened_bed <- unlist(filtered_bed)
  # 將結果保存為新的BED文件
  output_file <- paste0(output_path, basename(bed_file))
  export(flattened_bed, output_file, format = "bed")
  
  # 打印每個文件的轉換統計信息
  cat(paste(basename(bed_file), ":", length(flattened_bed), "regions lifted\n"))
}



# Import the chain file
#chain <- import.chain(chain_path)
#chain2 <- import.chain(chain_path2)
# 使用liftOver進行轉換
#lifted_bed <- liftOver(bed_file, chain)
#from hg38 to hg19
#lifted_bed2 <- liftOver(lifted_bed, chain2)
# 過濾掉空的GRanges對象
#filtered_bed <- lifted_bed2[sapply(lifted_bed2, function(x) length(x) > 0)]

# 將GRangesList轉換為一個GRanges物件
#flattened_bed <- unlist(filtered_bed)
#nrow(flattened_bed)
# 儲存轉換後的BED文件
#export(flattened_bed, "test01.bed")




################
peaks_new<-peaks_df
#######製作re-analysis用的peak file ####
# 指定要寫入的列
cols <- c("seqnames", "start", "end")
# 建立檔案名稱（這裡假設是 "Control.bed"）
filename <- "Control"
# 使用 write.table 函數將資料框寫入 .bed 檔案
fasta_order <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "X", "Y", "MT", "MU150189.1", "JACYVU010000493.1", "MU150193.1", "MU150196.1", "JACYVU010000706.1", "JACYVU010000238.1", "JACYVU010000319.1", "MU150220.1", "MU150222.1", "JACYVU010000731.1", "JACYVU010000732.1", "MU150223.1", "JACYVU010000738.1", "JACYVU010000744.1", "JACYVU010000754.1", "MU150200.1", "JACYVU010000619.1", "MU150203.1", "JACYVU010000589.1", "JACYVU010000634.1", "JACYVU010000642.1", "JACYVU010000653.1", "JACYVU010000315.1", "JACYVU010000587.1", "JACYVU010000665.1")
Control_df_unique <- Control_df_unique %>% arrange(match(seqnames, fasta_order), start)

write.table(Control_df_unique[, cols],
            paste0("./aggr_CT_peaks", filename, ".bed"),
            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE, 
            sep = "\t")

# 建立檔案名稱（這裡假設是 "Control.bed"）
filename <- "ROT"
# 使用 write.table 函數將資料框寫入 .bed 檔案
fasta_order <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "X", "Y", "MT", "MU150189.1", "JACYVU010000493.1", "MU150193.1", "MU150196.1", "JACYVU010000706.1", "JACYVU010000238.1", "JACYVU010000319.1", "MU150220.1", "MU150222.1", "JACYVU010000731.1", "JACYVU010000732.1", "MU150223.1", "JACYVU010000738.1", "JACYVU010000744.1", "JACYVU010000754.1", "MU150200.1", "JACYVU010000619.1", "MU150203.1", "JACYVU010000589.1", "JACYVU010000634.1", "JACYVU010000642.1", "JACYVU010000653.1", "JACYVU010000315.1", "JACYVU010000587.1", "JACYVU010000665.1")
ROT_df_unique <- ROT_df_unique %>% arrange(match(seqnames, fasta_order), start)

write.table(ROT_df_unique[, cols],
            paste0("./aggr_CT_peaks", filename, ".bed"),
            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE, 
            sep = "\t")
########
# 拆分 'peak_called_in' 列的值為單獨的行
peaks_df_expanded <- peaks_df %>%
  separate_rows(peak_called_in, sep = ",") %>%
  separate(peak_called_in, into = c("cluster", "phenotype"), sep = "_")

# 創建只有 'ROT' 的數據框
ROT_df <- peaks_df_expanded %>%
  dplyr::filter(phenotype == "ROT") %>%
  dplyr::select(seqnames, start, end, width, cluster)
ROT_df_unique <- ROT_df[!duplicated(ROT_df[, c("seqnames", "start", "end")]), ]

# 創建只有 'Control' 的數據框
Control_df <- peaks_df_expanded %>%
  dplyr::filter(phenotype == "control") %>%
  dplyr::select(seqnames, start, end, width, cluster)
# 移除重複的行
Control_df_unique <- Control_df[!duplicated(Control_df[, c("seqnames", "start", "end")]), ]

# 如果 "seqnames"、"start" 和 "end" 的組合可以唯一識別每個峰值，則使用此方法創建 'Common_df'
Common_df <- peaks_df_expanded %>%
  dplyr::filter(paste(seqnames, start, end, sep = "_") %in% paste(ROT_df$seqnames, ROT_df$start, ROT_df$end, sep = "_") &
                  paste(seqnames, start, end, sep = "_") %in% paste(Control_df$seqnames, Control_df$start, Control_df$end, sep = "_")) %>%
  dplyr::select(seqnames, start, end, width, cluster)

ROT_df$seqnames <- paste("chr", ROT_df$seqnames, sep = "")
Control_df$seqnames <- paste("chr", Control_df$seqnames, sep = "")
Common_df$seqnames <- paste("chr", Common_df$seqnames, sep = "")




write_bed <- function(peaks_df_expanded, cell_type, filename) {
  df_filtered <- peaks_df_expanded[peaks_df_expanded['cluster'] == cell_type, ]
  df_filtered$seqnames <- paste0("", df_filtered$seqnames)
  write.table(df_filtered[, c("seqnames", "start", "end")],
              paste0("./LDSC_aggr_bed/", filename, ".bed"),
              quote = F, row.names = F, col.names = F, sep = "\t")
}

cell_types <- c("Astrocyte", "Microglia", "OPC", "Olig", "Drd1-MSN","Drd2-MSN", "Glutamatergic", "Interneuron")
dataframes <- list(ROT = ROT_df, Control = Control_df, Common = Common_df)

for(cell_type in cell_types) {
  for(group in names(dataframes)) {
    write_bed(dataframes[[group]], cell_type, paste0(cell_type, "_", group))
  }
}









