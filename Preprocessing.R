# Load required libraries
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

# Clear workspace
rm(list = ls())

# Read 10X data
data <- Read10X(data.dir = "/Volumes/Seagate/multiome_cp/sample_aggr/filtered_feature_bc_matrix/")
rna_counts <- data$`Gene Expression`
atac_counts <- data$Peaks

# Create a Seurat object
data <- CreateSeuratObject(counts = rna_counts, project="all")
numCells <- nrow(data@meta.data)

# Add ATAC assay and perform initial preprocessing
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

# Add phenotype information based on the last digit of the cell identifier
library(dplyr)
library(stringr)
orig.ident <-rownames(data@meta.data)
last.digit <- as.numeric(str_extract(orig.ident, "-\\d+$"))
data$phenotype <- ifelse(last.digit %in% c(-1, -2), "control", "ROT")

# Filtering, normalization, scaling, and dimensionality reduction
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


# Doublet filtering
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

# Further preprocessing and clustering
data <- SCTransform(data,assay="SCT",vars.to.regress = c("percent.mt"),verbose=F) %>% RunPCA(ndims=30) %>% 
  FindNeighbors(dims = 1:30) %>% RunUMAP(dims = 1:30, reduction.name="umap.rna",reduction.key = "rnaUMAP_")
data<-FindMultiModalNeighbors(data, reduction.list=list("pca","lsi"), dims.list=list(1:30,2:50))
data<-RunUMAP(data,nn.name="weighted.nn", reduction.name="wnn.umap",reduction.key="wnnUMAP_")
data<-FindClusters(data, graph.name="wsnn", algorithm=3)

pdf("rna_atac.umap.pdf", width=14, height=6)
p1=DimPlot(data, reduction = "umap.rna") + ggtitle("RNA")+theme(legend.position="none")
p2=DimPlot(data, reduction = "umap.atac") + ggtitle("ATAC")+theme(legend.position="none")
p3=DimPlot(data, reduction="wnn.umap", label = TRUE)+ ggtitle("WNN")
p1+p2+p3
dev.off()
Idents(data)<-data$phenotype
DimPlot(data, reduction="wnn.umap", label = TRUE)


# Find markers for each cluster
DefaultAssay(data) <- "RNA"
Idents(data)<-data$seurat_clusters
all.markers<-FindAllMarkers(data, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5, test.use = "MAST")
dim(all.markers)
table(all.markers$cluster)
top3_markers <- as.data.frame(all.markers %>% group_by(cluster)%>%top_n(n=3, wt = avg_log2FC))
top3_markers

# Find conserved markers across different conditions within specific clusters
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

# Naming the clusters

###############################
###### naming the cluster #####
###############################

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


#########################
###### Peak-calling #####
#########################

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
