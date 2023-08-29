
# This script performs data normalization, differential gene and accessibility analysis, and heatmap visualization.

#########################
# Set the default assay to RNA and normalize the data
DefaultAssay(data)<-"RNA"
data<-NormalizeData(data, assay="RNA") %>% ScaleData(features=rownames(data))

# Compute average expression per cell type and phenotype
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

# Perform differential gene expression analysis between ROT and control for each cell type

########################################
########## DO DEG ROT/Control ##########
########################################

DefaultAssay(data)<-"RNA"
all.genes<-rownames(data)
data<-NormalizeData(data) 

# Subset data by cell type
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
#significant<-read.csv("mutiome_DEGs_MAST_ROTCtrl.csv")

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

# Perform differential accessibility analysis between ROT and control for each cell type

########################################
########## DO DAE ROT/Control ##########
########################################

DefaultAssay(data)<-"ATAC"
all.genes<-rownames(data)
data<-NormalizeData(data) 

# Subset data by cell type
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

noDup<-noDup[order(noDup$celltype, noDup$avg_log2FC),]
#lab<-which(noDup$gene %in% c(ad_related$Gene,"PLCG2","MAPT","ARL17B","SAMD4A","PTPRG","MDGA2","GPR158") & abs(noDup$avg_log2FC)>0.5)
ht=Heatmap(mat, cluster_rows=F,show_row_dend=F,cluster_columns=F,col=colorRamp2(c(-1,0,1),c("blue","white","orange")), name="MAST-adj log2FC",  show_column_names=T, show_row_names=F, column_title=NULL,row_names_gp = gpar(fontsize = 3), top_annotation=ha)
ht=Heatmap(mat, cluster_rows=F,show_row_dend=F,cluster_columns=F,col=colorRamp2(c(-1,0,1),c("blue","white","orange")), name="MAST-adj log2FC",  show_column_names=T, show_row_names=F, column_title=NULL,row_names_gp = gpar(fontsize = 3))

# Generate heatmaps for differential gene expression and accessibility
pdf("ROT_DAGs_MAST_heatmap_foldCs.pdf", width=6, height=8)
ht
dev.off()

