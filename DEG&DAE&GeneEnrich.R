
# This script performs data normalization, differential gene and accessibility analysis,GO funcitonal pathway analysis and heatmap visualization.

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
bar1<-HeatmapAnnotation(up=anno_barplot(tab[,2], gp=gpar(fill="red"), ylim=c(0,280)),
                        down=anno_barplot((tab[,1] * -1), gp=gpar(fill="blue"), ylim=c(-280,0),
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

# Additional Figure
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

#######################################
#########Gene set enrichment###########
#######################################

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

