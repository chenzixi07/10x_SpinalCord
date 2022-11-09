##### 20210807_s4hV6Sub_Monocle3_monPN.pdf
monpn.col <- c(
  "#6ab662",
  "#377EB8",
  "#cab3d6",
  "#2ea37e",
  "#ebc4b5",
  "#e11517",
  "#efd092",
  "#c47e70",
  "#b8e294",
  "#774ea3",
  "#ffd92f",
  "#b5d7e8",
  "#685f7e"
)

p <- monocle3::plot_cells(cds.monpn, 
                     label_cell_groups = F, label_branch_points = F,
                     label_principal_points = F, label_roots = F, 
                     label_leaves = F,show_trajectory_graph = F) +
  scale_color_manual(values = monpn.col) + coord_fixed()

ggsave("20211213_s4hV6Sub_Monocle3_monPN.pdf",p,width = 7.5,height = 4.5)

##### 1.	heatmap top5做细胞热图，cluster按照顺序：
#5，9，12，23，11，1，3，7，  10，19，  22， 6， 14，13，2， 21，4，15， 18，17， 25，16， 8，24，20
final_col <- c(
  "#cdd3b4",
  "#1f78b3",
  "#a7be8e",
  "#349f2d",
  "#ffe0ed",
  "#e2191b",
  "#acb564",
  "#cbb2d5",
  "#f2abd2",
  "#e6ab02",
  "#feff9a",
  "#f077a3",
  "#a6cee4",
  "#387fb8",
  "#65a61f",
  "#fdbf6f",
  "#a7761d",
  "#e6ab02",
  "#a7761d",
  "#7571b3",
  "#549336",
  "#b05928",
  "#dc0077",
  "#984da2",
  "#ff7f00"
)

Seurat::DimPlot(project)
project@active.ident
project.level1 <- levels(project@active.ident)
project@active.ident <- factor(project@active.ident, 
                               levels = c(5,9,12,23,11,1,3,
                                          7,10,19,22,6,14,13,2,
                                          21,4,15,18,17,25,16,8,24,20))

marker.project <- read.table("../20210707_s4h_v6/20210707.s4hV6.AllMarker_wilcox.xls",header = T,sep = "\t")
marker.project$cluster <- factor(marker.project$cluster, levels = c(5,9,12,23,11,1,3,
                                                                    7,10,19,22,6,14,13,2,
                                                                    21,4,15,18,17,25,16,8,24,20))
heatmarker <- marker.project  %>% group_by(cluster) %>% top_n(5,avg_log2FC) 
dim(heatmarker)

pos <- c()
for (i in (match(c(5,9,12,23,11,1,3,7,10,19,22,6,14,13,2,21,4,15,18,17,25,16,8,24,20),
                heatmarker$cluster))){
  pos <- c(pos,c(i:(i+4)))
}
pos

heatmarker <- heatmarker[pos,]
library(viridis)

col <- final_col[c(5,9,12,23,11,1,3,7,10,19,22,6,14,13,2,21,4,15,18,17,25,16,8,24,20)]

pdf("20211213.s4hV6.heatmap.pdf",width = 15,height = 15)
  Seurat::DoHeatmap(project,features = heatmarker$gene, group.colors = col) + scale_fill_viridis(option = "D") 
dev.off()


pdf("20211221.s4hV6.heatmap.pdf",width = 15,height = 6)
  #Seurat::DoHeatmap(project,features = heatmarker$gene, group.colors = col) + scale_fill_viridis(option = "D")
  Seurat::DoHeatmap(project,features = rev(heatmarker$gene), group.colors = col) + scale_fill_viridis(option = "D") 
dev.off()

###### PG曲线，然后颜色和cluster吻合。Cluster群颜色的文件名：20211206.s4hV6Sub.PGsplit.Cluster.res0.2.pdf
library(RColorBrewer)
cols=c(brewer.pal(12,"Paired"),brewer.pal(8,"Dark2"),
       brewer.pal(9,"Set1"),brewer.pal(8,"Set2"),
       brewer.pal(12,"Set3"),brewer.pal(8,"Pastel2"),
       brewer.pal(9,"Pastel1"),brewer.pal(8,"Accent"))
col=unique(cols)
color <- col[c(1:8)][c(1,6,3,4,5,2,7,8)]

marker <- read.table("20211213.marker1.txt", header = F)
cds.plot <- cds.pgsplit[rowData(cds.pgsplit)$gene_short_name %in% marker$V1,]
p<-monocle3::plot_genes_in_pseudotime(cds.plot, color_cells_by = "Seurat_res0.2", ncol = 4, cell_size = 1) + 
  scale_color_manual(values = color)
p
ggsave("20211213.s4hV6Sub.PGsplit.curve.pdf",p, width = 16, height = ceiling(length(marker$V1)/4)*2)

##### 202112015 monpn + mouse ref
marker <- read.table("20211206.GSEA.ref.txt", header = F)

p<- monocle3::plot_genes_by_group(cds.monpn, unique(marker$V1))
p
ggsave("20211215.s4hV6Sub.monpn.mouseRef.pdf",p, width = 16, height = ceiling(length(marker$V1)/4),limitsize = F)

##### D+V project.insplit to cds.in.0912
project.insplit@meta.data$DV <- "Mix"
project.insplit@meta.data$DV[project.insplit@meta.data$merge.split %in% c(20,16,17,8,21,2,6)] <- "Dorsal"
project.insplit@meta.data$DV[project.insplit@meta.data$merge.split %in% c(18,1,28,12,19)] <- "Ventral"

cds.in.0912@colData$DV <- "pIN"
Dorsal_id <- rownames(project.insplit@meta.data)[project.insplit@meta.data$merge.split %in% c(20,16,17,8,21,2,6)] 
Ventral_id <- rownames(project.insplit@meta.data)[project.insplit@meta.data$merge.split %in% c(18,1,28,12,19)]
Mix_id <- rownames(project.insplit@meta.data)[project.insplit@meta.data$merge.split %in% c(22)]

cds.in.0912@colData$DV[rownames(cds.in.0912@colData) %in% Dorsal_id] <- "Dorsal"
cds.in.0912@colData$DV[rownames(cds.in.0912@colData) %in% Ventral_id] <- "Ventral"
cds.in.0912@colData$DV[rownames(cds.in.0912@colData) %in% Mix_id] <- "Mix"


p1<-plot_cells(cds.in.0912,color_cells_by = "DV",
           label_cell_groups = F, label_branch_points = F,
           label_principal_points = F, label_roots = F, 
           label_leaves = F,show_trajectory_graph = F) + coord_fixed() + scale_color_manual(values = monpn.col[c(2,3,5,6)])
p1 
p2 <-plot_cells(cds.in.0912,color_cells_by = "cluster",
           label_cell_groups = F, label_branch_points = F,
           label_principal_points = F, label_roots = F, 
           label_leaves = F,show_trajectory_graph = F) + coord_fixed() 

p1+p2

##### D+V project.insplit SCENIC

##### 20211220 monpn
##### 20210807_s4hV6Sub_Monocle3_monPN.pdf
monpn.col <- c(
  "#6ab662",
  "#377EB8",
  "#cab3d6",
  "#2ea37e",
  "#ebc4b5",
  "#e11517",
  "#efd092",
  "#c47e70",
  "#b8e294",
  "#774ea3",
  "#ffd92f",
  "#b5d7e8",
  "#685f7e"
)

#
monpn.col1 <- rep("#D9D9D9",13)
monpn.col1[c(10,3,13,2,12)] <- monpn.col[c(10,3,13,2,12)]

p1 <- monocle3::plot_cells(cds.monpn, 
                          label_cell_groups = F, label_branch_points = F,
                          label_principal_points = F, label_roots = F, 
                          label_leaves = F,show_trajectory_graph = F) +
  scale_color_manual(values = monpn.col1) + coord_fixed()

#
monpn.col2 <- rep("#D9D9D9",13)
monpn.col2[c(10,6,4)] <- monpn.col[c(10,6,4)]

p2 <- monocle3::plot_cells(cds.monpn, 
                           label_cell_groups = F, label_branch_points = F,
                           label_principal_points = F, label_roots = F, 
                           label_leaves = F,show_trajectory_graph = F) +
  scale_color_manual(values = monpn.col2) + coord_fixed()
p2

#
monpn.col3 <- rep("#D9D9D9",13)
monpn.col3[c(11,1,9,5,7,8)] <- monpn.col[c(11,1,9,5,7,8)]

p3 <- monocle3::plot_cells(cds.monpn, 
                           label_cell_groups = F, label_branch_points = F,
                           label_principal_points = F, label_roots = F, 
                           label_leaves = F,show_trajectory_graph = F) +
  scale_color_manual(values = monpn.col3) + coord_fixed()
p3

pdf("20211220_s4hV6Sub_Monocle3_monPN.pdf",width = 7.5,height = 4.5)
print(p1)
print(p2)
print(p3)
dev.off()

##### 20211221 PGsplit featureplot
marker <- read.table("20211221.PGsplit.feature.marker.txt", header = F, sep = "\t")

pdf("20211221.s4hV6Sub.PGsplit.FeaturePlot.pdf", width = 24, height = ceiling(length(marker$V1)))
Seurat::FeaturePlot(project.pgsplit, reduction = "umap",
                    features = as.vector(unique(marker$V1)),
                    keep.scale = "all", ncol = 4,
                    cols = c("lightgrey", "red"))
dev.off()

##### 20211221 PGspslit monocle2 pseudotime heatmap
library(ComplexHeatmap)
library(circlize)

#pMN
df <- data.frame(Cluster = colData(cds.pgsplit)$Seurat_res0.2[match(names(pseudotime(cds.pgsplit)[order(pseudotime(cds.pgsplit))]),
                                                                    names(colData(cds.pgsplit)$Seurat_res0.2))], 
                 Pseudotime = pseudotime(cds.pgsplit)[order(pseudotime(cds.pgsplit))])
df <- df[df$Cluster %in% c(2,3,4,5,7,8),]

ClusterCol <- brewer.pal(8,"Paired")[c(6,3,4,5,7,8)]
#names(ClusterCol) <- sort(unique(colData(cds.pgsplit)$Seurat_res0.2))
names(ClusterCol) <- c(2,3,4,5,7,8)
ClusterCol

max_pseudo <- ceiling(max(pseudotime(cds.pgsplit)))

ha <- HeatmapAnnotation(df = df,col = list(Cluster = ClusterCol, 
                                           Pseudotime = colorRamp2(seq(from=0,to=max_pseudo,length=11),rev(brewer.pal(11, "RdBu")))),
                        annotation_name_side = "left")


df1 <- subset(df, Cluster %in% c(2,3,4,5,7,8), select=Cluster )
df2 <- subset(df, Cluster %in% c(2,3,4,5,7,8), select=Pseudotime )

ha1 <- HeatmapAnnotation(df = df1, col = list(Cluster = ClusterCol), 
                         annotation_name_side = "left")
ha2 <- HeatmapAnnotation(df = df2, col = list(Pseudotime = colorRamp2(seq(from=0,to=max_pseudo,length=11),rev(brewer.pal(11, "RdBu")))), 
                         annotation_name_side = "left")


pMN.genes <- read.table("20211221.pseudo.pMN.txt", header = T, sep = "\t")

for (tf in (unique(pMN.genes$TF))){
  genes <- subset(pMN.genes, TF == tf, select = Gene)$Gene
  genes <- genes[genes %in% rownames(rowData(cds.pgsplit))]
  
  pdf(paste("20211221_s4hV6Sub_Monocle3PGsplit.pMN.",tf,".heatmap.pdf",sep=""),width=12,height = length(genes)/5)
  pt.matrix <- exprs(cds.pgsplit)[match(genes,rownames(rowData(cds.pgsplit))),
                                  order(pseudotime(cds.pgsplit)[names(pseudotime(cds.pgsplit)) %in% rownames(df)])]
  pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
  pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
  rownames(pt.matrix) <- genes

#Ward.D2 Hierarchical Clustering
hthc <- Heatmap(
  as.matrix(pt.matrix),
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  top_annotation               = ha1,
  bottom_annotation            = ha2
)

hthc
  print(hthc)
dev.off()
}

#pIN
df <- data.frame(Cluster = colData(cds.pgsplit)$Seurat_res0.2[match(names(pseudotime(cds.pgsplit)[order(pseudotime(cds.pgsplit))]),
                                                                    names(colData(cds.pgsplit)$Seurat_res0.2))], 
                 Pseudotime = pseudotime(cds.pgsplit)[order(pseudotime(cds.pgsplit))])
df <- df[df$Cluster %in% c(1,3,4,5,6,7,8),]

ClusterCol <- brewer.pal(8,"Paired")[c(1,3,4,5,2,7,8)]
#names(ClusterCol) <- sort(unique(colData(cds.pgsplit)$Seurat_res0.2))
names(ClusterCol) <- c(1,3,4,5,6,7,8)
ClusterCol

max_pseudo <- ceiling(max(pseudotime(cds.pgsplit)))

ha <- HeatmapAnnotation(df = df,col = list(Cluster = ClusterCol, 
                                           Pseudotime = colorRamp2(seq(from=0,to=max_pseudo,length=11),rev(brewer.pal(11, "RdBu")))),
                        annotation_name_side = "left")


df1 <- subset(df, Cluster %in% c(1,3,4,5,6,7,8), select=Cluster )
df2 <- subset(df, Cluster %in% c(1,3,4,5,6,7,8), select=Pseudotime )

ha1 <- HeatmapAnnotation(df = df1, col = list(Cluster = ClusterCol), 
                         annotation_name_side = "left")
ha2 <- HeatmapAnnotation(df = df2, col = list(Pseudotime = colorRamp2(seq(from=0,to=max_pseudo,length=11),rev(brewer.pal(11, "RdBu")))), 
                         annotation_name_side = "left")


pIN.genes <- read.table("20211221.pseudo.pIN.txt", header = T, sep = "\t")

for (tf in (unique(pIN.genes$TF))){
  genes <- subset(pIN.genes, TF == tf, select = Gene)$Gene
  genes <- genes[genes %in% rownames(rowData(cds.pgsplit))]
  
  pdf(paste("20211221_s4hV6Sub_Monocle3PGsplit.pIN.",tf,".heatmap.pdf",sep=""),width=12,height = length(genes)/5)
  
  pt.matrix <- exprs(cds.pgsplit)[match(genes,rownames(rowData(cds.pgsplit))),
                                  order(pseudotime(cds.pgsplit)[names(pseudotime(cds.pgsplit)) %in% rownames(df)])]
  pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
  pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
  rownames(pt.matrix) <- genes
  
  #Ward.D2 Hierarchical Clustering
  hthc <- Heatmap(
    as.matrix(pt.matrix),
    name                         = "z-score",
    col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
    show_row_names               = TRUE,
    show_column_names            = FALSE,
    row_names_gp                 = gpar(fontsize = 6),
    clustering_method_rows = "ward.D2",
    clustering_method_columns = "ward.D2",
    row_title_rot                = 0,
    cluster_rows                 = TRUE,
    cluster_row_slices           = FALSE,
    cluster_columns              = FALSE,
    top_annotation               = ha1,
    bottom_annotation            = ha2
  )
  
  hthc
  print(hthc)
  dev.off()
}


##### 20211221 Curve repeat
marker <- unique(read.table("20211029_PGsplit_marker1.txt", header = F)$V1)
pt.matrix <- exprs(cds.pgsplit)[match(marker,rownames(rowData(cds.pgsplit))),]
#pt.matrix <- project.pgsplit@assays$RNA@scale.data[rownames(project.pgsplit@assays$RNA@scale.data) %in% marker$V1,]
pt.matrix <- as.data.frame(t(pt.matrix))
dim(pt.matrix)
head(pt.matrix)

pt.matrix$pseudotime <- pseudotime(cds.pgsplit)
pt.matrix$cluster <- colData(cds.pgsplit)$Seurat_res0.2
#pt.matrix$sample <- pData(cds.pgsplit)$orig.ident
pt.matrix$celltype <- "Progenitor"
pt.matrix$celltype[pt.matrix$cluster %in% c(1,6)] <- "pIN"
pt.matrix$celltype[pt.matrix$cluster == 2] <- "pMN"

#pt.matrix <- reshape2::melt(pt.matrix)
head(pt.matrix)

#data <- pt.matrix[pt.matrix$variable %in% marker$V1,]
data <- pt.matrix[pt.matrix$celltype %in% c("pIN","pMN"),]
data <- reshape2::melt(data, id=c("pseudotime", "cluster","celltype"))
data$celltype <- factor(data$celltype, levels = c("pMN","pIN"))

p<-ggplot(data) + geom_point(aes(x=pseudotime, y=value, color = celltype),size = 0.1) +
  geom_smooth(aes(x=pseudotime, y=value, color = celltype),size = 0.5) + scale_color_brewer(palette = "Set1") +
  theme_bw()+  facet_wrap(.~variable,ncol=4,scales = "free_y") +  scale_y_log10() +
  theme(panel.grid = element_blank(), strip.background = element_blank()) +
  xlab("Pesudotime") + ylab("Expression/log10")
p

pdf("20211029.PGsplit.pseudotime_curve1.pdf",width=10,height = length(marker)/4)
print(p)
dev.off()

##### 20211029 marker2 curve
marker <- unique(read.table("20211029_PGsplit_marker2.txt", header = F)$V1)
pt.matrix <- exprs(cds.pgsplit)[match(marker,rownames(rowData(cds.pgsplit))),]
#pt.matrix <- project.pgsplit@assays$RNA@scale.data[rownames(project.pgsplit@assays$RNA@scale.data) %in% marker$V1,]
pt.matrix <- as.data.frame(t(pt.matrix))
dim(pt.matrix)
head(pt.matrix)

pt.matrix$pseudotime <- pseudotime(cds.pgsplit)
pt.matrix$cluster <- colData(cds.pgsplit)$Seurat_res0.2
#pt.matrix$sample <- pData(cds.pgsplit)$orig.ident
pt.matrix$celltype <- "Progenitor"
pt.matrix$celltype[pt.matrix$cluster %in% c(1,6)] <- "pIN"
pt.matrix$celltype[pt.matrix$cluster == 2] <- "pMN"

#pt.matrix <- reshape2::melt(pt.matrix)
head(pt.matrix)

#data <- pt.matrix[pt.matrix$variable %in% marker$V1,]
data <- pt.matrix[pt.matrix$celltype %in% c("pIN","pMN"),]
data <- reshape2::melt(data, id=c("pseudotime", "cluster","celltype"))
data$celltype <- factor(data$celltype, levels = c("pMN","pIN"))

p<-ggplot(data) + geom_point(aes(x=pseudotime, y=value, color = celltype),size = 0.1) +
  geom_smooth(aes(x=pseudotime, y=value, color = celltype),size = 0.5) + scale_color_brewer(palette = "Set1") +
  theme_bw()+  facet_wrap(.~variable,ncol=4,scales = "free_y") +  scale_y_log10() +
  theme(panel.grid = element_blank(), strip.background = element_blank()) +
  xlab("Pesudotime") + ylab("Expression/log10")

pdf("20211029.PGsplit.pseudotime_curve2.pdf",width=10,height = length(marker)/4)
print(p)
dev.off()

##### 20211221 ggalluvial
library(ggalluvial)
library(export)

# pMN
data <- read.table("20211221.pMN.circle.txt", header = T, sep = "\t")

p1<-ggplot(data,
       aes(y = Count, axis1 = TF, axis2 = Gene)) +
  geom_alluvium(aes(fill = Source), width = 1/12) +
  geom_stratum(width = 1/10) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),size=3) +
  scale_x_discrete(limits = c("TF", "Gene"), expand = c(.05, .05)) +
  scale_fill_brewer(type = "qual", palette = "Set1") + theme_minimal() + ylab(NULL) +
  theme(panel.grid = element_blank(), axis.text.y = element_blank(), axis.text.x = element_text(color="black"))+
  ggtitle("Switch genes driving the differentiation in MN")

p2<-ggplot(data,
           aes(y = Count, axis1 = TF, axis2 = Gene)) +
  geom_alluvium(aes(fill = Source), width = 1/12) +
  geom_stratum(width = 1/10) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),size=0) +
  scale_x_discrete(limits = c("TF", "Gene"), expand = c(.05, .05)) +
  scale_fill_brewer(type = "qual", palette = "Set1") + theme_minimal() + ylab(NULL) +
  theme(panel.grid = element_blank(), axis.text.y = element_blank(), axis.text.x = element_text(color="black"))+
  ggtitle("Switch genes driving the differentiation in MN")

p1+p2

export::graph2ppt(file="20211221.pMN.alluvial.txt",width=10, height=10)

# pIN
data <- read.table("20211221.pIN.circle.txt", header = T, sep = "\t")

p1<-ggplot(data,
           aes(y = Count, axis1 = TF, axis2 = Gene)) +
  geom_alluvium(aes(fill = Source), width = 1/12) +
  geom_stratum(width = 1/10) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),size=3) +
  scale_x_discrete(limits = c("TF", "Gene"), expand = c(.05, .05)) +
  scale_fill_brewer(type = "qual", palette = "Set1") + theme_minimal() + ylab(NULL) +
  theme(panel.grid = element_blank(), axis.text.y = element_blank(), axis.text.x = element_text(color="black"))+
  ggtitle("Switch genes driving the differentiation in IN")

p2<-ggplot(data,
           aes(y = Count, axis1 = TF, axis2 = Gene)) +
  geom_alluvium(aes(fill = Source), width = 1/12) +
  geom_stratum(width = 1/10) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),size=0) +
  scale_x_discrete(limits = c("TF", "Gene"), expand = c(.05, .05)) +
  scale_fill_brewer(type = "qual", palette = "Set1") + theme_minimal() + ylab(NULL) +
  theme(panel.grid = element_blank(), axis.text.y = element_blank(), axis.text.x = element_text(color="black"))+
  ggtitle("Switch genes driving the differentiation in IN")

p1+p2

export::graph2ppt(file="20211221.pIN.alluvial.txt",width=10, height=10)

##### 20211223 pMN pIN rainbow heatmap
library(ComplexHeatmap)
library(circlize)
library(monocle3)
library(RColorBrewer)

#pMN
df <- data.frame(Cluster = colData(cds.pgsplit)$Seurat_res0.2[match(names(pseudotime(cds.pgsplit)[order(pseudotime(cds.pgsplit))]),
                                                                    names(colData(cds.pgsplit)$Seurat_res0.2))], 
                 Pseudotime = pseudotime(cds.pgsplit)[order(pseudotime(cds.pgsplit))])
df <- df[df$Cluster %in% c(2,3,4,5,7,8),]

ClusterCol <- brewer.pal(8,"Paired")[c(6,3,4,5,7,8)]
#names(ClusterCol) <- sort(unique(colData(cds.pgsplit)$Seurat_res0.2))
names(ClusterCol) <- c(2,3,4,5,7,8)
ClusterCol

max_pseudo <- ceiling(max(pseudotime(cds.pgsplit)))

ha <- HeatmapAnnotation(df = df,col = list(Cluster = ClusterCol, 
                                           Pseudotime = colorRamp2(seq(from=0,to=max_pseudo,length=11),rev(brewer.pal(11, "RdBu")))),
                        annotation_name_side = "left")


df1 <- subset(df, Cluster %in% c(2,3,4,5,7,8), select=Cluster )
df2 <- subset(df, Cluster %in% c(2,3,4,5,7,8), select=Pseudotime )

ha1 <- HeatmapAnnotation(df = df1, col = list(Cluster = ClusterCol), 
                         annotation_name_side = "left")
ha2 <- HeatmapAnnotation(df = df2, col = list(Pseudotime = colorRamp2(seq(from=0,to=max_pseudo,length=11),rev(brewer.pal(11, "RdBu")))), 
                         annotation_name_side = "left")


pMN.genes <- read.table("20211223.pMN.rainbow.heatmap.txt", header = T, sep = "\t")

for (tf in (unique(pMN.genes$TF))){
  genes <- subset(pMN.genes, TF == tf, select = Gene)$Gene
  genes <- genes[genes %in% rownames(rowData(cds.pgsplit))]
  
  pdf(paste("20211221_s4hV6Sub_Monocle3PGsplit.pMN.",tf,".heatmap.pdf",sep=""),width=12,height = length(genes)/3)
  pt.matrix <- exprs(cds.pgsplit)[match(genes,rownames(rowData(cds.pgsplit))),
                                  order(pseudotime(cds.pgsplit)[names(pseudotime(cds.pgsplit)) %in% rownames(df)])]
  pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
  pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
  rownames(pt.matrix) <- genes

#Ward.D2 Hierarchical Clustering
hthc <- Heatmap(
  as.matrix(pt.matrix),
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  top_annotation               = ha1,
  bottom_annotation            = ha2
)

#hthc
  print(hthc)
dev.off()
}

#pIN
df <- data.frame(Cluster = colData(cds.pgsplit)$Seurat_res0.2[match(names(pseudotime(cds.pgsplit)[order(pseudotime(cds.pgsplit))]),
                                                                    names(colData(cds.pgsplit)$Seurat_res0.2))], 
                 Pseudotime = pseudotime(cds.pgsplit)[order(pseudotime(cds.pgsplit))])
df <- df[df$Cluster %in% c(1,3,4,5,6,7,8),]

ClusterCol <- brewer.pal(8,"Paired")[c(1,3,4,5,2,7,8)]
#names(ClusterCol) <- sort(unique(colData(cds.pgsplit)$Seurat_res0.2))
names(ClusterCol) <- c(1,3,4,5,6,7,8)
ClusterCol

max_pseudo <- ceiling(max(pseudotime(cds.pgsplit)))

ha <- HeatmapAnnotation(df = df,col = list(Cluster = ClusterCol, 
                                           Pseudotime = colorRamp2(seq(from=0,to=max_pseudo,length=11),rev(brewer.pal(11, "RdBu")))),
                        annotation_name_side = "left")


df1 <- subset(df, Cluster %in% c(1,3,4,5,6,7,8), select=Cluster )
df2 <- subset(df, Cluster %in% c(1,3,4,5,6,7,8), select=Pseudotime )

ha1 <- HeatmapAnnotation(df = df1, col = list(Cluster = ClusterCol), 
                         annotation_name_side = "left")
ha2 <- HeatmapAnnotation(df = df2, col = list(Pseudotime = colorRamp2(seq(from=0,to=max_pseudo,length=11),rev(brewer.pal(11, "RdBu")))), 
                         annotation_name_side = "left")


pIN.genes <- read.table("20211223.pIN.rainbow.heatmap.txt", header = T, sep = "\t")

for (tf in (unique(pIN.genes$TF))){
  genes <- subset(pIN.genes, TF == tf, select = Gene)$Gene
  genes <- genes[genes %in% rownames(rowData(cds.pgsplit))]
  
  pdf(paste("20211221_s4hV6Sub_Monocle3PGsplit.pIN.",tf,".heatmap.pdf",sep=""),width=12,height = length(genes)/4)
  
  pt.matrix <- exprs(cds.pgsplit)[match(genes,rownames(rowData(cds.pgsplit))),
                                  order(pseudotime(cds.pgsplit)[names(pseudotime(cds.pgsplit)) %in% rownames(df)])]
  pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
  pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
  rownames(pt.matrix) <- genes
  
  #Ward.D2 Hierarchical Clustering
  hthc <- Heatmap(
    as.matrix(pt.matrix),
    name                         = "z-score",
    col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
    show_row_names               = TRUE,
    show_column_names            = FALSE,
    row_names_gp                 = gpar(fontsize = 6),
    clustering_method_rows = "ward.D2",
    clustering_method_columns = "ward.D2",
    row_title_rot                = 0,
    cluster_rows                 = TRUE,
    cluster_row_slices           = FALSE,
    cluster_columns              = FALSE,
    top_annotation               = ha1,
    bottom_annotation            = ha2
  )
  
  #hthc
  print(hthc)
  dev.off()
}

##### 20211223 pIN curve
marker <- unique(read.table("20211223.pIN.rainbow.heatmap.txt", header = T, sep = "\t")$Gene)
pt.matrix <- exprs(cds.pgsplit)[match(marker,rownames(rowData(cds.pgsplit))),]
#pt.matrix <- project.pgsplit@assays$RNA@scale.data[rownames(project.pgsplit@assays$RNA@scale.data) %in% marker$V1,]
pt.matrix <- as.data.frame(t(pt.matrix))
dim(pt.matrix)
head(pt.matrix)

pt.matrix$pseudotime <- pseudotime(cds.pgsplit)
pt.matrix$cluster <- colData(cds.pgsplit)$Seurat_res0.2
#pt.matrix$sample <- pData(cds.pgsplit)$orig.ident
pt.matrix$celltype <- "Progenitor"
pt.matrix$celltype[pt.matrix$cluster %in% c(1,6)] <- "pIN"
pt.matrix$celltype[pt.matrix$cluster == 2] <- "pMN"

#pt.matrix <- reshape2::melt(pt.matrix)
head(pt.matrix)

#data <- pt.matrix[pt.matrix$variable %in% marker$V1,]
data <- pt.matrix[pt.matrix$celltype %in% c("pIN","pMN"),]
data <- reshape2::melt(data, id=c("pseudotime", "cluster","celltype"))
data$celltype <- factor(data$celltype, levels = c("pMN","pIN"))

p<-ggplot(data) + geom_point(aes(x=pseudotime, y=value, color = celltype),size = 0.5) +
  geom_smooth(aes(x=pseudotime, y=value, color = celltype)) + scale_color_brewer(palette = "Set1") +
  theme_bw()+  facet_wrap(.~variable,ncol=4,scales = "free_y") +  scale_y_log10() +
  theme(panel.grid = element_blank(), strip.background = element_blank()) +
  xlab("Pesudotime") + ylab("Expression/log10")

pdf("20211223.PGsplit.pIN.pseudotime_curve1.pdf",width=10,height = length(marker)/4)
print(p)
dev.off()

##### 20211226 monpn 10，3，13，2，12 curve
marker <- read.table("20211226.monpn.curve1.txt", header = F)$V1
pt.matrix <- exprs(cds.mn.0912)[match(marker,rownames(rowData(cds.mn.0912))),]
#pt.matrix <- project.pgsplit@assays$RNA@scale.data[rownames(project.pgsplit@assays$RNA@scale.data) %in% marker$V1,]
pt.matrix <- as.data.frame(t(pt.matrix))
dim(pt.matrix)
head(pt.matrix)

pt.matrix$pseudotime <- pseudotime(cds.mn.0912)
pt.matrix$cluster <- colData(cds.mn.0912)$clusters
data <- reshape2::melt(pt.matrix, id=c("pseudotime", "cluster"))
data$cluster <- factor(data$cluster, levels = c(2,3,10,12,13))
head(data)

col <- monpn.col[c(2,3,10,12,13)]

p<-ggplot(data) + geom_point(aes(x = pseudotime, y = value, color = cluster),size = 0.5) +
  geom_smooth(aes(x=pseudotime, y=value, color = cluster)) + scale_color_manual(values = col) +
  theme_bw()+  facet_wrap(.~variable,ncol=4,scales = "free_y") +  scale_y_log10() +
  theme(panel.grid = element_blank(), strip.background = element_blank()) +
  xlab("Pesudotime") + ylab("Expression/log10")
p

pdf("20211226.monpn.mn.pseudotime_curve1.pdf",width=10,height = length(marker)/4)
  print(p)
dev.off()

##### 20211226 mn0912 heatmap
marker <- read.table("20211226.MN.marker1.txt", header = F)$V1
cluster_mean.mn0912[c(1:5),c(1:5)]
data <- t(cluster_mean.mn0912)
data <- data[rownames(data) %in% marker,c(5,3,4,2,1)]
head(data)

pdf("20211226.MN.marker1.pdf",height = length(marker)/7, width = 3)
pheatmap::pheatmap(data, cluster_cols = F, scale = "row",
                   cellwidth = 10, cellheight = 10)
dev.off()

##### MN heatmap
alldata <- read.table("20211226.MN.heatmap.txt", header = F, sep = "\t")
colnames(alldata) <- c("gene","group")

for (gp in unique(alldata$group)){
  marker <- subset(alldata, group == gp,select=gene)$gene
  data <- t(cluster_mean.mn0912)
  data <- data[rownames(data) %in% marker,c(5,3,4,2,1)]
  head(data)
  data <- data[match(marker,rownames(data)),]
  
  pdf(paste("20211226.MN.group.",gp,".pdf",sep = ""),height = 6, width = 3)
  pheatmap::pheatmap(data, cluster_cols = F, scale = "row", cluster_rows = F,
                     cellwidth = 10, cellheight = 10)
  dev.off()
}

##### 20211227 
marker <- unique(alldata$gene)
pt.matrix <- exprs(cds.mn.0912)[match(marker,rownames(rowData(cds.mn.0912))),]
#pt.matrix <- project.pgsplit@assays$RNA@scale.data[rownames(project.pgsplit@assays$RNA@scale.data) %in% marker$V1,]
pt.matrix <- as.data.frame(t(pt.matrix))
dim(pt.matrix)
head(pt.matrix)

pt.matrix$pseudotime <- pseudotime(cds.mn.0912)
pt.matrix$cluster <- colData(cds.mn.0912)$clusters
data <- reshape2::melt(pt.matrix, id=c("pseudotime", "cluster"))
data$cluster <- factor(data$cluster, levels = c(2,3,10,12,13))
head(data)

col <- monpn.col[c(2,3,10,12,13)]

p<-ggplot(data) + geom_point(aes(x = pseudotime, y = value, color = cluster),size = 0.5) +
  geom_smooth(aes(x=pseudotime, y=value, color = cluster)) + scale_color_manual(values = col) +
  theme_bw()+  facet_wrap(.~variable,ncol=4,scales = "free_y") +  scale_y_log10() +
  theme(panel.grid = element_blank(), strip.background = element_blank()) +
  xlab("Pesudotime") + ylab("Expression/log10")
p

pdf("20211226.monpn.mn.pseudotime_curve2.pdf",width=10,height = length(marker)/4)
  print(p)
dev.off()

##### 20211227 INsplit dotplot
library(magrittr)
library(Seurat)

source("/data2/chenzixi/20210419_WJ_SC/seurat_fetch_dotplot_data.r")
marker <- read.table("20211227.INsplit.marker1.txt", header = T, sep = "\t")

# use RenameIdents and SetIdents
project.insplit@active.ident <- project.insplit@meta.data$RNA_snn_res.3
idents<-as.factor(as.numeric(as.vector(project.insplit@active.ident))+1)
names(idents) <- rownames(project.insplit@meta.data)
project.insplit@active.ident <- idents
project.insplit@active.ident

cells.use <- WhichCells(project.insplit, idents = c(17,4,24))
project.insplit <- SetIdent(project.insplit, cells = cells.use, value = 17)

cells.use <- WhichCells(project.insplit, idents = c(8,25,10,7,9))
project.insplit <- SetIdent(project.insplit, cells = cells.use, value = 8)

cells.use <- WhichCells(project.insplit, idents = c(21,13))
project.insplit <- SetIdent(project.insplit, cells = cells.use, value = 21)

cells.use <- WhichCells(project.insplit, idents = c(2,3,14))
project.insplit <- SetIdent(project.insplit, cells = cells.use, value = 2)

cells.use <- WhichCells(project.insplit, idents = c(6,15))
project.insplit <- SetIdent(project.insplit, cells = cells.use, value = 6)

cells.use <- WhichCells(project.insplit, idents = c(22,26))
project.insplit <- SetIdent(project.insplit, cells = cells.use, value = 22)

cells.use <- WhichCells(project.insplit, idents = c(18,27))
project.insplit <- SetIdent(project.insplit, cells = cells.use, value = 18)

cells.use <- WhichCells(project.insplit, idents = c(1,5,11,23))
project.insplit <- SetIdent(project.insplit, cells = cells.use, value = 1)

project.insplit@active.ident
plotdata <- fetch_dotplot_data(project.insplit, features = as.vector(unique(marker$gene)), 
                               dot.scale = 6) 

head(plotdata)

plotdata$id <- factor(plotdata$id, levels = rev(c(20,16,17,8,21,2,6,22,18,1,28,12,19)))
plotdata$na.avg.exp <- plotdata$avg.exp
plotdata$type[plotdata$id %in% c(20,16,17,8,21,2,6)] <- "dorsal"
plotdata$type[plotdata$id %in% c(18,1,28,12,19)] <- "ventral"
plotdata$type[plotdata$id %in% c(22)] <- "mix"


p <- ggplot(plotdata) + geom_point(aes(features.plot, id, fill= type, size = avg.exp),shape = 21, color = "grey") + theme_bw() +
  scale_fill_manual(values = col[c(5,2,3)], name = "Celltype") + RotatedAxis() + scale_size_area(name = "Average expression")+
  xlab(NULL) + ylab("Subcluster") + theme(axis.text = element_text(color = "black"))
p

ggsave("20211227.INtype.DotPlot.pdf",p,width=30,height=4)

##### 20211228 gene1
source("/data2/chenzixi/20210419_WJ_SC/seurat_fetch_dotplot_data.r")
marker <- read.table("20211228.IN.gene1.txt", header = F, sep = "\t")

# use RenameIdents and SetIdents
project.insplit@active.ident <- project.insplit@meta.data$RNA_snn_res.3
idents<-as.factor(as.numeric(as.vector(project.insplit@active.ident))+1)
names(idents) <- rownames(project.insplit@meta.data)
project.insplit@active.ident <- idents
project.insplit@active.ident

cells.use <- WhichCells(project.insplit, idents = c(17,4,24))
project.insplit <- SetIdent(project.insplit, cells = cells.use, value = 17)

cells.use <- WhichCells(project.insplit, idents = c(8,25,10,7,9))
project.insplit <- SetIdent(project.insplit, cells = cells.use, value = 8)

cells.use <- WhichCells(project.insplit, idents = c(21,13))
project.insplit <- SetIdent(project.insplit, cells = cells.use, value = 21)

cells.use <- WhichCells(project.insplit, idents = c(2,3,14))
project.insplit <- SetIdent(project.insplit, cells = cells.use, value = 2)

cells.use <- WhichCells(project.insplit, idents = c(6,15))
project.insplit <- SetIdent(project.insplit, cells = cells.use, value = 6)

cells.use <- WhichCells(project.insplit, idents = c(22,26))
project.insplit <- SetIdent(project.insplit, cells = cells.use, value = 22)

cells.use <- WhichCells(project.insplit, idents = c(18,27))
project.insplit <- SetIdent(project.insplit, cells = cells.use, value = 18)

cells.use <- WhichCells(project.insplit, idents = c(1,5,11,23))
project.insplit <- SetIdent(project.insplit, cells = cells.use, value = 1)

project.insplit@active.ident
plotdata <- fetch_dotplot_data(project.insplit, features = as.vector(unique(marker$V1)), 
                               dot.scale = 6) 

head(plotdata)

plotdata$id <- factor(plotdata$id, levels = rev(c(20,16,17,8,21,2,6,22,18,1,28,12,19)))
plotdata$na.avg.exp <- plotdata$avg.exp
plotdata$type[plotdata$id %in% c(20,16,17,8,21,2,6)] <- "dorsal"
plotdata$type[plotdata$id %in% c(18,1,28,12,19)] <- "ventral"
plotdata$type[plotdata$id %in% c(22)] <- "mix"


p <- ggplot(plotdata) + geom_point(aes(features.plot, id, fill= type, size = avg.exp),shape = 21, color = "grey") + theme_bw() +
  scale_fill_manual(values = col[c(5,2,3)], name = "Celltype") + RotatedAxis() + scale_size_area(name = "Average expression")+
  xlab(NULL) + ylab("Subcluster") + theme(axis.text = element_text(color = "black"))
p

ggsave("20211228.INtype.gene1.DotPlot.pdf",p,width=15,height=4)

##### 20211228 TF1
source("/data2/chenzixi/20210419_WJ_SC/seurat_fetch_dotplot_data.r")
marker <- read.table("20211228.IN.TF1.txt", header = F, sep = "\t")

# use RenameIdents and SetIdents
project.insplit@active.ident <- project.insplit@meta.data$RNA_snn_res.3
idents<-as.factor(as.numeric(as.vector(project.insplit@active.ident))+1)
names(idents) <- rownames(project.insplit@meta.data)
project.insplit@active.ident <- idents
project.insplit@active.ident

cells.use <- WhichCells(project.insplit, idents = c(17,4,24))
project.insplit <- SetIdent(project.insplit, cells = cells.use, value = 17)

cells.use <- WhichCells(project.insplit, idents = c(8,25,10,7,9))
project.insplit <- SetIdent(project.insplit, cells = cells.use, value = 8)

cells.use <- WhichCells(project.insplit, idents = c(21,13))
project.insplit <- SetIdent(project.insplit, cells = cells.use, value = 21)

cells.use <- WhichCells(project.insplit, idents = c(2,3,14))
project.insplit <- SetIdent(project.insplit, cells = cells.use, value = 2)

cells.use <- WhichCells(project.insplit, idents = c(6,15))
project.insplit <- SetIdent(project.insplit, cells = cells.use, value = 6)

cells.use <- WhichCells(project.insplit, idents = c(22,26))
project.insplit <- SetIdent(project.insplit, cells = cells.use, value = 22)

cells.use <- WhichCells(project.insplit, idents = c(18,27))
project.insplit <- SetIdent(project.insplit, cells = cells.use, value = 18)

cells.use <- WhichCells(project.insplit, idents = c(1,5,11,23))
project.insplit <- SetIdent(project.insplit, cells = cells.use, value = 1)

project.insplit@active.ident
plotdata <- fetch_dotplot_data(project.insplit, features = as.vector(unique(marker$V1)), 
                               dot.scale = 6) 

head(plotdata)

plotdata$id <- factor(plotdata$id, levels = rev(c(20,16,17,8,21,2,6,22,18,1,28,12,19)))
plotdata$na.avg.exp <- plotdata$avg.exp
plotdata$type[plotdata$id %in% c(20,16,17,8,21,2,6)] <- "dorsal"
plotdata$type[plotdata$id %in% c(18,1,28,12,19)] <- "ventral"
plotdata$type[plotdata$id %in% c(22)] <- "mix"


p <- ggplot(plotdata) + geom_point(aes(features.plot, id, fill= type, size = avg.exp),shape = 21, color = "grey") + theme_bw() +
  scale_fill_manual(values = col[c(5,2,3)], name = "Celltype") + RotatedAxis() + scale_size_area(name = "Average expression")+
  xlab(NULL) + ylab("Subcluster") + theme(axis.text = element_text(color = "black"))
p

ggsave("20211228.INtype.TF1.DotPlot.pdf",p,width=15,height=4)

##### 20211228 gene2
source("/data2/chenzixi/20210419_WJ_SC/seurat_fetch_dotplot_data.r")
marker <- read.table("20211228.IN.gene2.txt", header = F, sep = "\t")

# use RenameIdents and SetIdents
project.insplit@active.ident <- project.insplit@meta.data$RNA_snn_res.3
idents<-as.factor(as.numeric(as.vector(project.insplit@active.ident))+1)
names(idents) <- rownames(project.insplit@meta.data)
project.insplit@active.ident <- idents
project.insplit@active.ident

cells.use <- WhichCells(project.insplit, idents = c(17,4,24))
project.insplit <- SetIdent(project.insplit, cells = cells.use, value = 17)

cells.use <- WhichCells(project.insplit, idents = c(8,25,10,7,9))
project.insplit <- SetIdent(project.insplit, cells = cells.use, value = 8)

cells.use <- WhichCells(project.insplit, idents = c(21,13))
project.insplit <- SetIdent(project.insplit, cells = cells.use, value = 21)

cells.use <- WhichCells(project.insplit, idents = c(2,3,14))
project.insplit <- SetIdent(project.insplit, cells = cells.use, value = 2)

cells.use <- WhichCells(project.insplit, idents = c(6,15))
project.insplit <- SetIdent(project.insplit, cells = cells.use, value = 6)

cells.use <- WhichCells(project.insplit, idents = c(22,26))
project.insplit <- SetIdent(project.insplit, cells = cells.use, value = 22)

cells.use <- WhichCells(project.insplit, idents = c(18,27))
project.insplit <- SetIdent(project.insplit, cells = cells.use, value = 18)

cells.use <- WhichCells(project.insplit, idents = c(1,5,11,23))
project.insplit <- SetIdent(project.insplit, cells = cells.use, value = 1)

project.insplit@active.ident
plotdata <- fetch_dotplot_data(project.insplit, features = as.vector(unique(marker$V1)), 
                               dot.scale = 6) 

head(plotdata)

plotdata$id <- factor(plotdata$id, levels = rev(c(20,16,17,8,21,2,6,22,18,1,28,12,19)))
plotdata$na.avg.exp <- plotdata$avg.exp
plotdata$type[plotdata$id %in% c(20,16,17,8,21,2,6)] <- "dorsal"
plotdata$type[plotdata$id %in% c(18,1,28,12,19)] <- "ventral"
plotdata$type[plotdata$id %in% c(22)] <- "mix"


p <- ggplot(plotdata) + geom_point(aes(features.plot, id, fill= type, size = avg.exp),shape = 21, color = "grey") + theme_bw() +
  scale_fill_manual(values = col[c(5,2,3)], name = "Celltype") + RotatedAxis() + scale_size_area(name = "Average expression")+
  xlab(NULL) + ylab("Subcluster") + theme(axis.text = element_text(color = "black"))
p

ggsave("20211228.INtype.gene2.DotPlot.pdf",p,width=10,height=4)

##### 20211229 MN 
marker <- unique(read.table("20211229.MN.marker1.txt", header = F)$V1)

data <- t(cluster_mean.mn0912)
data <- data[rownames(data) %in% marker,c(5,3,4,2,1)]
head(data)
data <- data[match(marker[marker %in% rownames(data)],rownames(data)),]

pdf(paste("20211229.MN.group.gene-MMC.pdf",sep = ""),height = 6, width = 3)
pheatmap::pheatmap(data, cluster_cols = F, scale = "row", cluster_rows = F,
                   cellwidth = 10, cellheight = 10)
dev.off()

marker <- unique(read.table("20211229.MN.marker2.txt", header = F)$V1)

data <- t(cluster_mean.mn0912)
data <- data[rownames(data) %in% marker,c(5,3,4,2,1)]
head(data)
data <- data[match(marker[marker %in% rownames(data)],rownames(data)),]

pdf(paste("20211229.MN.group.gene-MMC2.pdf",sep = ""),height = 6, width = 3)
pheatmap::pheatmap(data, cluster_cols = F, scale = "row", cluster_rows = F,
                   cellwidth = 10, cellheight = 10)
dev.off()

##### 20211229 MN curve1
marker <- read.table("20211229.MN.curve.marker1.txt", header = F, sep = "\t")$V1

pt.matrix <- exprs(cds.mn.0912)[match(marker,rownames(rowData(cds.mn.0912))),]
#pt.matrix <- project.pgsplit@assays$RNA@scale.data[rownames(project.pgsplit@assays$RNA@scale.data) %in% marker$V1,]
pt.matrix <- as.data.frame(t(pt.matrix))
dim(pt.matrix)
head(pt.matrix)

pt.matrix$pseudotime <- pseudotime(cds.mn.0912)
pt.matrix$cluster <- colData(cds.mn.0912)$clusters
data <- reshape2::melt(pt.matrix, id=c("pseudotime", "cluster"))
data$cluster <- factor(data$cluster, levels = c(2,3,10,12,13))
head(data)

col <- monpn.col[c(2,3,10,12,13)]

p<-ggplot(data) + geom_point(aes(x = pseudotime, y = value, color = cluster),size = 0.5) +
  geom_smooth(aes(x=pseudotime, y=value),color="grey30") + scale_color_manual(values = col) +
  theme_bw()+  facet_wrap(.~variable,ncol=4,scales = "free_y") +  scale_y_log10() +
  theme(panel.grid = element_blank(), strip.background = element_blank()) +
  xlab("Pesudotime") + ylab("Expression/log10")
p

pdf("20211229.monpn.mn.pseudotime_curve1.pdf",width=10,height = length(marker)/4)
  print(p)
dev.off()

p<-ggplot(data) + geom_point(aes(x = pseudotime, y = value, color = cluster),size = 0.5) +
  geom_smooth(aes(x=pseudotime, y=value, color = cluster)) + scale_color_manual(values = col) +
  theme_bw()+  facet_wrap(.~variable,ncol=4,scales = "free_y") +  scale_y_log10() +
  theme(panel.grid = element_blank(), strip.background = element_blank()) +
  xlab("Pesudotime") + ylab("Expression/log10")
p

##### 20211230
#cluster_mean.monpn
#expr <- as.data.frame(t(project.pgsplit@assays$RNA@scale.data))
#expr[1:5,1:5]
expr <- as.data.frame(t(exprs(cds.monpn)))
expr[1:5,1:5]
expr$cluster <- clusters(cds.monpn)[match(names(clusters(cds.monpn)),rownames(expr))]
attach(expr)
cluster_mean.cdsmonpn <-  as.data.frame(aggregate(expr, list(cluster), mean))
rownames(cluster_mean.cdsmonpn) <- cluster_mean.cdsmonpn$Group.1
cluster_mean.cdsmonpn[1:5,1:5]
cluster_mean.cdsmonpn <- cluster_mean.cdsmonpn[,-1]
cluster_mean.cdsmonpn[1:5,1:5]

marker <- read.table("20211230.marker1.txt", header = F)$V1
cluster_mean.cdsmonpn[c(1:5),c(1:5)]
data <- t(cluster_mean.cdsmonpn)
data[c(1:5),c(1:5)]
data <- data[rownames(data) %in% marker,c(6,4,10,11,1,9,5,7,8,3,13,2,12)]
head(data)

pdf("20211230.monpn.marker1.heatmap.pdf",height = 8, width = 4)
pheatmap(data,cluster_cols = F, scale = "row", cluster_rows = F,
         cellwidth = 10, cellheight = 10)
dev.off()

colData(cds.monpn)$clusters <- factor(colData(cds.monpn)$clusters, levels = c(6,4,10,11,1,9,5,7,8,3,13,2,12))

pdf("20211230.monpn.dotplot.pdf",height = 20, width = 8)
  plot_genes_by_group(cds.monpn,rev(marker),ordering_type = "none",group_cells_by = "clusters")
dev.off()

##### 20211230 curve
marker <- read.table("20211230.marker1.txt", header = F)$V1

cds.monpn <- order_cells(cds.monpn, root_pr_nodes= c("Y_25","Y_49"))
plot_cells(cds.monpn, trajectory_graph_segment_size = 0.25,
                   color_cells_by = "pseudotime",
                   label_cell_groups=FALSE,
                   label_leaves=FALSE,
                   label_branch_points=FALSE,
                   graph_label_size=1.5)+ggtitle("Root Y_25 Y_49")

pt.matrix <- exprs(cds.monpn)[match(marker,rownames(rowData(cds.monpn))),]
#pt.matrix <- project.pgsplit@assays$RNA@scale.data[rownames(project.pgsplit@assays$RNA@scale.data) %in% marker$V1,]
pt.matrix <- as.data.frame(t(pt.matrix))
dim(pt.matrix)
head(pt.matrix)

pt.matrix$pseudotime <- pseudotime(cds.monpn)
pt.matrix$cluster <- colData(cds.monpn)$clusters
data <- reshape2::melt(pt.matrix, id=c("pseudotime", "cluster"))
data$cluster <- factor(data$cluster, levels = c(1:13))
head(data)


# mn0912
data.plot <- data[data$cluster %in% c(2,3,10,12,13),]

col <- monpn.col[c(2,3,10,12,13)]

p<-ggplot(data.plot) + geom_point(aes(x = pseudotime, y = value, color = cluster),size = 0.5) +
  geom_smooth(aes(x=pseudotime, y=value),color="grey30") + scale_color_manual(values = col) +
  theme_bw()+  facet_wrap(.~variable,ncol=4,scales = "free_y") +  scale_y_log10() +
  theme(panel.grid = element_blank(), strip.background = element_blank()) +
  xlab("Pesudotime") + ylab("Expression/log10")
p

pdf("20211230.monpn.mn.pseudotime_curve1.pdf",width=10,height = length(marker)/4)
  print(p)
dev.off()

### in0912
data.plot <- data[data$cluster %in% c(1,5,7,8,9,11),]
col <- monpn.col[c(1,5,7,8,9,11)]

p<-ggplot(data.plot) + geom_point(aes(x = pseudotime, y = value, color = cluster),size = 0.5) +
  geom_smooth(aes(x=pseudotime, y=value),color="grey30") + scale_color_manual(values = col) +
  theme_bw()+  facet_wrap(.~variable,ncol=4,scales = "free_y") +  scale_y_log10() +
  theme(panel.grid = element_blank(), strip.background = element_blank()) +
  xlab("Pesudotime") + ylab("Expression/log10")
#p

pdf("20211230.monpn.in.pseudotime_curve1.pdf",width=10,height = length(marker)/4)
  print(p)
dev.off()

### pg0912
data.plot <- data[data$cluster %in% c(4,6,10),]
col <- monpn.col[c(4,6,10)]

p<-ggplot(data.plot) + geom_point(aes(x = pseudotime, y = value, color = cluster),size = 0.5) +
  geom_smooth(aes(x=pseudotime, y=value),color="grey30") + scale_color_manual(values = col) +
  theme_bw()+  facet_wrap(.~variable,ncol=4,scales = "free_y") +  scale_y_log10() +
  theme(panel.grid = element_blank(), strip.background = element_blank()) +
  xlab("Pesudotime") + ylab("Expression/log10")
p

pdf("20211230.monpn.pg.pseudotime_curve1.pdf",width=10,height = length(marker)/4)
  print(p)
dev.off()

##### 20220115 D vs V GSEA
library(magrittr)
library(tidyverse)
library(DOSE)
library(enrichplot)

GO_DATA_NRVC <- clusterProfiler:::get_GO_data("org.Hs.eg.db", "ALL", "SYMBOL")
PATHID2NAME <- as.data.frame(GO_DATA_NRVC$PATHID2NAME)
PATHID2NAME$GO <- rownames(PATHID2NAME)
colnames(PATHID2NAME) <- c("Description","GOTerm")

GO_BP <- GO_DATA_NRVC$PATHID2NAME[GO_DATA_NRVC$GO2ONT == "BP"]
#Wnt_NRVCGO<-names(GO_DATA_NRVC$PATHID2NAME[grep("Wnt", GO_DATA_NRVC$PATHID2NAME)])
genelist <- plyr::ldply(GO_DATA_NRVC$PATHID2EXTID, data.frame)
genelist.bp <- genelist[genelist[,1] %in% names(GO_BP),]

#genelist.sp <- genelist[genelist[,1] %in% c("GO:0030509","GO:0007224","GO:0060071"),]

# # use padj
# dor_vs_ven$gene <- rownames(dor_vs_ven)
# dv.genes <- dor_vs_ven %>% arrange(p_val_adj) %>% dplyr::select(gene, p_val_adj)
# dv.ranks <- rev(deframe(dv.genes[dv.genes$p_val_adj<0.1,]))
# dv.ranks
# 
# gsea.dv <- clusterProfiler::GSEA(dv.ranks, TERM2GENE = genelist.bp, verbose=FALSE, pvalueCutoff = 1, TERM2NAME = PATHID2NAME[,c(2,1)],
#                                  scoreType = "pos", eps = 0)
# write.table(gsea.dv,"20220115.dvgsea.xls", sep = "\t", quote = F, col.names = T, row.names = F)

## try log2FC
dv.genes <- dor_vs_ven %>% arrange(desc(avg_log2FC)) %>% dplyr::select(gene, avg_log2FC)
dv.ranks <- deframe(dv.genes)
dv.ranks

gsea.dv <- clusterProfiler::GSEA(dv.ranks, TERM2GENE = genelist.bp, TERM2NAME = PATHID2NAME[,c(2,1)],
                                 verbose=FALSE, pvalueCutoff = 1, eps = 0)
write.table(gsea.dv,"20220115.dvgsea.fc.xls", sep = "\t", quote = F, col.names = T, row.names = F)


gsea.dv.sp <- clusterProfiler::GSEA(dv.ranks, TERM2GENE = genelist.sp, TERM2NAME = PATHID2NAME[,c(2,1)],
                                 verbose=FALSE, pvalueCutoff = 1, eps = 0)

gseaplot2(gsea.dv.sp, gsea.dv.sp@result$Description ,color = "blue", subplots = 1:2)
       
#tyr gseGO   
gse_obj <- gseGO(geneList     = dv.ranks,
                 keyType = "SYMBOL",
                 OrgDb        = org.Hs.eg.db,
                 ont          = "BP",
                 exponent     = 1,
                 nPerm        = 10000,
                 minGSSize    = 5,
                 maxGSSize    = 500,
                 pvalueCutoff = 1,
                 pAdjustMethod= "BH",
                 verbose      = TRUE,
                 seed         = TRUE)
write.table(gse_obj,"20220115.dvgsea.fc.gse_obj.xls", sep = "\t", quote = F, col.names = T, row.names = F)

#"GO:0030509","GO:0007224","GO:0060071"
gseaplot2(gse_obj, c("GO:0030509","GO:0007224","GO:0060071"), subplots = 1:2)
ggsave("20220115.dv.gsea.pdf")

##### 20220116 try use significant genes for gsea
#p0.05
dor_vs_ven.p0.05 <- dor_vs_ven[dor_vs_ven$p_val < 0.05,]
dv.genes <- dor_vs_ven.p0.05 %>% arrange(desc(avg_log2FC)) %>% dplyr::select(gene, avg_log2FC)
dv.ranks <- deframe(dv.genes)
dv.ranks

gse_obj.sig <- gseGO(geneList     = dv.ranks,
                 keyType = "SYMBOL",
                 OrgDb        = org.Hs.eg.db,
                 ont          = "BP",
                 pvalueCutoff = 1,
                 eps = 0,
                 pAdjustMethod= "BH",
                 verbose      = TRUE,
                 seed         = TRUE)
write.table(gse_obj.sig,"20220115.dvgsea.fc.gse_obj.p0.05.xls", sep = "\t", quote = F, col.names = T, row.names = F)

#padj 0.05
dor_vs_ven.padj0.05 <- dor_vs_ven[dor_vs_ven$p_val_adj < 0.05,]
dv.genes <- dor_vs_ven.padj0.05 %>% arrange(desc(avg_log2FC)) %>% dplyr::select(gene, avg_log2FC)
dv.ranks <- deframe(dv.genes)
dv.ranks

gse_obj.sig <- gseGO(geneList     = dv.ranks,
                     keyType = "SYMBOL",
                     OrgDb        = org.Hs.eg.db,
                     ont          = "BP",
                     pvalueCutoff = 1,
                     eps = 0,
                     pAdjustMethod= "BH",
                     verbose      = TRUE,
                     seed         = TRUE)
write.table(gse_obj.sig,"20220115.dvgsea.fc.gse_obj.padj0.05.xls", sep = "\t", quote = F, col.names = T, row.names = F)

#all
dor_vs_ven.all <- dor_vs_ven
dv.genes <- dor_vs_ven.all %>% arrange(desc(avg_log2FC)) %>% dplyr::select(gene, avg_log2FC)
dv.ranks <- deframe(dv.genes)
dv.ranks

gse_obj.sig <- gseGO(geneList     = dv.ranks,
                     keyType = "SYMBOL",
                     OrgDb        = org.Hs.eg.db,
                     ont          = "BP",
                     pvalueCutoff = 1,
                     eps = 0,
                     pAdjustMethod= "BH",
                     verbose      = TRUE,
                     seed         = TRUE)
write.table(gse_obj.sig,"20220115.dvgsea.fc.gse_obj.all.xls", sep = "\t", quote = F, col.names = T, row.names = F)

gse_obj.sig.all <- gse_obj.sig
saveRDS(gse_obj.sig.all,"gse_obj.sig.all.RDS")

# Wnt
wnt <- c("GO:2000050",
         "GO:0007223",
         "GO:0090090",
         "GO:0016055",
         "GO:0198738",
         "GO:0030178",
         "GO:0030177",
         "GO:0035567",
         "GO:0060070",
         "GO:0030111",
         "GO:0090263",
         "GO:0060828",
         "GO:0060071")

gseaplot2(gse_obj.sig, wnt[wnt %in% gse_obj.sig@result$ID], subplots = 1:2)
# BMP
bmp <- c("GO:0030509",
         "GO:0071772",
         "GO:0071773",
         "GO:0030510",
         "GO:0030513",
         "GO:0030514")

gseaplot2(gse_obj.sig, bmp[bmp %in% gse_obj.sig@result$ID], subplots = 1:2)

# Shh
shh <- c("GO:0008589",
         "GO:0045880",
         "GO:0007224",
         "GO:0045879")

pdf("20220117.gsea.PG.pdf",width = 6,height = 4)
  gseaplot2(gse_obj.sig, c("GO:0030509","GO:0007224"), subplots = 1:2, color = c("red","blue")) 
  gseaplot2(gse_obj.sig, c("GO:0030509","GO:0008589"), subplots = 1:2, color = c("red","blue")) 
dev.off()
##### 20220117 PGsplit SCENIC AUC
library(monocle3)
library(AUCell)
load("/data2/chenzixi/20210419_WJ_SC/20210831_s4hV6Sub_SCENIC/PG0912/20210912_s4hSubV6_PG_SCENIC.Final.RData")
dim(aucell_regulonAUC)
dim(cds.pgsplit)

plot_cells(cds.pgsplit)
umap.coord <- cds.pgsplit@int_colData@listData$reducedDims$UMAP
plot_cells()
exprMat[1:5,1:5]

marker <- read.table("20220117.PGsplit.pIN.txt", header = F)$V1

pdf("20220117.PGsplit.pIN.Regulon.AUC.pdf",height = 6,width = 8)
AUCell::AUCell_plotTSNE(umap.coord,exprMat, 
                        aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))
                                          [marker],], 
                        plots="AUC",cex = 0.75)
dev.off()

pdf("20220117.PGsplit.pIN.Regulon.Expression.pdf",height = 6,width = 8)
AUCell::AUCell_plotTSNE(umap.coord,exprMat, 
                        aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))
                                          [marker],], 
                        plots="expression",cex = 0.75)
dev.off()


marker <- read.table("20220117.PGsplit.pMN.txt", header = F)$V1

pdf("20220117.PGsplit.pMN.Regulon.AUC.pdf",height = 6,width = 8)
AUCell::AUCell_plotTSNE(umap.coord,exprMat, 
                        aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))
                                          [marker],], 
                        plots="AUC",cex = 0.75)
dev.off()

pdf("20220117.PGsplit.pMN.Regulon.Expression.pdf",height = 6,width = 8)
AUCell::AUCell_plotTSNE(umap.coord,exprMat, 
                        aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))
                                          [marker],], 
                        plots="expression",cex = 0.75)
dev.off()

### all heatmap
marker1 <- read.table("20220117.PGsplit.pIN.txt", header = F)$V1
marker2 <- read.table("20220117.PGsplit.pMN.txt", header = F)$V1
marker <- c(marker1,marker2)

cellInfo <- data.frame(cellid=names(project.pgsplit@active.ident),CellType=project.pgsplit@active.ident)
#cellInfo <- data.frame(cellid=names(cds.pgsplit@clusters$UMAP$clusters),CellType=cds.pgsplit@clusters$UMAP$clusters)
cellInfo <- data.frame(cellInfo)
cellInfo <- cellInfo[rownames(cellInfo) %in% colnames(exprMat),]
head(cellInfo)

cellInfo <- data.frame(CellType=cellInfo$CellType, row.names = cellInfo$cellid)
head(cellInfo)
cellInfo$CellType <- factor(cellInfo$CellType, levels = sort(unique(cellInfo$CellType)))

#
#regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC))[marker],]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$CellType),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
#regulonActivity_byCellType <- regulonActivity_byCellType[,c(3,5,2,1,4)]

regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
#regulonActivity_byCellType_Scaled <- regulonActivity_byCellType_Scaled[,c(3,5,2,1,4)]
#regulonActivity_byCellType <- regulonActivity_byCellType[,c(3,5,2,1,4)]

ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled[,c(9,7,8,3,6,10,2,4,5,1)], name="Regulon activity\nz-score", 
                        cluster_rows = F, cluster_columns = F )
p <- ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled[,c(6,1,8,7,4,3,5,2)], name="Regulon activity\nz-score", 
                        cluster_rows = F, cluster_columns = F )


pdf("20220117.PGsplit.SCENIC.heatmap.pdf",height = 10,width = 8)
print(p)
dev.off()
## ggplot2 
library(ggplot2)
library(viridis)

regulons <- rownames(p@matrix)
regulon_auc <- aucell_regulonAUC[rownames(aucell_regulonAUC) %in% regulons,]
dim(regulon_auc)
data <- as.data.frame(t(as.data.frame(regulon_auc@assays@data$AUC)))
data[c(1:5),c(1:5)]

umap.coord <- cds.pgsplit@int_colData@listData$reducedDims$UMAP
umap.coord[,1]
data$umap1 <- as.vector(umap.coord[,1])
data$umap2 <- as.vector(umap.coord[,2])

test <- data[,c(1,67,68)]

pdf("20220117.PG.auc.pdf")
for (i in c(1:length(regulons))){
  test <- data[,c(i,67,68)]
#test$`PAX6_extended (33g)` <- scale(test$`PAX6_extended (33g)`)
colnames(test) <- c("AUCell","umap1","umap2")
#plotdata <- reshape2::melt(data,id=c("umap1", "umap2"))
p<-ggplot(test) + geom_point(aes(umap1,umap2, color=AUCell),size=0.25,alpha=0.75) + theme_classic() + coord_fixed() +
  scale_color_gradient2(low="grey95",mid = "lightgrey", high="red",midpoint = (min(test$AUCell)+max(test$AUCell))/5) +
  ggtitle(regulons[i]) + xlab("UMAP_1") + ylab("UMAP_2")
print(p)
}
dev.off()

ggplot(test) + geom_point(aes(umap1,umap2, color=AUCell),size=0.25,alpha=0.75) + theme_classic() + coord_fixed() +
  scale_color_gradient2(low="grey95",mid = "lightgrey", high="red",midpoint = (min(test$AUCell)+max(test$AUCell))/5) +
  ggtitle(regulons[i])

##### 20220117 MN scenic heatmap
library(ComplexHeatmap)
library(pheatmap)
load("/data2/chenzixi/20210419_WJ_SC/20210831_s4hV6Sub_SCENIC/MN0912/20210912_s4hSubV6_MN_SCENIC.Final.RData")
marker <- read.table("20220117.MN0912.scenic.txt", header = F)$V1

cellInfo <- data.frame(cellid=names(cds.mn.0912@clusters$UMAP$clusters),CellType=cds.mn.0912@clusters$UMAP$clusters)
cellInfo <- data.frame(cellInfo)
cellInfo <- cellInfo[rownames(cellInfo) %in% colnames(exprMat),]
head(cellInfo)

cellInfo <- data.frame(CellType=cellInfo$CellType, row.names = cellInfo$cellid)
head(cellInfo)
cellInfo$CellType <- factor(cellInfo$CellType, levels = sort(unique(cellInfo$CellType)))

#
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC))[marker],]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$CellType),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
#regulonActivity_byCellType <- regulonActivity_byCellType[,c(3,5,2,1,4)]

regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
regulonActivity_byCellType_Scaled <- regulonActivity_byCellType_Scaled[,c(3,5,2,1,4)]
regulonActivity_byCellType <- regulonActivity_byCellType[,c(3,5,2,1,4)]

pdf("20220117.MN0912.scenic.regulon.pdf",width = 6, height = 4)
ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity\nz-score", 
                        cluster_rows = F, cluster_columns = F )
dev.off()

##### 20220117  MN0912 show_trajectory_graph
monpn.col <- c(
  "#6ab662",
  "#377EB8",
  "#cab3d6",
  "#2ea37e",
  "#ebc4b5",
  "#e11517",
  "#efd092",
  "#c47e70",
  "#b8e294",
  "#774ea3",
  "#ffd92f",
  "#b5d7e8",
  "#685f7e"
)

p<- monocle3::plot_cells(cds.mn.0912, 
                     label_cell_groups = F, label_branch_points = F,
                     label_principal_points = F, label_roots = F, 
                     label_leaves = F,show_trajectory_graph = T) +
  scale_color_manual(values = monpn.col[c(2,3,10,12,13)]) + coord_fixed()

ggsave("20220117_s4hV6Sub_Monocle3_MN0912.pdf",p,width = 7.5,height = 4.5)

##### 20220117 MN0912 curve and rainbow-heatmap
## curve
# 3,10,13
marker <- read.table("20220117.MN0912.curve.txt", header = F)$V1

pt.matrix <- exprs(cds.mn.0912)[match(marker,rownames(rowData(cds.mn.0912))),]
#pt.matrix <- project.pgsplit@assays$RNA@scale.data[rownames(project.pgsplit@assays$RNA@scale.data) %in% marker$V1,]
pt.matrix <- as.data.frame(t(pt.matrix))
dim(pt.matrix)
head(pt.matrix)

pt.matrix$pseudotime <- pseudotime(cds.mn.0912)
pt.matrix$cluster <- colData(cds.mn.0912)$cluster 

head(pt.matrix)

#data <- pt.matrix[pt.matrix$variable %in% marker$V1,]
data <- pt.matrix[pt.matrix$cluster %in% c(10,3,13),]
data <- reshape2::melt(data, id=c("pseudotime", "cluster"))
data$celltype <- factor(data$cluster, levels = c(3,10,13))

p<-ggplot(data) + geom_point(aes(x=pseudotime, y=value, color = cluster),size = 0.1) +
  geom_smooth(aes(x=pseudotime, y=value),size = 0.5) + scale_color_manual(values = monpn.col[c(3,10,13)]) +
  theme_bw()+  facet_wrap(.~variable,ncol=4,scales = "free_y") +  
  theme(panel.grid = element_blank(), strip.background = element_blank()) +
  xlab("Pesudotime") + ylab("Expression")
p

pdf("20220117.MN0912.to13.pdf",width=10,height = length(marker)/4)
  print(p)
dev.off()

# 2,3,10,13
marker <- read.table("20220117.MN0912.curve.txt", header = F)$V1

pt.matrix <- exprs(cds.mn.0912)[match(marker,rownames(rowData(cds.mn.0912))),]
#pt.matrix <- project.pgsplit@assays$RNA@scale.data[rownames(project.pgsplit@assays$RNA@scale.data) %in% marker$V1,]
pt.matrix <- as.data.frame(t(pt.matrix))
dim(pt.matrix)
head(pt.matrix)

pt.matrix$pseudotime <- pseudotime(cds.mn.0912)
pt.matrix$cluster <- colData(cds.mn.0912)$cluster 

head(pt.matrix)

data <- pt.matrix[pt.matrix$cluster %in% c(10,2,3,12),]
data <- reshape2::melt(data, id=c("pseudotime", "cluster"))
data$celltype <- factor(data$cluster, levels = c(2,3,10,12))

p<-ggplot(data) + geom_point(aes(x=pseudotime, y=value, color = cluster),size = 0.1) +
  geom_smooth(aes(x=pseudotime, y=value),size = 0.5) + scale_color_manual(values = monpn.col[c(2,3,10,12)]) +
  theme_bw()+  facet_wrap(.~variable,ncol=4,scales = "free_y") +  
  theme(panel.grid = element_blank(), strip.background = element_blank()) +
  xlab("Pesudotime") + ylab("Expression")
p

pdf("20220117.MN0912.to12.pdf",width=10,height = length(marker)/4)
  print(p)
dev.off()

## heatmap
library(ComplexHeatmap)
library(circlize)

pseudotime <- pseudotime(cds.mn.0912)

# 10,3,13,2,12
df.all <- data.frame(Cluster = colData(cds.mn.0912)$cluster[match(names(pseudotime(cds.mn.0912)[order(pseudotime(cds.mn.0912))]),
                                                                    names(colData(cds.mn.0912)$cluster))], 
                 Pseudotime = pseudotime(cds.mn.0912)[order(pseudotime(cds.mn.0912))])


ClusterCol <- monpn.col[c(10,3,13,2,12)]
#names(ClusterCol) <- sort(unique(colData(cds.mn.0912)$Seurat_res0.2))
names(ClusterCol) <- c(10,3,13,2,12)
ClusterCol

max_pseudo <- ceiling(max(df.all$Pseudotime))
# 
# ha <- HeatmapAnnotation(df = df,col = list(Cluster = ClusterCol, 
#                                            Pseudotime = colorRamp2(seq(from=0,to=max_pseudo,length=11),rev(brewer.pal(11, "RdBu")))),
#                         annotation_name_side = "left")


df1 <- subset(df.all, Cluster %in% c(10,3,13,2,12), select=Cluster )
df2 <- subset(df.all, Cluster %in% c(10,3,13,2,12), select=Pseudotime )

ha1 <- HeatmapAnnotation(df = df1, col = list(Cluster = ClusterCol), 
                         annotation_name_side = "left")
ha2 <- HeatmapAnnotation(df = df2, col = list(Pseudotime = colorRamp2(seq(from=0,to=max_pseudo,length=11),rev(brewer.pal(11, "RdBu")))), 
                         annotation_name_side = "left")


marker <- read.table("20220117.MN0912.curve.txt", header = F)$V1

genes <- marker[marker %in% rownames(rowData(cds.mn.0912))]
  
#
pt.matrix <- exprs(cds.mn.0912)[match(genes,rownames(rowData(cds.mn.0912))),
                                match(rownames(df.all),colnames(exprs(cds.mn.0912)))]
#                                order(pseudotime(cds.mn.0912)[names(pseudotime(cds.mn.0912)) %in% rownames(df)])]
pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes
  
#Ward.D2 Hierarchical Clustering
pdf(paste("20220117_MN0912.all.rainbowheatmap.pdf",sep=""),width=12,height = length(genes)/3)
  hthc <- Heatmap(
    as.matrix(pt.matrix),
    name                         = "z-score",
    col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
    show_row_names               = TRUE,
    show_column_names            = FALSE,
    row_names_gp                 = gpar(fontsize = 6),
    clustering_method_rows = "ward.D2",
    clustering_method_columns = "ward.D2",
    row_title_rot                = 0,
    cluster_rows                 = TRUE,
    cluster_row_slices           = FALSE,
    cluster_columns              = FALSE,
    top_annotation               = ha1,
    bottom_annotation            = ha2
  )
  
  #hthc
  print(hthc)
dev.off()

# 10,3,13
# df <- data.frame(Cluster = colData(cds.mn.0912)$cluster[match(names(pseudotime(cds.mn.0912)[order(pseudotime(cds.mn.0912))]),
#                                                               names(colData(cds.mn.0912)$cluster))], 
#                  Pseudotime = pseudotime(cds.mn.0912)[order(pseudotime(cds.mn.0912))])
df <- df.all[df.all$Cluster %in% c(10,3,13),]

ClusterCol <- monpn.col[c(10,3,13)]
#names(ClusterCol) <- sort(unique(colData(cds.mn.0912)$Seurat_res0.2))
names(ClusterCol) <- c(10,3,13)
ClusterCol

max_pseudo <- ceiling(max(df$Pseudotime))

#ha <- HeatmapAnnotation(df = df,col = list(Cluster = ClusterCol, 
#                                           Pseudotime = colorRamp2(seq(from=0,to=max_pseudo,length=11),rev(brewer.pal(11, "RdBu")))),
#                        annotation_name_side = "left")


df1 <- subset(df, Cluster %in% c(10,3,13), select=Cluster )
df2 <- subset(df, Cluster %in% c(10,3,13), select=Pseudotime )

ha1 <- HeatmapAnnotation(df = df1, col = list(Cluster = ClusterCol), 
                         annotation_name_side = "left")
ha2 <- HeatmapAnnotation(df = df2, col = list(Pseudotime = colorRamp2(seq(from=0,to=max_pseudo,length=11),rev(brewer.pal(11, "RdBu")))), 
                         annotation_name_side = "left")


marker <- read.table("20220117.MN0912.curve.txt", header = F)$V1

genes <- marker[marker %in% rownames(rowData(cds.mn.0912))]

#
pt.matrix <- exprs(cds.mn.0912)[match(genes,rownames(rowData(cds.mn.0912))),
                                match(rownames(df),colnames(exprs(cds.mn.0912)))]
#                                order(pseudotime(cds.mn.0912)[names(pseudotime(cds.mn.0912)) %in% rownames(df)])]
pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes

#Ward.D2 Hierarchical Clustering
pdf(paste("20220117_MN0912.to13.rainbowheatmap.pdf",sep=""),width=12,height = length(genes)/3)
hthc <- Heatmap(
  as.matrix(pt.matrix),
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  top_annotation               = ha1,
  bottom_annotation            = ha2
)

#hthc
print(hthc)
dev.off()

# 10,2,3,12
#df <- data.frame(Cluster = colData(cds.mn.0912)$cluster[match(names(pseudotime(cds.mn.0912)[order(pseudotime(cds.mn.0912))]),
#                                                             names(colData(cds.mn.0912)$cluster))], 
#                 Pseudotime = pseudotime(cds.mn.0912)[order(pseudotime(cds.mn.0912))])
df <- df.all[df.all$Cluster %in% c(10,2,3,12),]

ClusterCol <- monpn.col[c(10,2,3,12)]
#names(ClusterCol) <- sort(unique(colData(cds.mn.0912)$Seurat_res0.2))
names(ClusterCol) <- c(10,2,3,12)
ClusterCol

max_pseudo <- ceiling(max(df$Pseudotime))

ha <- HeatmapAnnotation(df = df,col = list(Cluster = ClusterCol, 
                                           Pseudotime = colorRamp2(seq(from=0,to=max_pseudo,length=11),rev(brewer.pal(11, "RdBu")))),
                        annotation_name_side = "left")


df1 <- subset(df, Cluster %in% c(10,2,3,12), select=Cluster )
df2 <- subset(df, Cluster %in% c(10,2,3,12), select=Pseudotime )

ha1 <- HeatmapAnnotation(df = df1, col = list(Cluster = ClusterCol), 
                         annotation_name_side = "left")
ha2 <- HeatmapAnnotation(df = df2, col = list(Pseudotime = colorRamp2(seq(from=0,to=max_pseudo,length=11),rev(brewer.pal(11, "RdBu")))), 
                         annotation_name_side = "left")


marker <- read.table("20220117.MN0912.curve.txt", header = F)$V1

genes <- marker[marker %in% rownames(rowData(cds.mn.0912))]

#
pt.matrix <- exprs(cds.mn.0912)[match(genes,rownames(rowData(cds.mn.0912))),
                                match(rownames(df),colnames(exprs(cds.mn.0912)))]
#                                order(pseudotime(cds.mn.0912)[names(pseudotime(cds.mn.0912)) %in% rownames(df)])]
pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes

#Ward.D2 Hierarchical Clustering
pdf(paste("20220117_MN0912.to12.rainbowheatmap.pdf",sep=""),width=12,height = length(genes)/3)
hthc <- Heatmap(
  as.matrix(pt.matrix),
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  top_annotation               = ha1,
  bottom_annotation            = ha2
)

#hthc
print(hthc)
dev.off()

##### 20220119 IN dotplot try violin
# use RenameIdents and SetIdents
project.insplit@active.ident <- project.insplit@meta.data$RNA_snn_res.3
idents<-as.factor(as.numeric(as.vector(project.insplit@active.ident))+1)
names(idents) <- rownames(project.insplit@meta.data)
project.insplit@active.ident <- idents
project.insplit@active.ident

cells.use <- WhichCells(project.insplit, idents = c(17,4,24))
project.insplit <- SetIdent(project.insplit, cells = cells.use, value = 17)

cells.use <- WhichCells(project.insplit, idents = c(8,25,10,7,9))
project.insplit <- SetIdent(project.insplit, cells = cells.use, value = 8)

cells.use <- WhichCells(project.insplit, idents = c(21,13))
project.insplit <- SetIdent(project.insplit, cells = cells.use, value = 21)

cells.use <- WhichCells(project.insplit, idents = c(2,3,14))
project.insplit <- SetIdent(project.insplit, cells = cells.use, value = 2)

cells.use <- WhichCells(project.insplit, idents = c(6,15))
project.insplit <- SetIdent(project.insplit, cells = cells.use, value = 6)

cells.use <- WhichCells(project.insplit, idents = c(22,26))
project.insplit <- SetIdent(project.insplit, cells = cells.use, value = 22)

cells.use <- WhichCells(project.insplit, idents = c(18,27))
project.insplit <- SetIdent(project.insplit, cells = cells.use, value = 18)

cells.use <- WhichCells(project.insplit, idents = c(1,5,11,23))
project.insplit <- SetIdent(project.insplit, cells = cells.use, value = 1)

project.insplit@active.ident <- factor(project.insplit@active.ident, levels = c(20,16,17,8,21,2,6,22,18,1,28,12,19))

### marker1 
marker <- read.table("20220119.IN.marker1.txt", header = F, sep = "\t")
plot.data <- project.insplit@assays$RNA@scale.data[rownames(project.insplit@assays$RNA@scale.data) %in% marker$V1,]
dim(plot.data)
plot.data <- as.data.frame(plot.data)
plot.data[c(1:5),c(1:5)]

plot.data <- t(plot.data)
plot.data <- as.data.frame(plot.data)
plot.data[c(1:5),c(1:5)]
dim(plot.data)

##### merge cluster
length(project.insplit@active.ident)
plot.data.final <- plot.data

plot.data.final$cluster <- as.factor(project.insplit@active.ident)
plot.data.final$barcode <- rownames(plot.data.final)
plot.data.final$sample <- project.insplit@meta.data$Sample
dim(plot.data.final)

plot.data.final <- reshape2::melt(plot.data.final,id.vars=c("barcode","cluster","sample"))
head(plot.data.final)
colnames(plot.data.final) <- c("barcode","cluster","sample","symbol","expression")
plot.data.final$symbol <- factor(plot.data.final$symbol, levels = marker$V1)
plot.data.final$cluster <- factor(plot.data.final$cluster, 
                                  levels = rev(c(20,16,17,8,21,2,6,22,18,1,28,12,19)))

plot.data.final$type[plot.data.final$cluster %in% c(20,16,17,8,21,2,6)] <- "Dorsal"
plot.data.final$type[plot.data.final$cluster %in% c(18,1,28,12,19)] <- "Ventral"
plot.data.final$type[plot.data.final$cluster %in% c(22)] <- "Mix"

#color <- col[c(1:8)][c(3,4,5,7,8,2,1,6)]

p<-ggplot(plot.data.final)+
  geom_violin(aes(x=expression,y=cluster,fill = type),scale="width")+
  theme_bw()+ ylab("")+xlab("")+
  facet_grid(.~symbol) +
  theme(plot.margin=unit(x=c(0.05,0,0,0),units="inches"),strip.background = element_blank()) +
  theme(strip.text.x = element_text(angle = 90)) +
  theme(axis.line.y = element_line(color = "black"),axis.title.y=element_text(size=8)) +
  theme(axis.line.x = element_blank(),axis.ticks.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid=element_blank()) +
  scale_fill_manual(values = col[c(5,2,3)], name = "Celltype")
p
ggsave("20220119.IN.miniviolin1.pdf",p,width = 30,height = 10)

### marker2
marker <- read.table("20220119.IN.marker2.txt", header = F, sep = "\t")
plot.data <- project.insplit@assays$RNA@scale.data[rownames(project.insplit@assays$RNA@scale.data) %in% marker$V1,]
dim(plot.data)
plot.data <- as.data.frame(plot.data)
plot.data[c(1:5),c(1:5)]

plot.data <- t(plot.data)
plot.data <- as.data.frame(plot.data)
plot.data[c(1:5),c(1:5)]
dim(plot.data)

##### merge cluster
length(project.insplit@active.ident)
plot.data.final <- plot.data

plot.data.final$cluster <- as.factor(project.insplit@active.ident)
plot.data.final$barcode <- rownames(plot.data.final)
plot.data.final$sample <- project.insplit@meta.data$Sample
dim(plot.data.final)

plot.data.final <- reshape2::melt(plot.data.final,id.vars=c("barcode","cluster","sample"))
head(plot.data.final)
colnames(plot.data.final) <- c("barcode","cluster","sample","symbol","expression")
plot.data.final$symbol <- factor(plot.data.final$symbol, levels = marker$V1)
plot.data.final$cluster <- factor(plot.data.final$cluster, 
                                  levels = rev(c(20,16,17,8,21,2,6,22,18,1,28,12,19)))

plot.data.final$type[plot.data.final$cluster %in% c(20,16,17,8,21,2,6)] <- "Dorsal"
plot.data.final$type[plot.data.final$cluster %in% c(18,1,28,12,19)] <- "Ventral"
plot.data.final$type[plot.data.final$cluster %in% c(22)] <- "Mix"

#color <- col[c(1:8)][c(3,4,5,7,8,2,1,6)]

p<-ggplot(plot.data.final)+
  geom_violin(aes(x=expression,y=cluster,fill = type),scale="width")+
  theme_bw()+ ylab("")+xlab("")+
  facet_grid(.~symbol) +
  theme(plot.margin=unit(x=c(0.05,0,0,0),units="inches"),strip.background = element_blank()) +
  theme(strip.text.x = element_text(angle = 90)) +
  theme(axis.line.y = element_line(color = "black"),axis.title.y=element_text(size=8)) +
  theme(axis.line.x = element_blank(),axis.ticks.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid=element_blank()) +
  scale_fill_manual(values = col[c(5,2,3)], name = "Celltype")
p
ggsave("20220119.IN.miniviolin2.pdf",p,width = 30,height = 10)

### marker3
marker <- read.table("20211227.INsplit.marker1.txt", header = T, sep = "\t")
plot.data <- project.insplit@assays$RNA@scale.data[rownames(project.insplit@assays$RNA@scale.data) %in% marker$gene,]
dim(plot.data)
plot.data <- as.data.frame(plot.data)
plot.data[c(1:5),c(1:5)]

plot.data <- t(plot.data)
plot.data <- as.data.frame(plot.data)
plot.data[c(1:5),c(1:5)]
dim(plot.data)

##### merge cluster
length(project.insplit@active.ident)
plot.data.final <- plot.data

plot.data.final$cluster <- as.factor(project.insplit@active.ident)
plot.data.final$barcode <- rownames(plot.data.final)
plot.data.final$sample <- project.insplit@meta.data$Sample
dim(plot.data.final)

plot.data.final <- reshape2::melt(plot.data.final,id.vars=c("barcode","cluster","sample"))
head(plot.data.final)
colnames(plot.data.final) <- c("barcode","cluster","sample","symbol","expression")
plot.data.final$symbol <- factor(plot.data.final$symbol, levels = unique(marker$gene))
plot.data.final$cluster <- factor(plot.data.final$cluster, 
                                  levels = rev(c(20,16,17,8,21,2,6,22,18,1,28,12,19)))

plot.data.final$type[plot.data.final$cluster %in% c(20,16,17,8,21,2,6)] <- "Dorsal"
plot.data.final$type[plot.data.final$cluster %in% c(18,1,28,12,19)] <- "Ventral"
plot.data.final$type[plot.data.final$cluster %in% c(22)] <- "Mix"

#color <- col[c(1:8)][c(3,4,5,7,8,2,1,6)]

p<-ggplot(plot.data.final)+
  geom_violin(aes(x=expression,y=cluster,fill = type),scale="width")+
  theme_bw()+ ylab("")+xlab("")+
  facet_grid(.~symbol) +
  theme(plot.margin=unit(x=c(0.05,0,0,0),units="inches"),strip.background = element_blank()) +
  theme(strip.text.x = element_text(angle = 90)) +
  theme(axis.line.y = element_line(color = "black"),axis.title.y=element_text(size=8)) +
  theme(axis.line.x = element_blank(),axis.ticks.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid=element_blank()) +
  scale_fill_manual(values = col[c(5,2,3)], name = "Celltype")
p
ggsave("20220119.IN.miniviolin3.pdf",p,width = 120,height = 10, limitsize = F)

### marker4
marker <- read.table("20220119.IN.marker4.txt", header = F, sep = "\t")
plot.data <- project.insplit@assays$RNA@scale.data[rownames(project.insplit@assays$RNA@scale.data) %in% marker$V1,]
dim(plot.data)
plot.data <- as.data.frame(plot.data)
plot.data[c(1:5),c(1:5)]

plot.data <- t(plot.data)
plot.data <- as.data.frame(plot.data)
plot.data[c(1:5),c(1:5)]
dim(plot.data)

##### merge cluster
length(project.insplit@active.ident)
plot.data.final <- plot.data

plot.data.final$cluster <- as.factor(project.insplit@active.ident)
plot.data.final$barcode <- rownames(plot.data.final)
plot.data.final$sample <- project.insplit@meta.data$Sample
dim(plot.data.final)

plot.data.final <- reshape2::melt(plot.data.final,id.vars=c("barcode","cluster","sample"))
head(plot.data.final)
colnames(plot.data.final) <- c("barcode","cluster","sample","symbol","expression")
plot.data.final$symbol <- factor(plot.data.final$symbol, levels = rev(unique(marker$V1)))
plot.data.final$cluster <- factor(plot.data.final$cluster, 
                                  levels = rev(c(20,16,17,8,21,2,6,22,18,1,28,12,19)))

plot.data.final$type[plot.data.final$cluster %in% c(20,16,17,8,21,2,6)] <- "Dorsal"
plot.data.final$type[plot.data.final$cluster %in% c(18,1,28,12,19)] <- "Ventral"
plot.data.final$type[plot.data.final$cluster %in% c(22)] <- "Mix"

#color <- col[c(1:8)][c(3,4,5,7,8,2,1,6)]

p<-ggplot(plot.data.final)+
  geom_violin(aes(x=expression,y=cluster,fill = type),scale="width")+
  theme_bw()+ ylab("")+xlab("")+
  facet_grid(.~symbol) +
  theme(plot.margin=unit(x=c(0.05,0,0,0),units="inches"),strip.background = element_blank()) +
  theme(strip.text.x = element_text(angle = 270, face = "italic")) +
  theme(axis.line.y = element_line(color = "black")) +
  theme(axis.line.x = element_blank(),axis.ticks = element_blank()) +
  theme(axis.text.x = element_blank(),axis.text.y =element_text(size=8,color = "black", angle = 270, hjust = 0.5,vjust = 0.5)) +
  theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid=element_blank()) +
  scale_fill_manual(values = col[c(5,2,3)], name = "Celltype")
p
ggsave("20220119.IN.miniviolin4.pdf",p,width = 30,height = 10, limitsize = F)

##### 20220119 PG scenic heatmap
### all heatmap
marker1 <- read.table("20220117.PGsplit.pIN.txt", header = F)$V1
marker2 <- read.table("20220117.PGsplit.pMN.txt", header = F)$V1
marker <- c(marker1,marker2)[c(1:53)]

cellInfo <- data.frame(cellid=names(project.pgsplit@active.ident),CellType=project.pgsplit@active.ident)
#cellInfo <- data.frame(cellid=names(cds.pgsplit@clusters$UMAP$clusters),CellType=cds.pgsplit@clusters$UMAP$clusters)
cellInfo <- data.frame(cellInfo)
cellInfo <- cellInfo[rownames(cellInfo) %in% colnames(exprMat),]
head(cellInfo)

cellInfo <- data.frame(CellType=cellInfo$CellType, row.names = cellInfo$cellid)
head(cellInfo)
cellInfo$CellType <- factor(cellInfo$CellType, levels = sort(unique(cellInfo$CellType)))

#
#regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC))[marker],]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$CellType),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
#regulonActivity_byCellType <- regulonActivity_byCellType[,c(3,5,2,1,4)]

regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
#regulonActivity_byCellType_Scaled <- regulonActivity_byCellType_Scaled[,c(3,5,2,1,4)]
#regulonActivity_byCellType <- regulonActivity_byCellType[,c(3,5,2,1,4)]

#ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled[,c(9,7,8,3,6,10,2,4,5,1)], name="Regulon activity\nz-score", 
#                        cluster_rows = F, cluster_columns = F )
p <- ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled[,c(6,1,8,7,4,3,5,2)], name="Regulon activity\nz-score", 
                             cluster_rows = F, cluster_columns = F )


pdf("20220119.PGsplit.SCENIC.heatmap.pdf",height = 10,width = 8)
print(p)
dev.off()

##### 20220119 PGsplit pseudotime curve
marker <- unique(read.table("20220119.PGsplit.marker1.txt", header = F, sep = "\t")$V1)

pt.matrix <- exprs(cds.pgsplit)[match(marker,rownames(rowData(cds.pgsplit))),]
#pt.matrix <- project.pgsplit@assays$RNA@scale.data[rownames(project.pgsplit@assays$RNA@scale.data) %in% marker$V1,]
pt.matrix <- as.data.frame(t(pt.matrix))
dim(pt.matrix)
head(pt.matrix)

pt.matrix$pseudotime <- pseudotime(cds.pgsplit)
pt.matrix$cluster <- colData(cds.pgsplit)$Seurat_res0.2
#pt.matrix$sample <- pData(cds.pgsplit)$orig.ident
pt.matrix$celltype <- "Progenitor"
pt.matrix$celltype[pt.matrix$cluster %in% c(1,6)] <- "pIN"
pt.matrix$celltype[pt.matrix$cluster == 2] <- "pMN"

#pt.matrix <- reshape2::melt(pt.matrix)
head(pt.matrix)

#data <- pt.matrix[pt.matrix$variable %in% marker$V1,]
data <- pt.matrix[pt.matrix$celltype %in% c("pIN","pMN"),]
data <- reshape2::melt(data, id=c("pseudotime", "cluster","celltype"))
data$celltype <- factor(data$celltype, levels = c("pMN","pIN"))

p<-ggplot(data) + geom_point(aes(x=pseudotime, y=value, color = celltype),size = 0.1) +
  geom_smooth(aes(x=pseudotime, y=value, color = celltype),size = 0.5) + scale_color_brewer(palette = "Set1") +
  theme_bw()+  facet_wrap(.~variable,ncol=4,scales = "free_y") +  scale_y_log10() +
  theme(panel.grid = element_blank(), strip.background = element_blank()) +
  xlab("Pesudotime") + ylab("Expression/log10")
p

pdf("20220119.PGsplit.pseudotime_curve1.pdf",width=10,height = length(marker)/4)
  print(p)
dev.off()

##### 20220129 GO1
library(clusterProfiler)
library(org.Hs.eg.db)

marker <- read.table("20220129.GO1.txt", header = F, sep = "\t")
go1 <- enrichGO(marker$V1, OrgDb='org.Hs.eg.db', ont = "BP", keyType = "SYMBOL")
write.table(go1, "20220129.GO1.xls", sep = "\t", quote = F, row.names = F)

##### 20220129 GO2
marker <- read.table("20220129.GO2.txt", header = F, sep = "\t")
go2 <- enrichGO(marker$V1, OrgDb='org.Hs.eg.db', ont = "BP", keyType = "SYMBOL")
write.table(go2, "20220129.GO2.xls", sep = "\t", quote = F, row.names = F)

##### SCENIC Heatmap
library(ComplexHeatmap)
library(pheatmap)
library(SCENIC)
library(AUCell)
library(RColorBrewer)
library(circlize)
load("/data2/chenzixi/20210419_WJ_SC/20210831_s4hV6Sub_SCENIC/MN0912/20210912_s4hSubV6_MN_SCENIC.Final.RData")
marker <- read.table("20220129.SCENIC.heatmapTF.txt", header = T)

cellInfo <- data.frame(cellid=names(cds.mn.0912@clusters$UMAP$clusters),CellType=cds.mn.0912@clusters$UMAP$clusters)
cellInfo <- data.frame(cellInfo)
cellInfo <- cellInfo[rownames(cellInfo) %in% colnames(exprMat),]
head(cellInfo)

cellInfo <- data.frame(CellType=cellInfo$CellType, row.names = cellInfo$cellid)
head(cellInfo)
cellInfo$CellType <- factor(cellInfo$CellType, levels = sort(unique(cellInfo$CellType)))

regulonAUC.all <- regulonAUC

# HMC
regulonAUC <- regulonAUC.all
#regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC))[marker$TF[marker$Type == "HMC"]],]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$CellType),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
#regulonActivity_byCellType <- regulonActivity_byCellType[,c(3,5,2,1,4)]

regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
regulonActivity_byCellType_Scaled <- regulonActivity_byCellType_Scaled[,c(3,5,2,1,4)]
#regulonActivity_byCellType <- regulonActivity_byCellType[,c(3,5,2,1,4)]

pdf("20220129.MN0912.HMC.regulon.pdf",width = 5, height = dim(regulonActivity_byCellType_Scaled)[2]*0.75)
ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity\nz-score", 
                        cluster_rows = F, cluster_columns = F)
dev.off()

# LMC
regulonAUC <- regulonAUC.all
#regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC))[marker$TF[marker$Type == "LMC"]],]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$CellType),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
#regulonActivity_byCellType <- regulonActivity_byCellType[,c(3,5,2,1,4)]

regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
regulonActivity_byCellType_Scaled <- regulonActivity_byCellType_Scaled[,c(3,5,2,1,4)]
#regulonActivity_byCellType <- regulonActivity_byCellType[,c(3,5,2,1,4)]

pdf("20220129.MN0912.LMC.regulon.pdf",width = 5, height = dim(regulonActivity_byCellType_Scaled)[2]*0.75)
ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity\nz-score", 
                        cluster_rows = F, cluster_columns = F)
dev.off()

# MMC
regulonAUC <- regulonAUC.all
#regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC))[marker$TF[marker$Type == "MMC"]],]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$CellType),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
#regulonActivity_byCellType <- regulonActivity_byCellType[,c(3,5,2,1,4)]

regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
regulonActivity_byCellType_Scaled <- regulonActivity_byCellType_Scaled[,c(3,5,2,1,4)]
#regulonActivity_byCellType <- regulonActivity_byCellType[,c(3,5,2,1,4)]

pdf("20220129.MN0912.MMC.regulon.pdf",width = 5, height = dim(regulonActivity_byCellType_Scaled)[2]*0.75)
ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity\nz-score", 
                        cluster_rows = F, cluster_columns = F)
dev.off()

# All
regulonAUC <- regulonAUC.all
#regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC))[marker$TF],]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$CellType),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
#regulonActivity_byCellType <- regulonActivity_byCellType[,c(3,5,2,1,4)]

regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
regulonActivity_byCellType_Scaled <- regulonActivity_byCellType_Scaled[,c(3,5,2,1,4)]
#regulonActivity_byCellType <- regulonActivity_byCellType[,c(3,5,2,1,4)]

pdf("20220129.MN0912.All.regulon.pdf",width = 7, height = dim(regulonActivity_byCellType_Scaled)[2]*0.8)
ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, 
                        name="Regulon activity\nz-score", 
                        col = colorRamp2(c(-2, 0, 2), c("darkseagreen", "white", "red")),
                        cluster_rows = F, cluster_columns = F)
dev.off()

### Violin1
marker <- read.table("20220129.violin1.txt", header = F, sep = "\t")
plot.data <- project.insplit@assays$RNA@scale.data[rownames(project.insplit@assays$RNA@scale.data) %in% marker$V1,]
dim(plot.data)
plot.data <- as.data.frame(plot.data)
plot.data[c(1:5),c(1:5)]

plot.data <- t(plot.data)
plot.data <- as.data.frame(plot.data)
plot.data[c(1:5),c(1:5)]
dim(plot.data)

##### merge cluster
length(project.insplit@active.ident)
plot.data.final <- plot.data

plot.data.final$cluster <- as.factor(project.insplit@active.ident)
plot.data.final$barcode <- rownames(plot.data.final)
plot.data.final$sample <- project.insplit@meta.data$Sample
dim(plot.data.final)

plot.data.final <- reshape2::melt(plot.data.final,id.vars=c("barcode","cluster","sample"))
head(plot.data.final)
colnames(plot.data.final) <- c("barcode","cluster","sample","symbol","expression")
plot.data.final$symbol <- factor(plot.data.final$symbol, levels = rev(unique(marker$V1)))
plot.data.final$cluster <- factor(plot.data.final$cluster, 
                                  levels = rev(c(20,16,17,8,21,2,6,22,18,1,28,12,19)))

plot.data.final$type[plot.data.final$cluster %in% c(20,16,17,8,21,2,6)] <- "Dorsal"
plot.data.final$type[plot.data.final$cluster %in% c(18,1,28,12,19)] <- "Ventral"
plot.data.final$type[plot.data.final$cluster %in% c(22)] <- "Mix"

#color <- col[c(1:8)][c(3,4,5,7,8,2,1,6)]

p<-ggplot(plot.data.final)+
  geom_violin(aes(x=expression,y=cluster,fill = type),scale="width")+
  theme_bw()+ ylab("")+xlab("")+
  facet_grid(.~symbol) +
  theme(plot.margin=unit(x=c(0.05,0,0,0),units="inches"),strip.background = element_blank()) +
  theme(strip.text.x = element_text(angle = 270, face = "italic")) +
  theme(axis.line.y = element_line(color = "black")) +
  theme(axis.line.x = element_blank(),axis.ticks = element_blank()) +
  theme(axis.text.x = element_blank(),axis.text.y =element_text(size=8,color = "black", angle = 270, hjust = 0.5,vjust = 0.5)) +
  theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid=element_blank()) +
  scale_fill_manual(values = col[c(5,2,3)], name = "Celltype")
p
ggsave("20220129.IN.miniviolin1.pdf",p,width = 6,height = 4, limitsize = F)

### Violin2
plot.data <- project.insplit@assays$RNA@scale.data[rownames(project.insplit@assays$RNA@scale.data) %in% marker$V1,]
dim(plot.data)
plot.data <- as.data.frame(plot.data)
plot.data[c(1:5),c(1:5)]

plot.data <- t(plot.data)
plot.data <- as.data.frame(plot.data)
plot.data[c(1:5),c(1:5)]
dim(plot.data)

##### merge cluster
length(project.insplit@active.ident)
plot.data.final <- plot.data

plot.data.final$cluster <- as.factor(project.insplit@active.ident)
plot.data.final$barcode <- rownames(plot.data.final)
plot.data.final$sample <- project.insplit@meta.data$Sample
dim(plot.data.final)

plot.data.final <- reshape2::melt(plot.data.final,id.vars=c("barcode","cluster","sample"))
head(plot.data.final)
colnames(plot.data.final) <- c("barcode","cluster","sample","symbol","expression")
plot.data.final$symbol <- factor(plot.data.final$symbol, levels = rev(unique(marker$V1)))
plot.data.final$cluster <- factor(plot.data.final$cluster, 
                                  levels = rev(c(20,16,17,8,21,2,6,22,18,1,28,12,19)))

plot.data.final$type[plot.data.final$cluster %in% c(20,16,17,8,21,2,6)] <- "Dorsal"
plot.data.final$type[plot.data.final$cluster %in% c(18,1,28,12,19)] <- "Ventral"
plot.data.final$type[plot.data.final$cluster %in% c(22)] <- "Mix"

#color <- col[c(1:8)][c(3,4,5,7,8,2,1,6)]

p<-ggplot(plot.data.final)+
  geom_violin(aes(x=expression,y=cluster,fill = type),scale="width")+
  theme_bw()+ ylab("")+xlab("")+
  facet_grid(.~symbol) +
  theme(plot.margin=unit(x=c(0.05,0,0,0),units="inches"),strip.background = element_blank()) +
  theme(strip.text.x = element_text(angle = 270, face = "italic")) +
  theme(axis.line.y = element_line(color = "black")) +
  theme(axis.line.x = element_blank(),axis.ticks = element_blank()) +
  theme(axis.text.x = element_blank(),axis.text.y =element_text(size=8,color = "black", angle = 270, hjust = 0.5,vjust = 0.5)) +
  theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid=element_blank()) +
  scale_fill_manual(values = col[c(5,2,3)], name = "Celltype")
p
ggsave("20220129.IN.miniviolin2.pdf",p,width = 6,height = 4, limitsize = F)

#### cds.monpn gsea
monpn.wilcox.diff <- read.table("20211207.wilcoxauc.cds.monpn.cluster.xls", header = T, sep = "\t")

clusters.monpn <- unique(monpn.wilcox.diff$group)
clusters.monpn

pdf("20220130.cdsmonpn.gsea.selected.pdf",width = 6,height = 4)
for (cluster in clusters.monpn){
  data <- monpn.wilcox.diff %>% subset(group == cluster) %>% arrange(desc(logFC)) %>% dplyr::select(feature, logFC)
  ranks <- deframe(data)
  ranks
  
  ranks.gsea <- gseGO(geneList     = ranks,
                       keyType = "SYMBOL",
                       OrgDb        = org.Hs.eg.db,
                       ont          = "BP",
                       pvalueCutoff = 1,
                       eps = 0,
                       pAdjustMethod= "BH",
                       verbose      = TRUE,
                       seed         = TRUE)
  write.table(ranks.gsea,paste("20220130.cdsmonpn.gsea",cluster,"xls",sep = "."), sep = "\t", quote = F, col.names = T, row.names = F)
  gseaplot2(ranks.gsea, c("GO:0030509","GO:0007224","GO:0060071","GO:0008589"), subplots = 1:2,title = cluster)
}
dev.off()

## try "GO:0030509","GO:0007224" only
GO_DATA_NRVC <- clusterProfiler:::get_GO_data("org.Hs.eg.db", "ALL", "SYMBOL")
PATHID2NAME <- as.data.frame(GO_DATA_NRVC$PATHID2NAME)
PATHID2NAME$GO <- rownames(PATHID2NAME)
colnames(PATHID2NAME) <- c("Description","GOTerm")

#GO_BP <- GO_DATA_NRVC$PATHID2NAME[GO_DATA_NRVC$GO2ONT == "BP"]
#Wnt_NRVCGO<-names(GO_DATA_NRVC$PATHID2NAME[grep("Wnt", GO_DATA_NRVC$PATHID2NAME)])
genelist <- plyr::ldply(GO_DATA_NRVC$PATHID2EXTID, data.frame)
#genelist.bp <- genelist[genelist[,1] %in% names(GO_BP),]

genelist.spec <- genelist[genelist[,1] %in% c("GO:0030509","GO:0007224","GO:0060071","GO:0008589"),]

pdf("20220130.cdsmonpn.gsea.spec.pdf",width = 6,height = 4)
for (cluster in clusters.monpn){
  data <- monpn.wilcox.diff %>% subset(group == cluster) %>% arrange(desc(logFC)) %>% dplyr::select(feature, logFC)
  ranks <- deframe(data)
  ranks
  
  ranks.gsea.spec <- clusterProfiler::GSEA(ranks, TERM2GENE = genelist.spec, TERM2NAME = PATHID2NAME[,c(2,1)],
                                           pvalueCutoff = 1,
                                           eps = 0,
                                           pAdjustMethod= "BH",
                                           verbose      = TRUE,
                                           seed         = TRUE)
  write.table(ranks.gsea.spec,paste("20220130.cdsmonpn.gsea",cluster,"spec.xls",sep = "."), sep = "\t", quote = F, col.names = T, row.names = F)
  p<-gseaplot2(ranks.gsea.spec, c("GO:0030509","GO:0007224","GO:0060071","GO:0008589"), subplots = 1:2,title = cluster)
  print(p)
}
dev.off()

##### MN0912 scenic heatmap
library(ComplexHeatmap)
library(pheatmap)
load("/data2/chenzixi/20210419_WJ_SC/20210831_s4hV6Sub_SCENIC/MN0912/20210912_s4hSubV6_MN_SCENIC.Final.RData")

cellInfo <- data.frame(cellid=names(cds.mn.0912@clusters$UMAP$clusters),CellType=cds.mn.0912@clusters$UMAP$clusters)
cellInfo <- data.frame(cellInfo)
cellInfo <- cellInfo[rownames(cellInfo) %in% colnames(exprMat),]
head(cellInfo)

cellInfo <- data.frame(CellType=cellInfo$CellType, row.names = cellInfo$cellid)
head(cellInfo)
cellInfo$CellType <- factor(cellInfo$CellType, levels = sort(unique(cellInfo$CellType)))

regulonAUC.all <- regulonAUC


regulonAUC <- regulonAUC.all
#regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC))]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$CellType),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
#regulonActivity_byCellType <- regulonActivity_byCellType[,c(3,5,2,1,4)]

regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
#regulonActivity_byCellType_Scaled <- regulonActivity_byCellType_Scaled[,c(3,5,2,1,4)]
#regulonActivity_byCellType <- regulonActivity_byCellType[,c(3,5,2,1,4)]

pdf("20220130.MN0912.heatmap.pdf",width = 5, height = dim(regulonActivity_byCellType_Scaled)[2]*2)
ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity\nz-score", 
                        cluster_rows = T, cluster_columns = T,
                        row_names_gp = gpar(fontsize = 6),
                        column_names_gp = gpar(fontsize = 6))
pheatmap::pheatmap(regulonActivity_byCellType_Scaled, 
                   color=colorRampPalette(c("blue","white","red"))(100), 
                   breaks=seq(-3, 3, length.out = 100),treeheight_row=10,
                   treeheight_col=10, border_color=NA,fontsize_col = 6, fontsize_row = 6)
dev.off()

pdf("20220130.MN0912.heatmap.noname.pdf",width = 5, height = dim(regulonActivity_byCellType_Scaled)[2]*2)
ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity\nz-score", 
                        cluster_rows = T, cluster_columns = T,
                        column_names_gp = gpar(fontsize = 6) , show_row_names = FALSE)

pheatmap::pheatmap(regulonActivity_byCellType_Scaled, 
                   color=colorRampPalette(c("blue","white","red"))(100), 
                   breaks=seq(-3, 3, length.out = 100),treeheight_row=10,
                   treeheight_col=10, border_color=NA,fontsize_col = 6, show_rownames = FALSE,)
dev.off()

##### IN0912 scenic heatmap
library(ComplexHeatmap)
library(pheatmap)
load("/data2/chenzixi/20210419_WJ_SC/20210831_s4hV6Sub_SCENIC/IN0912/20210912_s4hSubV6_IN_SCENIC.Final.RData")

cellInfo <- data.frame(cellid=names(cds.in.0912@clusters$UMAP$clusters),CellType=cds.in.0912@clusters$UMAP$clusters)
cellInfo <- data.frame(cellInfo)
cellInfo <- cellInfo[rownames(cellInfo) %in% colnames(exprMat),]
head(cellInfo)

cellInfo <- data.frame(CellType=cellInfo$CellType, row.names = cellInfo$cellid)
head(cellInfo)
cellInfo$CellType <- factor(cellInfo$CellType, levels = sort(unique(cellInfo$CellType)))

regulonAUC.all <- regulonAUC


regulonAUC <- regulonAUC.all
#regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC))]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$CellType),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
#regulonActivity_byCellType <- regulonActivity_byCellType[,c(3,5,2,1,4)]

regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
#regulonActivity_byCellType_Scaled <- regulonActivity_byCellType_Scaled[,c(3,5,2,1,4)]
#regulonActivity_byCellType <- regulonActivity_byCellType[,c(3,5,2,1,4)]

pdf("20220130.IN0912.heatmap.pdf",width = 5, height = dim(regulonActivity_byCellType_Scaled)[2]*2)
ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity\nz-score", 
                        cluster_rows = T, cluster_columns = T,
                        row_names_gp = gpar(fontsize = 6),
                        column_names_gp = gpar(fontsize = 6))
pheatmap::pheatmap(regulonActivity_byCellType_Scaled, 
                   color=colorRampPalette(c("blue","white","red"))(100), 
                   breaks=seq(-3, 3, length.out = 100),treeheight_row=10,
                   treeheight_col=10, border_color=NA,fontsize_col = 6, fontsize_row = 6)
dev.off()

pdf("20220130.IN0912.heatmap.noname.pdf",width = 5, height = dim(regulonActivity_byCellType_Scaled)[2]*2)
ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity\nz-score", 
                        cluster_rows = T, cluster_columns = T,
                        column_names_gp = gpar(fontsize = 6) , show_row_names = FALSE)

pheatmap::pheatmap(regulonActivity_byCellType_Scaled, 
                   color=colorRampPalette(c("blue","white","red"))(100), 
                   breaks=seq(-3, 3, length.out = 100),treeheight_row=10,
                   treeheight_col=10, border_color=NA,fontsize_col = 6, show_rownames = FALSE,)
dev.off()

##### mouse ref GSEA cds.monpn
# mouse ref
genelist <- read.table("20211206.GSEA.ref.txt", header = F, sep = "\t")
genelist$Ont <- "Disease"
genelist <- genelist[c(2,1)]
colnames(genelist) <- c("Ont","Gene")
genelist

cds.monpn.mouse.gsea <- data.frame(cluster = clusters.monpn)
cds.monpn.mouse.gsea$NES <- NA
cds.monpn.mouse.gsea$pvalue <- NA

pdf("20220130.cdsmonpn.gsea.mouse.Disease.pdf",width = 6,height = 4)
for (cluster in clusters.monpn){
  data <- monpn.wilcox.diff %>% subset(group == cluster) %>% arrange(desc(logFC)) %>% dplyr::select(feature, logFC)
  ranks <- deframe(data)
  ranks
  
  ranks.gsea.disease <- clusterProfiler::GSEA(ranks, TERM2GENE = genelist, 
                                           pvalueCutoff = 1,
                                           eps = 0,
                                           pAdjustMethod= "BH",
                                           verbose      = TRUE,
                                           seed         = TRUE)
  
  cds.monpn.mouse.gsea$NES[cluster] <- ranks.gsea.disease@result$NES
  cds.monpn.mouse.gsea$pvalue[cluster] <- ranks.gsea.disease@result$pvalue
  
  write.table(ranks.gsea.disease,paste("20220130.cdsmonpn.gsea",cluster,"Disease.xls",sep = "."), sep = "\t", quote = F, col.names = T, row.names = F)
  p<- gseaplot2(ranks.gsea.disease, 1, subplots = 1:2,title = cluster)
  print(p)
}
dev.off()

write.table(cds.monpn.mouse.gsea,"20220130.cdsmonpn.gsea.all.Disease.xls",sep = "\t", quote = F, col.names = T, row.names = F)

cds.monpn.mouse.gsea$cluster <- factor(cds.monpn.mouse.gsea$cluster,
                                       levels = c(6,10,3,2,13,12,4,11,1,9,8,7,5))

monpn.col <- c(
  "#6ab662",
  "#377EB8",
  "#cab3d6",
  "#2ea37e",
  "#ebc4b5",
  "#e11517",
  "#efd092",
  "#c47e70",
  "#b8e294",
  "#774ea3",
  "#ffd92f",
  "#b5d7e8",
  "#685f7e"
)

gsea.col <- monpn.col[c(6,10,3,2,13,12,4,11,1,9,8,7,5)]
ggplot(cds.monpn.mouse.gsea) + geom_bar(aes(cluster, weight=NES, fill=cluster),width=0.75) + theme_classic() +
  xlab("Cluster") + ylab("NES value") + scale_fill_manual(values = gsea.col) + theme(legend.position = "none", axis.text = element_text(color="black"))
  
ggsave("20220130.cdsmonpn.mouse.gsea.NESbar.pdf",height = 3,width = 4)

#### monPN marker
marker <- read.table("marker1.txt",header = F)[,1]
p <- plot_cells(cds, genes = marker, 
                show_trajectory_graph=FALSE, 
                label_cell_groups=FALSE)+
  scale_color_gradient(low = "lightgrey",high= "firebrick")+ coord_fixed()
ggsave("20220130_s4hV6Sub_Monocle3_monPN.marker1.pdf",p,height=32,width = 32)

marker <- read.table("marker2.txt",header = F)[,1]
p<-plot_cells(cds, genes = marker, 
              show_trajectory_graph=FALSE, 
              label_cell_groups=FALSE)+
  scale_color_gradient(low = "lightgrey",high= "firebrick")+ coord_fixed()
ggsave("20220130_s4hV6Sub_Monocle3_monPN.marker2.pdf",p,height=24,width = 32)

marker <- read.table("marker3.txt",header = F)[,1]
p<-plot_cells(cds, genes = unique(marker), 
              show_trajectory_graph=FALSE, 
              label_cell_groups=FALSE)+
  scale_color_gradient(low = "lightgrey",high= "firebrick")+ coord_fixed()
ggsave("20220130_s4hV6Sub_Monocle3_monPN.marker3.pdf",p,height=48,width = 56,limitsize = FALSE)

marker <- read.table("marker4.txt",header = F)[,1]
p<-plot_cells(cds, genes = unique(marker), 
              show_trajectory_graph=FALSE, 
              label_cell_groups=FALSE)+
  scale_color_gradient(low = "lightgrey",high= "firebrick")+ coord_fixed()
ggsave("20220130_s4hV6Sub_Monocle3_monPN.marker4.pdf",p,height=32,width = 32)

#### gsea enrich genes heatmap "20220115.dvgsea.fc.gse_obj.all.xls"
project.insplit@meta.data$DV <- "Mix"
project.insplit@meta.data$DV[project.insplit@meta.data$merge.split %in% c(20,16,17,8,21,2,6)] <- "Dorsal"
project.insplit@meta.data$DV[project.insplit@meta.data$merge.split %in% c(18,1,28,12,19)] <- "Ventral"

expr <- as.data.frame(t(project.insplit@assays$RNA@data))
expr[1:5,1:5]
expr$DV <- project.insplit$DV[match(names(project.insplit$DV),rownames(expr))]
detach(expr)
attach(expr)
cluster_mean.insplit <-  as.data.frame(aggregate(expr, by=list(DV), mean))
rownames(cluster_mean.insplit) <- cluster_mean.insplit$Group.1
cluster_mean.insplit[1:5,1:5]
cluster_mean.insplit <- cluster_mean.insplit[,-1]
cluster_mean.insplit[,1:5]


marker <- read.table("20220130.INsplit.DV.gsea.enrichedGenes.txt", header = F)

plotdata <- cluster_mean.insplit[,colnames(cluster_mean.insplit) %in% marker$V1]
plotdata <- t(plotdata)
plotdata
#colnames(plotdata) <- c(1:8)
plotdata <- plotdata[match(marker$V1,rownames(plotdata)),]

pheatmap(plotdata, cluster_cols = F, cluster_rows = T, scale = "row", filename = "20211014.s4hV6Sub.PGsplit.scenic.heatmap.pdf")

##### 20220211 pMN GO
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(stringr)

data <- read.table("20220211.pMN.GO.txt", header = T, sep = "\t")
data <- data[order(data$Ratio),]
data$Description <- factor(data$Description, levels = data$Description)

ggplot(data) + geom_point(aes(Ratio, Description, fill=-log10(pvalue), size = Count), shape = 21) + theme_bw() +
  theme(axis.text = element_text(color = "black",size=10), axis.ticks.y = element_blank()) +
  scale_fill_gradientn(colors = brewer.pal(9, "YlGnBu")[2:5]) + 
  ylab(NULL)+   scale_y_discrete(labels=(function(x) str_wrap(x, width=25))) +
  coord_cartesian(expand = TRUE)

ggsave("20220211.pMN.GO.pdf", width = 6, height = 4)


##### 20220211 pIN GO
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(stringr)

data <- read.table("20220211.pIN.GO.txt", header = T, sep = "\t")
data <- data[order(data$Ratio),]
data$Description <- factor(data$Description, levels = data$Description)

ggplot(data) + geom_point(aes(Ratio, Description, fill=-log10(pvalue), size = Count), shape = 21) + theme_bw() +
  theme(axis.text = element_text(color = "black",size=10), axis.ticks.y = element_blank()) +
  scale_fill_gradientn(colors = brewer.pal(9, "YlGnBu")[2:5]) + 
  ylab(NULL)+   scale_y_discrete(labels=(function(x) str_wrap(x, width=25))) +
  coord_cartesian(expand = TRUE)

ggsave("20220211.pIN.GO.pdf", width = 6, height = 4)

#####  20220211 PGsplit monocle3
library(monocle3)

cds.pgsplit <- readRDS("./RDS/cds.pgsplit.RDS")

cols=c(brewer.pal(12,"Paired"),brewer.pal(8,"Dark2"),
       brewer.pal(9,"Set1"),brewer.pal(8,"Set2"),
       brewer.pal(12,"Set3"),brewer.pal(8,"Pastel2"),
       brewer.pal(9,"Pastel1"),brewer.pal(8,"Accent"))
col=unique(cols)
col.pgsplit <- col[1:8]
col.pgsplit[2] <- "#ea5151"
col.pgsplit[6] <- "#1F78B3"

p <- plot_cells(cds.pgsplit, reduction_method = "UMAP", graph_label_size=10, show_trajectory_graph = F,
                color_cells_by = "Seurat_res0.2", trajectory_graph_segment_size = 0.5,
                label_cell_groups = F, label_groups_by_cluster = F,
                label_branch_points = F, label_roots = F, label_leaves = F,
                label_principal_points = F) + ggtitle("Seurat_cluster")+
  scale_color_manual(values = col.pgsplit)+ coord_fixed()
p


pdf("20220211_s4hV6Sub2_Monocle3_PGsplit.pdf")
print(p)
dev.off()
