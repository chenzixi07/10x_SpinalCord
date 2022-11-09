library(Seurat)
library(monocle3)
library(patchwork)
library(ggplot2)
library(plyr)
setwd("D:/Projects/wj_sc_project2/20220429_plot")

### Fig 1C
project <- readRDS("D:/Projects/wj_sc_project2/RDS/project.RDS")

genelist <- c("PAX3","PAX6", "ISL1","ISL2", "BHLHE22", "POU4F1")

p1 <- FeaturePlot(project, reduction = "tsne",
                  features = "PAX3",
                  keep.scale = "all",
                  cols = c("Gainsboro", "firebrick")) + coord_fixed(ratio = 1)

p2 <- FeaturePlot(project, reduction = "tsne",
                  features = "PAX6",
                  keep.scale = "all",
                  cols = c("Gainsboro", "firebrick")) + coord_fixed(ratio = 1)

p3 <- FeaturePlot(project, reduction = "tsne",
                  features = "ISL1",
                  keep.scale = "all",
                  cols = c("Gainsboro", "firebrick")) + coord_fixed(ratio = 1)

p4 <- FeaturePlot(project, reduction = "tsne",
                  features = "ISL2",
                  keep.scale = "all",
                  cols = c("Gainsboro", "firebrick")) + coord_fixed(ratio = 1)

p5 <- FeaturePlot(project, reduction = "tsne",
                  features = "BHLHE22",
                  keep.scale = "all",
                  cols = c("Gainsboro", "firebrick")) + coord_fixed(ratio = 1)

p6 <- FeaturePlot(project, reduction = "tsne",
                  features = "POU4F1",
                  keep.scale = "all",
                  cols = c("Gainsboro", "firebrick")) + coord_fixed(ratio = 1)

pp1 <- p1 + p2 + plot_layout(ncol = 2)
pp1

pp2 <- p3 + p4 + plot_layout(ncol = 2)
pp2

pp3 <- p5 + p6 + plot_layout(ncol = 2)
pp3

p <- pp1/pp2/pp3 + plot_layout(ncol = 1, nrow = 3)
p

ggsave("Fig1C.pdf",p,width = 8,height = 12)

### Fig 1D
project <- readRDS("D:/Projects/wj_sc_project2/RDS/project.RDS")

col1 <- rep("Gainsboro", 25)
col1[c(1,3,5,7,9,11,12,13,21,23)] <- "firebrick"

project@active.ident <- factor(project@active.ident,
                               levels= c(1:25))

pdf("Fig1D.pdf")
DimPlot(project, reduction = "tsne" ,cols =col1) +coord_fixed(ratio = 1)
dev.off()

### Fig 1E
cds.monpn <- readRDS("D:/Projects/wj_sc_project2/RDS/cds.monpn.RDS")

col <- c("#4d8532",
         "#e8afd0",
         "#e17da2",
         "#3da436",
         "#abb562",
         "#eb6d6d",
         "#ccd2b4",
         "#a7be8e",
         "#64a61f",
         "#ed56a0",
         "#a3ca3e",
         "#fae1ed",
         "#e4b1ed")

p<-plot_cells(cds.monpn, reduction_method = "UMAP", graph_label_size=5, show_trajectory_graph = F,
           color_cells_by = "cluster", 
           label_cell_groups = F, label_groups_by_cluster = F,
           label_branch_points = F, label_roots = F, label_leaves = F,
           label_principal_points = F) +
  scale_color_manual(values = col) + coord_fixed()
p
ggsave("Fig1E.pdf",p)

##### sup 2 
genelist <- c("PAX3","PAX6","NEUROG2","OLIG2"," NEUROG1","NKX6-1","ISL1","ISL2","MNX1","LHX5","BHLHE22","LMX1B")

p<-plot_cells(cds.monpn, genes = genelist, 
           show_trajectory_graph=FALSE, 
           label_cell_groups=FALSE)+
  scale_color_gradient(low = "Gainsboro",high= "firebrick")+ coord_fixed()

ggsave("FigS2.pdf",p)

### Fig2A
cds.monpn <- readRDS("D:/Projects/wj_sc_project2/RDS/cds.monpn.RDS")

col <- c("#4d8532",
         "#e8afd0",
         "#e17da2",
         "#3da436",
         "#abb562",
         "#eb6d6d",
         "#ccd2b4",
         "#a7be8e",
         "#64a61f",
         "#ed56a0",
         "#a3ca3e",
         "#fae1ed",
         "#e4b1ed")

col2 <- rep("Gainsboro", 13)
col2[4] <- col[4]
col2[6] <- col[6]
col2[10] <- col[10]

p<-plot_cells(cds.monpn, reduction_method = "UMAP", graph_label_size=5, show_trajectory_graph = F,
              color_cells_by = "cluster", 
              label_cell_groups = F, label_groups_by_cluster = F,
              label_branch_points = F, label_roots = F, label_leaves = F,
              label_principal_points = F) +
  scale_color_manual(values = col2) + coord_fixed()
p
ggsave("Fig2A.pdf",p)

### Fig2B
project.pgsplit <- readRDS("D:/Projects/wj_sc_project2/RDS/project.pgsplit.RDS")

col <- c("#94cdf6",
         "#fb9a99",
         "#a3ca3e",
         "#6f892a",
         "#ff7f01",
         "#77a5ca",
         "#fdc884",
         "#ccfb4d")

p<-DimPlot(project.pgsplit, reduction = "umap" ,label=F , cols = col) +coord_fixed()
p
ggsave("Fig2B.pdf",p)

### Fig2E
project.pgsplit <- readRDS("D:/Projects/wj_sc_project2/RDS/project.pgsplit.RDS")
marker <- read.table("Fig2E.txt",sep = "\t", header = F)

plot.data <- project.pgsplit@assays$RNA@scale.data[rownames(project.pgsplit@assays$RNA@scale.data) %in% rev(marker$V1),]
dim(plot.data)
plot.data <- as.data.frame(plot.data)
plot.data[c(1:5),c(1:5)]

plot.data <- t(plot.data)
plot.data <- as.data.frame(plot.data)
plot.data[c(1:5),c(1:5)]
dim(plot.data)

##### merge cluster
length(project.pgsplit@active.ident)
plot.data.final <- plot.data

plot.data.final$cluster <- as.factor(project.pgsplit@active.ident)
plot.data.final$barcode <- rownames(plot.data.final)
plot.data.final$sample <-project.pgsplit@meta.data$Sample
dim(plot.data.final)

plot.data.final <- reshape2::melt(plot.data.final,id.vars=c("barcode","cluster","sample"))
head(plot.data.final)
colnames(plot.data.final) <- c("barcode","cluster","sample","symbol","expression")
plot.data.final$symbol <- factor(plot.data.final$symbol, levels = marker$V1)
plot.data.final$cluster <- factor(plot.data.final$cluster, 
                                  levels = c(3,4,5,7,8,2,1,6))

plot.data.final2 <- plot.data.final[plot.data.final$cluster %in% c(3,4,5,7,8),]
plot.data.final2$cluster <- factor(plot.data.final2$cluster, 
                                   levels = rev(c(5,7,3,4,8)))

col5 <- col[c(1:8)][rev(c(5,7,3,4,8))]

p2 <-ggplot(plot.data.final2)+
  geom_violin(aes(x=expression,y=cluster,fill = cluster),scale="width",width=0.75)+
  theme_bw()+ ylab("")+xlab("")+
  facet_grid(.~symbol) +
  theme(plot.margin=unit(x=c(0.05,0,0,0),units="inches"),strip.background = element_blank()) +
  theme(strip.text.x = element_text(angle = 270,face = "italic")) +
  theme(axis.line.y = element_blank(),axis.text.y=element_text(color = "black",angle = 270)) +
  theme(axis.line.x = element_blank(),axis.ticks = element_blank(), legend.position = "none") +
  theme(axis.text.x = element_blank(),) +
  theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid=element_blank()) +
  scale_fill_manual(values = col5) 
p2

ggsave("Fig2E.pdf",p2,width = 24,height = 8, units = "cm" )

### Fig2F
marker <- c("CLU","SLC1A3","BCAN","MT2A","CST3","HOPX","ASCL1","MKI67","ASPM","CKS1B","UBE2T","CENPF","UBE2C","CCNB1","CCNB2","SPC25","HES6")
cds.plot <- cds.pgsplit[rowData(cds.pgsplit)$gene_short_name %in% marker,]
p<-monocle3::plot_genes_in_pseudotime(cds.plot, color_cells_by = "Seurat_res0.2", ncol = 4, cell_size = 1) + 
  scale_color_manual(values = col)
p
ggsave("Fig2F.pdf",p, width = 16, height = ceiling(length(marker)/4)*2)

### FigS3

p1 <- FeaturePlot(project.pgsplit, reduction = "umap",
                  features = "PAX6",
                  keep.scale = "all",
                  cols = c("Gainsboro", "firebrick")) + coord_fixed(ratio = 1)

p2 <- FeaturePlot(project.pgsplit, reduction = "umap",
                  features = "PAX7",
                  keep.scale = "all",
                  cols = c("Gainsboro", "firebrick")) + coord_fixed(ratio = 1)

p3 <- FeaturePlot(project.pgsplit, reduction = "umap",
                  features = "NKX6-1",
                  keep.scale = "all",
                  cols = c("Gainsboro", "firebrick")) + coord_fixed(ratio = 1)

p4 <- FeaturePlot(project.pgsplit, reduction = "umap",
                  features = "OLIG2",
                  keep.scale = "all",
                  cols = c("Gainsboro", "firebrick")) + coord_fixed(ratio = 1)

p5 <- FeaturePlot(project.pgsplit, reduction = "umap",
                  features = "NEUROG1",
                  keep.scale = "all",
                  cols = c("Gainsboro", "firebrick")) + coord_fixed(ratio = 1)

p6 <- FeaturePlot(project.pgsplit, reduction = "umap",
                  features = "NEUROG2",
                  keep.scale = "all",
                  cols = c("Gainsboro", "firebrick")) + coord_fixed(ratio = 1)

pp1 <- p1 + p2 + plot_layout(ncol = 2)
pp1

pp2 <- p3 + p4 + plot_layout(ncol = 2)
pp2

pp3 <- p5 + p6 + plot_layout(ncol = 2)
pp3

p <- pp1/pp2/pp3 + plot_layout(ncol = 1, nrow = 3)
p

ggsave("FigS3.pdf",p,width = 8,height = 12)

### Fig3A
cds.pgsplit <- readRDS("D:/Projects/wj_sc_project2/RDS/cds.pgsplit.RDS")

p <- plot_cells(cds.pgsplit, reduction_method = "UMAP", graph_label_size=10, show_trajectory_graph = F,
                color_cells_by = "Seurat_res0.2", trajectory_graph_segment_size = 0.5,
                label_cell_groups = F, label_groups_by_cluster = F,
                label_branch_points = F, label_roots = F, label_leaves = F,
                label_principal_points = F) + ggtitle("Seurat_cluster")+
  scale_color_manual(values = col)+ coord_fixed()
p
ggsave("Fig3A.pdf",p)

### cluster ratio 
project <- readRDS("D:/Projects/wj_sc_project2/RDS/project.RDS")

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

project@active.ident <- factor(project@active.ident,
                               levels= c(1:25))

p1 <- DimPlot(project,reduction = "tsne", cols = final_col) + coord_fixed()
p2 <- DimPlot(project,reduction = "umap", cols = final_col) + coord_fixed()
p1 + p2

input <- data.frame(Barcode=colnames(project),
                    Sample=project$Sample,
                    Cluster=project@active.ident)
data <- reshape2::melt(table(input$Sample,input$Cluster))
colnames(data) <- c("Sample","Cluster","Count")
data$Cluster <- ordered(data$Cluster)
data$Sample <- factor(data$Sample,levels = rev(c("SC8", "SC9", "SC10", "SC12")))
data_cluster<-ddply(data,"Sample",transform, Percent = Count / sum(Count))
data_sample<-ddply(data,"Cluster",transform, Percent = Count / sum(Count))


p<-ggplot(data_cluster, aes(x=Sample, y=Percent,fill=Cluster)) + 
  geom_bar(stat="identity",colour = "white", position = "stack",width = 0.8)+
  xlab("Samples")+ylab("Relative percent")+theme(legend.position="right")+
  theme(panel.background=element_blank(),axis.text = element_text(size = 11, colour = "black"), axis.ticks.x=element_blank())+
  scale_fill_manual(values=final_col) +
  theme(legend.key = element_blank())
p <- p + coord_flip()

ggsave("Cluster.sample.ratio.pdf",p,width=10,height=6)

#### pIN split sample violin
library(Seurat)
library(ggplot2)
library(patchwork)

marker <- read.table("pIN.marker.txt", header = F, sep = "\t")

project.insplit <- readRDS("D:/Projects/wj_sc_project2/RDS/project.insplit.RDS")

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

plot.data <- project.insplit@assays$RNA@scale.data[rownames(project.insplit@assays$RNA@scale.data) %in% rev(marker$V1),]
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
plot.data.final$sample <-project.insplit@meta.data$Sample
dim(plot.data.final)

plot.data.final <- reshape2::melt(plot.data.final,id.vars=c("barcode","cluster","sample"))
head(plot.data.final)
colnames(plot.data.final) <- c("barcode","cluster","sample","symbol","expression")
plot.data.final$symbol <- factor(plot.data.final$symbol, levels = marker$V1)
plot.data.final$cluster <- factor(plot.data.final$cluster, 
                                  levels = rev(c(20,16,17,8,21,2,6,22,18,1,28,12,19)))
plot.data.final$type[plot.data.final$cluster %in% c(20,16,17,8,21,2,6)] <- "dorsal"
plot.data.final$type[plot.data.final$cluster %in% c(18,1,28,12,19)] <- "ventral"
plot.data.final$type[plot.data.final$cluster %in% c(22)] <- "mix"

col <- c("#FA9A99", "#1F78B4", "#B1DF89")

p <-ggplot(plot.data.final)+
  geom_violin(aes(x=expression,y=cluster,fill = type),scale="width",width=0.75)+
  theme_bw()+ ylab("")+xlab("")+
  facet_grid(.~symbol) +
  theme(plot.margin=unit(x=c(0.05,0,0,0),units="inches"),strip.background = element_blank()) +
  theme(strip.text.x = element_text(angle = 270,face = "italic")) +
  theme(axis.line.y = element_blank(),axis.text.y=element_text(color = "black",angle = 270)) +
  theme(axis.line.x = element_blank(),axis.ticks = element_blank()) +
  theme(axis.text.x = element_blank(),) +
  theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid=element_blank()) +
  scale_fill_manual(values = col) + ggtitle("All")
 # theme(legend.position = "none")
p

ggsave("INsplit.violin.All.pdf",p,width=12,height=6)

p1 <- ggplot(plot.data.final[plot.data.final$sample == "SC8",])+
  geom_violin(aes(x=expression,y=cluster,fill = type),scale="width",width=0.75)+
  theme_bw()+ ylab("")+xlab("")+
  facet_grid(.~symbol) +
  theme(plot.margin=unit(x=c(0.05,0,0,0),units="inches"),strip.background = element_blank()) +
  theme(strip.text.x = element_text(angle = 270,face = "italic")) +
  theme(axis.line.y = element_blank(),axis.text.y=element_text(color = "black",angle = 270)) +
  theme(axis.line.x = element_blank(),axis.ticks = element_blank()) +
  theme(axis.text.x = element_blank(),) +
  theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid=element_blank()) +
  scale_fill_manual(values = col) +ggtitle("SC8")
p1

ggsave("INsplit.violin.SC8.pdf",p1,width=12,height=6)

p2 <- ggplot(plot.data.final[plot.data.final$sample == "SC9",])+
  geom_violin(aes(x=expression,y=cluster,fill = type),scale="width",width=0.75)+
  theme_bw()+ ylab("")+xlab("")+
  facet_grid(.~symbol) +
  theme(plot.margin=unit(x=c(0.05,0,0,0),units="inches"),strip.background = element_blank()) +
  theme(strip.text.x = element_text(angle = 270,face = "italic")) +
  theme(axis.line.y = element_blank(),axis.text.y=element_text(color = "black",angle = 270)) +
  theme(axis.line.x = element_blank(),axis.ticks = element_blank()) +
  theme(axis.text.x = element_blank(),) +
  theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid=element_blank()) +
  scale_fill_manual(values = col) +ggtitle("SC9")
p2

ggsave("INsplit.violin.SC9.pdf",p2,width=12,height=6)

p3 <- ggplot(plot.data.final[plot.data.final$sample == "SC10",])+
  geom_violin(aes(x=expression,y=cluster,fill = type),scale="width",width=0.75)+
  theme_bw()+ ylab("")+xlab("")+
  facet_grid(.~symbol) +
  theme(plot.margin=unit(x=c(0.05,0,0,0),units="inches"),strip.background = element_blank()) +
  theme(strip.text.x = element_text(angle = 270,face = "italic")) +
  theme(axis.line.y = element_blank(),axis.text.y=element_text(color = "black",angle = 270)) +
  theme(axis.line.x = element_blank(),axis.ticks = element_blank()) +
  theme(axis.text.x = element_blank(),) +
  theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid=element_blank()) +
  scale_fill_manual(values = col) +ggtitle("SC10")
p3

ggsave("INsplit.violin.SC10.pdf",p3,width=12,height=6)

p4 <- ggplot(plot.data.final[plot.data.final$sample == "SC12",])+
  geom_violin(aes(x=expression,y=cluster,fill = type),scale="width",width=0.75)+
  theme_bw()+ ylab("")+xlab("")+
  facet_grid(.~symbol) +
  theme(plot.margin=unit(x=c(0.05,0,0,0),units="inches"),strip.background = element_blank()) +
  theme(strip.text.x = element_text(angle = 270,face = "italic")) +
  theme(axis.line.y = element_blank(),axis.text.y=element_text(color = "black",angle = 270)) +
  theme(axis.line.x = element_blank(),axis.ticks = element_blank()) +
  theme(axis.text.x = element_blank(),) +
  theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid=element_blank()) +
  scale_fill_manual(values = col) +ggtitle("SC12")
p4

ggsave("INsplit.violin.SC12.pdf",p4,width=12,height=6)

p.all <- p/p1/p2/p3/p4 +plot_layout(ncol = 1, nrow = 5)

ggsave("INsplit.violin.sample.pdf",p.all,width=12,height=30)

##### 20000515 figure 5 split sample
#C
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

#fetch
marker <- read.table("20210926.marker1.txt", header = F, sep = "\t")

#SC8
subset <- subset(project.insplit, Sample == "SC8")
plotdata <- fetch_dotplot_data(subset, features = as.vector(unique(marker$V1)), 
                               dot.scale = 6) 

head(plotdata)

plotdata$id <- factor(plotdata$id, levels = rev(c(20,16,17,8,21,2,6,22,18,1,28,12,19)))
plotdata$na.avg.exp <- plotdata$avg.exp
plotdata$type[plotdata$id %in% c(20,16,17,8,21,2,6)] <- "dorsal"
plotdata$type[plotdata$id %in% c(18,1,28,12,19)] <- "ventral"
plotdata$type[plotdata$id %in% c(22)] <- "mix"

col <- c("#FA9A99", "#1F78B4", "#B1DF89")

p <- ggplot(plotdata) + geom_point(aes(features.plot, id, fill= type, size = avg.exp),shape = 21, color = "grey") + theme_bw() +
  scale_fill_manual(values = col, name = "Celltype") + RotatedAxis() + scale_size_area(name = "Average expression")+
  xlab(NULL) + ylab("Subcluster") + theme(axis.text = element_text(color = "black"))
p
ggsave("20210926.s4hV6Sub.INtype.DotPlot.SC8.pdf",p,width=8,height=4)

#SC9
subset <- subset(project.insplit, Sample == "SC9")
plotdata <- fetch_dotplot_data(subset, features = as.vector(unique(marker$V1)), 
                               dot.scale = 6) 

head(plotdata)

plotdata$id <- factor(plotdata$id, levels = rev(c(20,16,17,8,21,2,6,22,18,1,28,12,19)))
plotdata$na.avg.exp <- plotdata$avg.exp
plotdata$type[plotdata$id %in% c(20,16,17,8,21,2,6)] <- "dorsal"
plotdata$type[plotdata$id %in% c(18,1,28,12,19)] <- "ventral"
plotdata$type[plotdata$id %in% c(22)] <- "mix"

col <- c("#FA9A99", "#1F78B4", "#B1DF89")

p <- ggplot(plotdata) + geom_point(aes(features.plot, id, fill= type, size = avg.exp),shape = 21, color = "grey") + theme_bw() +
  scale_fill_manual(values = col, name = "Celltype") + RotatedAxis() + scale_size_area(name = "Average expression")+
  xlab(NULL) + ylab("Subcluster") + theme(axis.text = element_text(color = "black"))
p
ggsave("20210926.s4hV6Sub.INtype.DotPlot.SC9.pdf",p,width=8,height=4)

#SC10
subset <- subset(project.insplit, Sample == "SC10")
plotdata <- fetch_dotplot_data(subset, features = as.vector(unique(marker$V1)), 
                               dot.scale = 6) 

head(plotdata)

plotdata$id <- factor(plotdata$id, levels = rev(c(20,16,17,8,21,2,6,22,18,1,28,12,19)))
plotdata$na.avg.exp <- plotdata$avg.exp
plotdata$type[plotdata$id %in% c(20,16,17,8,21,2,6)] <- "dorsal"
plotdata$type[plotdata$id %in% c(18,1,28,12,19)] <- "ventral"
plotdata$type[plotdata$id %in% c(22)] <- "mix"

col <- c("#FA9A99", "#1F78B4", "#B1DF89")

p <- ggplot(plotdata) + geom_point(aes(features.plot, id, fill= type, size = avg.exp),shape = 21, color = "grey") + theme_bw() +
  scale_fill_manual(values = col, name = "Celltype") + RotatedAxis() + scale_size_area(name = "Average expression")+
  xlab(NULL) + ylab("Subcluster") + theme(axis.text = element_text(color = "black"))
p
ggsave("20210926.s4hV6Sub.INtype.DotPlot.SC10.pdf",p,width=8,height=4)

#SC12
subset <- subset(project.insplit, Sample == "SC12")
plotdata <- fetch_dotplot_data(subset, features = as.vector(unique(marker$V1)), 
                               dot.scale = 6) 

head(plotdata)

plotdata$id <- factor(plotdata$id, levels = rev(c(20,16,17,8,21,2,6,22,18,1,28,12,19)))
plotdata$na.avg.exp <- plotdata$avg.exp
plotdata$type[plotdata$id %in% c(20,16,17,8,21,2,6)] <- "dorsal"
plotdata$type[plotdata$id %in% c(18,1,28,12,19)] <- "ventral"
plotdata$type[plotdata$id %in% c(22)] <- "mix"

col <- c("#FA9A99", "#1F78B4", "#B1DF89")

p <- ggplot(plotdata) + geom_point(aes(features.plot, id, fill= type, size = avg.exp),shape = 21, color = "grey") + theme_bw() +
  scale_fill_manual(values = col, name = "Celltype") + RotatedAxis() + scale_size_area(name = "Average expression")+
  xlab(NULL) + ylab("Subcluster") + theme(axis.text = element_text(color = "black"))
p
ggsave("20210926.s4hV6Sub.INtype.DotPlot.SC12.pdf",p,width=8,height=4)

# F
marker <- read.table("20220129.violin1.txt", header = F, sep = "\t")

project.insplit <- readRDS("D:/Projects/wj_sc_project2/RDS/project.insplit.RDS")

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

plot.data <- project.insplit@assays$RNA@scale.data[rownames(project.insplit@assays$RNA@scale.data) %in% rev(marker$V1),]
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
plot.data.final$sample <-project.insplit@meta.data$Sample
dim(plot.data.final)

plot.data.final <- reshape2::melt(plot.data.final,id.vars=c("barcode","cluster","sample"))
head(plot.data.final)
colnames(plot.data.final) <- c("barcode","cluster","sample","symbol","expression")
plot.data.final$symbol <- factor(plot.data.final$symbol, levels = rev(marker$V1))
plot.data.final$cluster <- factor(plot.data.final$cluster, 
                                  levels = rev(c(20,16,17,8,21,2,6,22,18,1,28,12,19)))
plot.data.final$type[plot.data.final$cluster %in% c(20,16,17,8,21,2,6)] <- "dorsal"
plot.data.final$type[plot.data.final$cluster %in% c(18,1,28,12,19)] <- "ventral"
plot.data.final$type[plot.data.final$cluster %in% c(22)] <- "mix"

col <- c("#FA9A99", "#1F78B4", "#B1DF89")

p <-ggplot(plot.data.final)+
  geom_violin(aes(x=expression,y=cluster,fill = type),scale="width",width=0.75)+
  theme_bw()+ ylab("")+xlab("")+
  facet_grid(.~symbol) +
  theme(plot.margin=unit(x=c(0.05,0,0,0),units="inches"),strip.background = element_blank()) +
  theme(strip.text.x = element_text(angle = 270,face = "italic")) +
  theme(axis.line.y = element_blank(),axis.text.y=element_text(color = "black",angle = 270)) +
  theme(axis.line.x = element_blank(),axis.ticks = element_blank()) +
  theme(axis.text.x = element_blank(),) +
  theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid=element_blank()) +
  scale_fill_manual(values = col) + ggtitle("All")
# theme(legend.position = "none")
p

ggsave("20220129.IN.miniviolin1.pdf",p,width=6,height=6)

p1 <- ggplot(plot.data.final[plot.data.final$sample == "SC8",])+
  geom_violin(aes(x=expression,y=cluster,fill = type),scale="width",width=0.75)+
  theme_bw()+ ylab("")+xlab("")+
  facet_grid(.~symbol) +
  theme(plot.margin=unit(x=c(0.05,0,0,0),units="inches"),strip.background = element_blank()) +
  theme(strip.text.x = element_text(angle = 270,face = "italic")) +
  theme(axis.line.y = element_blank(),axis.text.y=element_text(color = "black",angle = 270)) +
  theme(axis.line.x = element_blank(),axis.ticks = element_blank()) +
  theme(axis.text.x = element_blank(),) +
  theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid=element_blank()) +
  scale_fill_manual(values = col) +ggtitle("SC8")
p1

ggsave("20220129.IN.miniviolin1.SC8.pdf",p1,width=6,height=6)

p2 <- ggplot(plot.data.final[plot.data.final$sample == "SC9",])+
  geom_violin(aes(x=expression,y=cluster,fill = type),scale="width",width=0.75)+
  theme_bw()+ ylab("")+xlab("")+
  facet_grid(.~symbol) +
  theme(plot.margin=unit(x=c(0.05,0,0,0),units="inches"),strip.background = element_blank()) +
  theme(strip.text.x = element_text(angle = 270,face = "italic")) +
  theme(axis.line.y = element_blank(),axis.text.y=element_text(color = "black",angle = 270)) +
  theme(axis.line.x = element_blank(),axis.ticks = element_blank()) +
  theme(axis.text.x = element_blank(),) +
  theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid=element_blank()) +
  scale_fill_manual(values = col) +ggtitle("SC9")
p2

ggsave("20220129.IN.miniviolin1.SC9.pdf",p2,width=6,height=6)

p3 <- ggplot(plot.data.final[plot.data.final$sample == "SC10",])+
  geom_violin(aes(x=expression,y=cluster,fill = type),scale="width",width=0.75)+
  theme_bw()+ ylab("")+xlab("")+
  facet_grid(.~symbol) +
  theme(plot.margin=unit(x=c(0.05,0,0,0),units="inches"),strip.background = element_blank()) +
  theme(strip.text.x = element_text(angle = 270,face = "italic")) +
  theme(axis.line.y = element_blank(),axis.text.y=element_text(color = "black",angle = 270)) +
  theme(axis.line.x = element_blank(),axis.ticks = element_blank()) +
  theme(axis.text.x = element_blank(),) +
  theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid=element_blank()) +
  scale_fill_manual(values = col) +ggtitle("SC10")
p3

ggsave("20220129.IN.miniviolin1.SC10.pdf",p3,width=6,height=6)

p4 <- ggplot(plot.data.final[plot.data.final$sample == "SC12",])+
  geom_violin(aes(x=expression,y=cluster,fill = type),scale="width",width=0.75)+
  theme_bw()+ ylab("")+xlab("")+
  facet_grid(.~symbol) +
  theme(plot.margin=unit(x=c(0.05,0,0,0),units="inches"),strip.background = element_blank()) +
  theme(strip.text.x = element_text(angle = 270,face = "italic")) +
  theme(axis.line.y = element_blank(),axis.text.y=element_text(color = "black",angle = 270)) +
  theme(axis.line.x = element_blank(),axis.ticks = element_blank()) +
  theme(axis.text.x = element_blank(),) +
  theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid=element_blank()) +
  scale_fill_manual(values = col) +ggtitle("SC12")
p4

ggsave("20220129.IN.miniviolin1.SC12.pdf",p4,width=6,height=6)

p.all <- p/p1/p2/p3/p4 +plot_layout(ncol = 1, nrow = 5)

ggsave("20220129.IN.miniviolin1.sample.pdf",p.all,width=6,height=30)

# G
marker <- read.table("20220129.violin2.txt", header = F, sep = "\t")

project.insplit <- readRDS("D:/Projects/wj_sc_project2/RDS/project.insplit.RDS")

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

plot.data <- project.insplit@assays$RNA@scale.data[rownames(project.insplit@assays$RNA@scale.data) %in% rev(marker$V1),]
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
plot.data.final$sample <-project.insplit@meta.data$Sample
dim(plot.data.final)

plot.data.final <- reshape2::melt(plot.data.final,id.vars=c("barcode","cluster","sample"))
head(plot.data.final)
colnames(plot.data.final) <- c("barcode","cluster","sample","symbol","expression")
plot.data.final$symbol <- factor(plot.data.final$symbol, levels = rev(marker$V1))
plot.data.final$cluster <- factor(plot.data.final$cluster, 
                                  levels = rev(c(20,16,17,8,21,2,6,22,18,1,28,12,19)))
plot.data.final$type[plot.data.final$cluster %in% c(20,16,17,8,21,2,6)] <- "dorsal"
plot.data.final$type[plot.data.final$cluster %in% c(18,1,28,12,19)] <- "ventral"
plot.data.final$type[plot.data.final$cluster %in% c(22)] <- "mix"

col <- c("#FA9A99", "#1F78B4", "#B1DF89")

p <-ggplot(plot.data.final)+
  geom_violin(aes(x=expression,y=cluster,fill = type),scale="width",width=0.75)+
  theme_bw()+ ylab("")+xlab("")+
  facet_grid(.~symbol) +
  theme(plot.margin=unit(x=c(0.05,0,0,0),units="inches"),strip.background = element_blank()) +
  theme(strip.text.x = element_text(angle = 270,face = "italic")) +
  theme(axis.line.y = element_blank(),axis.text.y=element_text(color = "black",angle = 270)) +
  theme(axis.line.x = element_blank(),axis.ticks = element_blank()) +
  theme(axis.text.x = element_blank(),) +
  theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid=element_blank()) +
  scale_fill_manual(values = col) + ggtitle("All")
# theme(legend.position = "none")
p

ggsave("20220129.IN.miniviolin2.pdf",p,width=6,height=6)

p1 <- ggplot(plot.data.final[plot.data.final$sample == "SC8",])+
  geom_violin(aes(x=expression,y=cluster,fill = type),scale="width",width=0.75)+
  theme_bw()+ ylab("")+xlab("")+
  facet_grid(.~symbol) +
  theme(plot.margin=unit(x=c(0.05,0,0,0),units="inches"),strip.background = element_blank()) +
  theme(strip.text.x = element_text(angle = 270,face = "italic")) +
  theme(axis.line.y = element_blank(),axis.text.y=element_text(color = "black",angle = 270)) +
  theme(axis.line.x = element_blank(),axis.ticks = element_blank()) +
  theme(axis.text.x = element_blank(),) +
  theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid=element_blank()) +
  scale_fill_manual(values = col) +ggtitle("SC8")
p1

ggsave("20220129.IN.miniviolin2.SC8.pdf",p1,width=6,height=6)

p2 <- ggplot(plot.data.final[plot.data.final$sample == "SC9",])+
  geom_violin(aes(x=expression,y=cluster,fill = type),scale="width",width=0.75)+
  theme_bw()+ ylab("")+xlab("")+
  facet_grid(.~symbol) +
  theme(plot.margin=unit(x=c(0.05,0,0,0),units="inches"),strip.background = element_blank()) +
  theme(strip.text.x = element_text(angle = 270,face = "italic")) +
  theme(axis.line.y = element_blank(),axis.text.y=element_text(color = "black",angle = 270)) +
  theme(axis.line.x = element_blank(),axis.ticks = element_blank()) +
  theme(axis.text.x = element_blank(),) +
  theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid=element_blank()) +
  scale_fill_manual(values = col) +ggtitle("SC9")
p2

ggsave("20220129.IN.miniviolin2.SC9.pdf",p2,width=6,height=6)

p3 <- ggplot(plot.data.final[plot.data.final$sample == "SC10",])+
  geom_violin(aes(x=expression,y=cluster,fill = type),scale="width",width=0.75)+
  theme_bw()+ ylab("")+xlab("")+
  facet_grid(.~symbol) +
  theme(plot.margin=unit(x=c(0.05,0,0,0),units="inches"),strip.background = element_blank()) +
  theme(strip.text.x = element_text(angle = 270,face = "italic")) +
  theme(axis.line.y = element_blank(),axis.text.y=element_text(color = "black",angle = 270)) +
  theme(axis.line.x = element_blank(),axis.ticks = element_blank()) +
  theme(axis.text.x = element_blank(),) +
  theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid=element_blank()) +
  scale_fill_manual(values = col) +ggtitle("SC10")
p3

ggsave("20220129.IN.miniviolin2.SC10.pdf",p3,width=6,height=6)

p4 <- ggplot(plot.data.final[plot.data.final$sample == "SC12",])+
  geom_violin(aes(x=expression,y=cluster,fill = type),scale="width",width=0.75)+
  theme_bw()+ ylab("")+xlab("")+
  facet_grid(.~symbol) +
  theme(plot.margin=unit(x=c(0.05,0,0,0),units="inches"),strip.background = element_blank()) +
  theme(strip.text.x = element_text(angle = 270,face = "italic")) +
  theme(axis.line.y = element_blank(),axis.text.y=element_text(color = "black",angle = 270)) +
  theme(axis.line.x = element_blank(),axis.ticks = element_blank()) +
  theme(axis.text.x = element_blank(),) +
  theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid=element_blank()) +
  scale_fill_manual(values = col) +ggtitle("SC12")
p4

ggsave("20220129.IN.miniviolin2.SC12.pdf",p4,width=6,height=6)

p.all <- p/p1/p2/p3/p4 +plot_layout(ncol = 1, nrow = 5)

ggsave("20220129.IN.miniviolin2.sample.pdf",p.all,width=6,height=30)

###  pgsplit cluster heatmap
library(pheatmap)

project.pgsplit <- readRDS("D:/Projects/wj_sc_project2/RDS/project.pgsplit.RDS")
# cluster mean heatmap
expr <- as.data.frame(t(project.pgsplit@assays$RNA@scale.data))
expr[1:5,1:5]
expr$cluster <- project.pgsplit@active.ident[match(names(project.pgsplit@active.ident),rownames(expr))]
attach(expr)
cluster_mean.pgsplit <-  as.data.frame(aggregate(expr, list(cluster), mean))
rownames(cluster_mean.pgsplit) <- cluster_mean.pgsplit$Group.1
cluster_mean.pgsplit[1:5,1:5]
cluster_mean.pgsplit <- cluster_mean.pgsplit[,-1]
cluster_mean.pgsplit[,1:5]
plot <- t(cluster_mean.pgsplit)
plot[1:5,1:5]

data <- readxl::read_xlsx("20220516.xlsx")
pdf("20220516.PGsplit.heatmap.pdf", height=10)
for (tf in unique(data$TF)){
  gene <- data$gene[data$TF == tf]
  plotdata <- plot[colnames(cluster_mean.pgsplit) %in% gene,]
  #plotdata <- t(plotdata)
  
  colnames(plotdata) <- c(1:8)
  plotdata <- plotdata[match(gene,rownames(plotdata)),match(c(5,7,4,3,8,2,1,6),colnames(plotdata))]
  pheatmap(plotdata, cluster_cols = F, cluster_rows = T, scale = "row", cellheight = 10,cellwidth = 20)
}

dev.off()

DimPlot(project.pgsplit)

### pgsplit alluvial
library(ggalluvial)
library(magrittr)
library(dplyr)

data <- readxl::read_xlsx("20220516.xlsx")
data <- data[,c(1,2,13)]

#data$group <- factor(data$group, levels = c("pIN1","pMN2","SVZ3","SVZ4",
#                                                 "VZ5","pIN6","VZ7","SVZ8"))

data <- arrange(data, cluster)

alldata <- group_by(data,cluster,TF,gene) %>% summarise(., count = n())
#alldata$cluster <- factor(alldata$cluster, levels = c(3,4,5,7,8,2,1,6))

data_long <- to_lodes_form(data.frame(alldata),
                           key = "Demographic",
                           axes = 2:3)
data_long
#data_long$cluster <- factor(data_long$cluster, levels = c(3,4,5,7,8))

level <- unique(c(unique(data$gene),c("ASCL1","HES6","LHX1","RAD21")))
level <- unique(c(data$gene,"ASCL1","HES6","LHX1","RAD21"))
data_long$stratum <- factor(data_long$stratum, levels = level)
data_long$cluster <- factor(data_long$cluster)
#data_long$group <- factor(data_long$group, levels = c("pIN1","pMN2","SVZ3","SVZ4",
#                                                      "VZ5","pIN6","VZ7","SVZ8"))

color <- c("#94cdf6",
         "#fb9a99",
         "#a3ca3e",
         "#6f892a",
         "#ff7f01",
         "#77a5ca",
         "#fdc884",
         "#ccfb4d")

ggplot(data = data_long,
       aes(x = Demographic, stratum = stratum, alluvium = alluvium,
           y = 1, label = stratum)) +
  geom_alluvium(aes(fill=cluster))+
  geom_stratum(aes(fill=cluster)) + geom_text(stat = "stratum") +
  theme_minimal() + scale_fill_manual(values = color) +
  theme(axis.title = element_blank(), axis.text = element_blank(), 
        panel.grid = element_blank())

ggsave("20220516.alluvial.pdf",width = 7,height=12)

## 20200519 cds.monpn violin
cds.monpn <- readRDS("../RDS/cds.monpn.RDS")

plot_cells(cds.monpn)

colData(cds.monpn)$clusters <- factor(colData(cds.monpn)$clusters,
                                      levels = c(10,3,13,2,12,6,4,11,1,9,8,7,5))
# sheet1
data <- readxl::read_xlsx("20220519NRURONVIOLIN.xlsx", sheet = 1)

cds_subset <- cds.monpn[rowData(cds.monpn)$gene_short_name %in% data$sheet1[c(1)],]
p <- plot_genes_violin(cds_subset, group_cells_by="clusters") + facet_wrap(~Sample, ncol = 4) +
  ggtitle(data$sheet1[c(1)])

for (i in c(2:54)){
  cds_subset <- cds.monpn[rowData(cds.monpn)$gene_short_name %in% data$sheet1[i],]
  p1 <- plot_genes_violin(cds_subset, group_cells_by="clusters") + facet_wrap(~Sample, ncol = 4) +
    ggtitle(data$sheet1[i])
  #print(data$sheet1[i])
  #print(p1)
  p <- p/p1
}
ggsave("20220519.neuronviolin1.pdf", p, width = 10, height = 100, limitsize = F)

# sheet2
data <- readxl::read_xlsx("20220519NRURONVIOLIN.xlsx", sheet = 2)

cds_subset <- cds.monpn[rowData(cds.monpn)$gene_short_name %in% data$sheet2[c(1)],]
p <- plot_genes_violin(cds_subset, group_cells_by="clusters") + facet_wrap(~Sample, ncol = 4) +
  ggtitle(data$sheet2[c(1)])

for (i in c(2:78)){
  cds_subset <- cds.monpn[rowData(cds.monpn)$gene_short_name %in% data$sheet2[i],]
  p1 <- plot_genes_violin(cds_subset, group_cells_by="clusters") + facet_wrap(~Sample, ncol = 4) +
    ggtitle(data$sheet2[i])
  #print(data$sheet2[i])
  #print(p1)
  p <- p/p1
}
ggsave("20220519.neuronviolin2.pdf", p, width = 10, height = 150, limitsize = F)

### 20220529 monPN
library(monocle3)
library(ggplot2)

cds.monpn <- readRDS("../RDS/cds.monpn.RDS")
col <- c("#D9D9D9",
         "#377EB8",
         "#CAB4D6",
         "#D9D9D9",
         "#D9D9D9",
         "#D9D9D9",
         "#D9D9D9",
         "#D9D9D9",
         "#D9D9D9",
         "#764DA3",
         "#D9D9D9",
         "#B5D7E8",
         "#fcafff")

p<-plot_cells(cds.monpn, reduction_method = "UMAP", graph_label_size=5, show_trajectory_graph = T,
              color_cells_by = "cluster", 
              label_cell_groups = F, label_groups_by_cluster = F,
              label_branch_points = F, label_roots = F, label_leaves = F,
              label_principal_points = F) +
  scale_color_manual(values = col) + coord_fixed()
p
ggsave("20220529.Fig4A.pdf",p)

col <-  c("#6ab662",
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
       "#fcafff")
p<-plot_cells(cds.monpn, reduction_method = "UMAP", graph_label_size=5, show_trajectory_graph = F,
              color_cells_by = "cluster", 
              label_cell_groups = F, label_groups_by_cluster = F,
              label_branch_points = F, label_roots = F, label_leaves = F,
              label_principal_points = F) +
  scale_color_manual(values = col) + coord_fixed()
p
ggsave("20220529.Fig1E.pdf",p)


col <- c("#D9D9D9",
         "#377EB8",
         "#CAB4D6",
         "#D9D9D9",
         "#D9D9D9",
         "#D9D9D9",
         "#D9D9D9",
         "#D9D9D9",
         "#D9D9D9",
         "#764DA3",
         "#D9D9D9",
         "#B5D7E8",
         "#fcafff")

p<-plot_cells(cds.monpn, reduction_method = "UMAP", graph_label_size=5, show_trajectory_graph = F,
              color_cells_by = "cluster", 
              label_cell_groups = F, label_groups_by_cluster = F,
              label_branch_points = F, label_roots = F, label_leaves = F,
              label_principal_points = F) +
  scale_color_manual(values = col) + coord_fixed() + facet_wrap(~Sample)
p
ggsave("20220529.Fig4A.splitSample.pdf",p)

col <-  c("#6ab662",
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
          "#fcafff")
p<-plot_cells(cds.monpn, reduction_method = "UMAP", graph_label_size=5, show_trajectory_graph = F,
              color_cells_by = "cluster", 
              label_cell_groups = F, label_groups_by_cluster = F,
              label_branch_points = F, label_roots = F, label_leaves = F,
              label_principal_points = F) +
  scale_color_manual(values = col) + coord_fixed() + facet_wrap(~Sample)
p
ggsave("20220529.Fig1E.splitSample.pdf",p)


### 20220529 fig 2E
col <- c("#94cdf6",
         "#fb9a99",
         "#a3ca3e",
         "#6f892a",
         "#ff7f01",
         "#77a5ca",
         "#fdc884",
         "#ccfb4d")

project.pgsplit <- readRDS("D:/Projects/wj_sc_project2/RDS/project.pgsplit.RDS")
marker <- read.table("Fig2E.txt",sep = "\t", header = F)

plot.data <- project.pgsplit@assays$RNA@scale.data[rownames(project.pgsplit@assays$RNA@scale.data) %in% rev(marker$V1),]
dim(plot.data)
plot.data <- as.data.frame(plot.data)
plot.data[c(1:5),c(1:5)]

plot.data <- t(plot.data)
plot.data <- as.data.frame(plot.data)
plot.data[c(1:5),c(1:5)]
dim(plot.data)

##### merge cluster
length(project.pgsplit@active.ident)
plot.data.final <- plot.data

plot.data.final$cluster <- as.factor(project.pgsplit@active.ident)
plot.data.final$barcode <- rownames(plot.data.final)
plot.data.final$sample <-project.pgsplit@meta.data$Sample
dim(plot.data.final)

plot.data.final <- reshape2::melt(plot.data.final,id.vars=c("barcode","cluster","sample"))
head(plot.data.final)
colnames(plot.data.final) <- c("barcode","cluster","sample","symbol","expression")
plot.data.final$symbol <- factor(plot.data.final$symbol, levels = rev(marker$V1))
plot.data.final$cluster <- factor(plot.data.final$cluster, 
                                  levels = c(3,4,5,7,8,2,1,6))

plot.data.final2 <- plot.data.final[plot.data.final$cluster %in% c(3,4,5,7,8),]
plot.data.final2$cluster <- factor(plot.data.final2$cluster, 
                                   levels = rev(c(5,7,3,4,8)))

col5 <- col[c(1:8)][rev(c(5,7,3,4,8))]

p2 <-ggplot(plot.data.final2)+
  geom_violin(aes(x=expression,y=cluster,fill = cluster),scale="width",width=0.75)+
  theme_bw()+ ylab("")+xlab("")+
  facet_grid(.~symbol) +
  theme(plot.margin=unit(x=c(0.05,0,0,0),units="inches"),strip.background = element_blank()) +
  theme(strip.text.x = element_text(angle = 270,face = "italic")) +
  theme(axis.line.y = element_blank(),axis.text.y=element_text(color = "black",angle = 270)) +
  theme(axis.line.x = element_blank(),axis.ticks = element_blank(), legend.position = "none") +
  theme(axis.text.x = element_blank(),) +
  theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid=element_blank()) +
  scale_fill_manual(values = col5) 
p2

ggsave("20220529.Fig2E.pdf",p2,width = 24,height = 8, units = "cm" )

### GO
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(stringr)

#data <- read.table("20220211.pIN.GO.txt", header = T, sep = "\t")
data <- readxl::read_xlsx("20220524作图对应数据.xlsx", sheet = 5)
data$Ratio <- data$Count/1596
data <- data[order(data$Ratio),]
data$Description <- factor(data$Description, levels = data$Description)

ggplot(data) + geom_point(aes(Ratio, Description, fill=-log10(pvalue), size = Count), shape = 21) + theme_bw() +
  theme(axis.text = element_text(color = "black",size=10), axis.ticks.y = element_blank()) +
  scale_fill_gradientn(colors = brewer.pal(9, "YlGnBu")[2:5]) +
  ylab(NULL)+   scale_y_discrete(labels=(function(x) str_wrap(x, width=25))) +
  coord_cartesian(expand = TRUE)

ggsave("20220529.pMN.GO.pdf", width = 6, height = 4)

data <- readxl::read_xlsx("20220524作图对应数据.xlsx", sheet = 6)
data$Ratio <- data$Count/628
data <- data[order(data$Ratio),]
data$Description <- factor(data$Description, levels = data$Description)

ggplot(data) + geom_point(aes(Ratio, Description, fill=-log10(pvalue), size = Count), shape = 21) + theme_bw() +
  theme(axis.text = element_text(color = "black",size=10), axis.ticks.y = element_blank()) +
  scale_fill_gradientn(colors = brewer.pal(9, "YlGnBu")[2:5]) +
  ylab(NULL)+   scale_y_discrete(labels=(function(x) str_wrap(x, width=25))) +
  coord_cartesian(expand = TRUE)

ggsave("20220529.pIN.GO.pdf", width = 6, height = 4)

### 20220529 20210923.s4hV6Sub.INsplit.mergeCluster.pdf split sample
library(Seurat)
project.insplit <- readRDS("../RDS/project.insplit.RDS")

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


in.col <- c("#fb9a99",
            "#f781bf",
            "#ec56a3",
            "#984ea3",
            "#cab2d6",
            "#a4a1cd",
            "#6d419c",
            "#f7bb55",
            "#4daf4a",
            "#1b9e77",
            "#66a61e",
            "#b8e294",
            "#33a02c")

project.insplit@active.ident <- factor(project.insplit@active.ident, levels = c(20,16,17,8,21,2,6,22,18,1,28,12,19))

pdf("20220529.s4hV6Sub.INsplit.mergeCluster.pdf",width = 15,height = 10)
DimPlot(project.insplit, reduction = "pca" ,cols =in.col,split.by = "Sample",ncol = 2) +coord_fixed()
DimPlot(project.insplit, reduction = "harmony",cols =in.col,split.by = "Sample",ncol = 2) +coord_fixed()
DimPlot(project.insplit, reduction = "tsne" ,cols =in.col,split.by = "Sample",ncol = 2) +coord_fixed()
DimPlot(project.insplit, reduction = "umap" ,cols =in.col ,split.by = "Sample",ncol = 2) +coord_fixed()
DimPlot(project.insplit, reduction = "tsne" ,label=T ,cols =in.col,split.by = "Sample",ncol = 2) +coord_fixed()
DimPlot(project.insplit, reduction = "umap" ,label=T ,cols =in.col,split.by = "Sample",ncol = 2) +coord_fixed()
dev.off()

### 20220605 cds.monpn ratio
library(monocle3)
library(ggplot2)

cds.monpn <- readRDS("../RDS/cds.monpn.RDS")

data <- data.frame(colData(cds.monpn)$Sample,colData(cds.monpn)$clusters)
colnames(data) <- c("Sample","Cluster")

p1 <- ggplot(data) + geom_bar(aes(Cluster, fill=Sample),position = "fill") +
  theme_classic() + ylab("Relative percent") +
  scale_fill_manual(values = c("#FFDFEC","#A7BD8E","#FCBE6E","#7570B3"))

p2 <- ggplot(data[data$Sample %in% c("SC8","SC9","SC12"),]) + geom_bar(aes(Cluster, fill=Sample),position = "fill") +
  theme_classic() + ylab("Relative percent") +
  scale_fill_manual(values = c("#FFDFEC","#A7BD8E","#7570B3"))

ggsave("20220605.monpn.Cluster_Ratio1.pdf",p1)
ggsave("20220605.monpn.Cluster_Ratio2.pdf",p2)
### 20220605 cds.monpn violin
marker <- c("PAX6",
            "PAX3",
            "OLIG2",
            "NKX6-1",
            "NEUROG1",
            "NEUROG2",
            "DBX2",
            "GSX1",
            "GSX2",
            "MNX1",
            "BHLHE22",
            "LHX5")
#plot.data <- project.pgsplit@assays$RNA@scale.data[rownames(project.pgsplit@assays$RNA@scale.data) 
#                                                   %in% rev(marker),]
#violin1 normalized_counts
plot.data <- normalized_counts(cds.monpn)[rownames(normalized_counts(cds.monpn)) %in% rev(marker),]

dim(plot.data)
plot.data <- as.data.frame(as.matrix(plot.data))
plot.data[c(1:5),c(1:5)]

plot.data <- t(plot.data)
plot.data <- as.data.frame(plot.data)
plot.data[c(1:5),c(1:5)]
dim(plot.data)

##### merge cluster
dim(colData(cds.monpn))
plot.data.final <- plot.data

#plot.data.final$cluster <- as.factor(project.pgsplit@active.ident)
plot.data.final$cluster <- as.factor(colData(cds.monpn)$clusters)
plot.data.final$barcode <- rownames(plot.data.final)
#plot.data.final$sample <-project.pgsplit@meta.data$Sample
plot.data.final$sample <- colData(cds.monpn)$Sample
dim(plot.data.final)

plot.data.final <- reshape2::melt(plot.data.final,id.vars=c("barcode","cluster","sample"))
head(plot.data.final)
colnames(plot.data.final) <- c("barcode","cluster","sample","symbol","expression")
plot.data.final$symbol <- factor(plot.data.final$symbol, levels = rev(marker))
plot.data.final$cluster <- factor(plot.data.final$cluster, 
                                  levels = rev(c(1:13)))

#plot.data.final2 <- plot.data.final[plot.data.final$cluster %in% c(8,9,12),]
#plot.data.final2$cluster <- factor(plot.data.final2$cluster, 
#                                   levels = rev(c(8,9,12)))
col <-  c("#6ab662",
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
          "#fcafff")
#col5 <- rev(col[c(8,9,12)])

p2 <-ggplot(plot.data.final)+
  geom_violin(aes(x=expression,y=cluster,fill = cluster),scale="width",width=0.75)+
  theme_bw()+ ylab("")+xlab("")+
  facet_grid(.~symbol) +
  theme(plot.margin=unit(x=c(0.05,0,0,0),units="inches"),strip.background = element_blank()) +
  theme(strip.text.x = element_text(angle = 270,face = "italic")) +
  theme(axis.line.y = element_blank(),axis.text.y=element_text(color = "black",angle = 270)) +
  theme(axis.line.x = element_blank(),axis.ticks = element_blank(), legend.position = "none") +
  theme(axis.text.x = element_blank(),) +
  theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid=element_blank()) +
  scale_fill_manual(values = rev(col)) + scale_x_log10()
p2

ggsave("20220605.violin1.pdf",p2,width = 18,height = 8, units = "cm" )

cds_subset <- cds.monpn[rowData(cds.monpn)$gene_short_name %in% marker,]
rowData(cds_subset)$gene_short_name <- factor(rowData(cds_subset)$gene_short_name,
                                              levels = marker)

### use monocle3 violin
library(patchwork)
library(monocle3)

cds.monpn <- readRDS("../RDS/cds.monpn.RDS")

cds.monpn <- cds.monpn[,pData(cds.monpn)$Sample %in% c("SC8","SC9","SC12")]
cluster <- clusters(cds.monpn)
cds.monpn@colData$clusters <-factor(cluster, levels = c(6,10,3,13,2,12,4,11,1,9,8,5,7))

#p1 <- plot_genes_violin(cds_subset, group_cells_by="clusters") 
#  scale_fill_manual(values = col[c(6,10,3,13,2,12,4,11,1,6,8,5,7)])
#p1

# split sample

cds_subset <- cds.monpn[rowData(cds.monpn)$gene_short_name %in% marker[1],]
p <- plot_genes_violin(cds_subset, group_cells_by="clusters") + facet_wrap(~Sample, ncol = 4) +
  ggtitle(marker[1]) + scale_fill_manual(values = col[c(6,10,3,13,2,12,4,11,1,9,8,5,7)]) + 
  theme(axis.line.x = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank()) + xlab(NULL)
p
for (i in c(2:(length(marker)-1))){
  cds_subset <- cds.monpn[rowData(cds.monpn)$gene_short_name %in% marker[i],]
  p1 <- plot_genes_violin(cds_subset, group_cells_by="clusters") + facet_wrap(~Sample, ncol = 4) +
    ggtitle(marker[i]) + scale_fill_manual(values = col[c(6,10,3,13,2,12,4,11,1,9,8,5,7)]) +
    theme(axis.line.x = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank()) + xlab(NULL)
  print(marker[i])
  print(p1)
  p <- p/p1
}

cds_subset <- cds.monpn[rowData(cds.monpn)$gene_short_name %in% marker[length(marker)],]
p1 <- plot_genes_violin(cds_subset, group_cells_by="clusters") + facet_wrap(~Sample, ncol = 4) +
  ggtitle(marker[length(marker)]) + scale_fill_manual(values = col[c(6,10,3,13,2,12,4,11,1,9,8,5,7)])
print(marker[length(marker)])
print(p1)
p <- p/p1

ggsave("20220605.neuronviolin1.pdf", p, width = 10, height = 20, limitsize = F)

### 20220605 rainbow heatmap
library(ComplexHeatmap)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(circlize)
library(monocle3)
library(pheatmap)

genes  <- c("PAX6",
                     "PAX3",
                     "OLIG2",
                     "NKX6-1",
                     "NEUROG1",
                     "DBX2",
                     "GSX1",
                     "GSX2",
                     "NEUROG2",
                     "MNX1",
                     "BHLHE22",
                     "LHX5")

cds.monpn <- order_cells(cds.monpn, root_pr_nodes = c("Y_25","Y_49"))
pt.matrix <- exprs(cds.monpn)[match(genes,rownames(rowData(cds.monpn))),order(pseudotime(cds.monpn))]
#Can also use "normalized_counts" instead of "exprs" to use various normalization methods, for example:
#normalized_counts(cds, norm_method = "log")

pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes

df <- data.frame(Cluster = clusters(cds.monpn)[match(names(pseudotime(cds.monpn)[order(pseudotime(cds.monpn))]),names(clusters(cds.monpn)))], 
                 Pseudotime = pseudotime(cds.monpn)[order(pseudotime(cds.monpn))])

col <-  c("#6ab662",
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
          "#fcafff")
ClusterCol <- col
names(ClusterCol) <- sort(unique(clusters(cds.monpn)))
ClusterCol

max_pseudo <- ceiling(max(pseudotime(cds.monpn)))

ha <- HeatmapAnnotation(df = df,col = list(Cluster = ClusterCol, 
                                           Pseudotime = colorRamp2(seq(from=0,to=max_pseudo,length=11),rev(brewer.pal(11, "RdBu")))),
                        annotation_name_side = "left")


df1 <- data.frame(Cluster = clusters(cds.monpn)[match(names(pseudotime(cds.monpn)[order(pseudotime(cds.monpn))]),names(clusters(cds.monpn)))])
df2 <- data.frame(Pseudotime = pseudotime(cds.monpn)[order(pseudotime(cds.monpn))])

ha1 <- HeatmapAnnotation(df = df1, col = list(Cluster = ClusterCol), 
                         annotation_name_side = "left")
ha2 <- HeatmapAnnotation(df = df2, col = list(Pseudotime = colorRamp2(seq(from=0,to=max_pseudo,length=11),rev(brewer.pal(11, "RdBu")))), 
                         annotation_name_side = "left")
#K means with 6 groups
htkm <- Heatmap(
  as.matrix(pt.matrix),
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  km = 6,
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  top_annotation               = ha
)

htkm
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
  cluster_rows                 = FALSE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  top_annotation               = ha1,
  bottom_annotation            = ha2
)

hthc

pdf("20220605_s4hV6Sub_Monocle3Neuron.heatmap.pdf",width=12,height = 4)
  print(htkm)
  print(hthc)
dev.off()

#### 20220605 pseudotime curve
marker  <- c("PAX6",
            "PAX3",
            "OLIG2",
            "NKX6-1",
            "NEUROG1",
            "DBX2",
            "GSX1",
            "GSX2",
            "NEUROG2",
            "MNX1",
            "BHLHE22",
            "LHX5")

cds.monpn <- order_cells(cds.monpn, root_pr_nodes= c("Y_25","Y_49"))
plot_cells(cds.monpn, trajectory_graph_segment_size = 0.25,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)+ggtitle("Root Y_25 Y_49")

pt.matrix <- exprs(cds.monpn)[match(marker,rownames(rowData(cds.monpn))),]
#pt.matrix <- project.pgsplit@assays$RNA@scale.data[rownames(project.pgsplit@assays$RNA@scale.data) %in% marker$V1,]
pt.matrix <- as.data.frame(t(as.matrix(pt.matrix)))
dim(pt.matrix)
head(pt.matrix)

pt.matrix$pseudotime <- pseudotime(cds.monpn)
pt.matrix$cluster <- colData(cds.monpn)$clusters
data <- reshape2::melt(pt.matrix, id=c("pseudotime", "cluster"))
data$cluster <- factor(data$cluster, levels = c(1:13))
head(data)

data.plot <- data
#data.plot <- data[data$cluster %in% c(2,3,10,12,13),]

#col <- monpn.col[c(2,3,10,12,13)]

p<-ggplot(data.plot) + geom_point(aes(x = pseudotime, y = value, color = cluster),size = 0.5) +
  geom_smooth(aes(x=pseudotime, y=value),color="grey30") + scale_color_manual(values = col) +
  theme_bw()+  facet_wrap(.~variable,ncol=4,scales = "free_y") +  scale_y_log10() +
  theme(panel.grid = element_blank(), strip.background = element_blank()) +
  xlab("Pesudotime") + ylab("Expression/log10")
p

pdf("20220605.monpn.pseudotime_curve1.pdf",width=10,height = length(marker)/3)
  print(p)
dev.off()

#### 20220607 rainbow heatmap split
library(ComplexHeatmap)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(circlize)
library(monocle3)
library(pheatmap)

genes  <- c("PAX6",
            "PAX3",
            "OLIG2",
            "NKX6-1",
            "NEUROG1",
            "DBX2",
            "GSX1",
            "GSX2",
            "NEUROG2",
            "MNX1",
            "BHLHE22",
            "LHX5")


cds.monpn <- readRDS("../RDS/cds.monpn.RDS")
cds.monpn <- order_cells(cds.monpn, root_pr_nodes = c("Y_25","Y_49"))

all_cluster <- clusters(cds.monpn)[match(names(pseudotime(cds.monpn)[order(pseudotime(cds.monpn))]),names(clusters(cds.monpn)))]
all_pseudotime <- pseudotime(cds.monpn)[order(pseudotime(cds.monpn))]

# part1 6,10,3,13,2,12
pt.matrix <- exprs(cds.monpn)[match(genes,rownames(rowData(cds.monpn))),order(pseudotime(cds.monpn))]
pt.matrix <- pt.matrix[,colnames(pt.matrix) %in% rownames(colData(cds.monpn))[colData(cds.monpn)$clusters %in% c(6,10,3,13,2,12)]]

id <- rownames(colData(cds.monpn))[colData(cds.monpn)$clusters %in% c(6,10,3,13,2,12)]
#Can also use "normalized_counts" instead of "exprs" to use various normalization methods, for example:
#normalized_counts(cds, norm_method = "log")

pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes

df <- data.frame(Cluster = all_cluster[names(all_cluster) %in% id], 
                 Pseudotime = all_pseudotime[names(all_pseudotime)  %in% id])

col <-  c("#6ab662",
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
          "#fcafff")

col <- col[sort(c(6,10,3,13,2,12))]
ClusterCol <- col
names(ClusterCol) <- sort(c(6,10,3,13,2,12))
ClusterCol

max_pseudo <- ceiling(max(pseudotime(cds.monpn)))

ha <- HeatmapAnnotation(df = df,col = list(Cluster = ClusterCol, 
                                           Pseudotime = colorRamp2(seq(from=0,to=max_pseudo,length=11),rev(brewer.pal(11, "RdBu")))),
                        annotation_name_side = "left")


df1 <- data.frame(Cluster = all_cluster[names(all_cluster) %in% id])
df2 <- data.frame(Pseudotime = all_pseudotime[names(all_pseudotime)  %in% id])

ha1 <- HeatmapAnnotation(df = df1, col = list(Cluster = ClusterCol), 
                         annotation_name_side = "left")
ha2 <- HeatmapAnnotation(df = df2, col = list(Pseudotime = colorRamp2(seq(from=0,to=max_pseudo,length=11),rev(brewer.pal(11, "RdBu")))), 
                         annotation_name_side = "left")
#K means with 6 groups
htkm <- Heatmap(
  as.matrix(pt.matrix),
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  km = 6,
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  top_annotation               = ha
)

htkm
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
  cluster_rows                 = FALSE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  top_annotation               = ha1,
  bottom_annotation            = ha2
)

hthc

pdf("20220607_s4hV6Sub_Monocle3Neuron.heatmap1.pdf",width=6,height = 2)
  print(htkm)
  print(hthc)
dev.off()

# part1 6,4,11,1,9,8,5,7
pt.matrix <- exprs(cds.monpn)[match(genes,rownames(rowData(cds.monpn))),order(pseudotime(cds.monpn))]
pt.matrix <- pt.matrix[,colnames(pt.matrix) %in% rownames(colData(cds.monpn))[colData(cds.monpn)$clusters %in% c(6,4,11,1,9,8,5,7)]]

id <- rownames(colData(cds.monpn))[colData(cds.monpn)$clusters %in% c(6,4,11,1,9,8,5,7)]
#Can also use "normalized_counts" instead of "exprs" to use various normalization methods, for example:
#normalized_counts(cds, norm_method = "log")

pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes

df <- data.frame(Cluster = all_cluster[names(all_cluster) %in% id], 
                 Pseudotime = all_pseudotime[names(all_pseudotime)  %in% id])

col <-  c("#6ab662",
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
          "#fcafff")

col <- col[sort(c(6,4,11,1,9,8,5,7))]
ClusterCol <- col
names(ClusterCol) <- sort(c(6,4,11,1,9,8,5,7))
ClusterCol

max_pseudo <- ceiling(max(pseudotime(cds.monpn)))

ha <- HeatmapAnnotation(df = df,col = list(Cluster = ClusterCol, 
                                           Pseudotime = colorRamp2(seq(from=0,to=max_pseudo,length=11),rev(brewer.pal(11, "RdBu")))),
                        annotation_name_side = "left")


df1 <- data.frame(Cluster = all_cluster[names(all_cluster) %in% id])
df2 <- data.frame(Pseudotime = all_pseudotime[names(all_pseudotime)  %in% id])

ha1 <- HeatmapAnnotation(df = df1, col = list(Cluster = ClusterCol), 
                         annotation_name_side = "left")
ha2 <- HeatmapAnnotation(df = df2, col = list(Pseudotime = colorRamp2(seq(from=0,to=max_pseudo,length=11),rev(brewer.pal(11, "RdBu")))), 
                         annotation_name_side = "left")
#K means with 6 groups
htkm <- Heatmap(
  as.matrix(pt.matrix),
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  km = 6,
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  top_annotation               = ha
)

htkm
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
  cluster_rows                 = FALSE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  top_annotation               = ha1,
  bottom_annotation            = ha2
)

hthc

pdf("20220607_s4hV6Sub_Monocle3Neuron.heatmap2.pdf",width=6,height = 2)
  print(htkm)
  print(hthc)
dev.off()

#### 20220607 rainbow heatmap split, calculate z-score before split data
library(ComplexHeatmap)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(circlize)
library(monocle3)
library(pheatmap)

genes  <- c("PAX6",
            "PAX3",
            "OLIG2",
            "NKX6-1",
            "NEUROG1",
            "DBX2",
            "GSX1",
            "GSX2",
            "NEUROG2",
            "MNX1",
            "BHLHE22",
            "LHX5")


cds.monpn <- readRDS("../RDS/cds.monpn.RDS")
cds.monpn <- order_cells(cds.monpn, root_pr_nodes = c("Y_25","Y_49"))

all_cluster <- clusters(cds.monpn)[match(names(pseudotime(cds.monpn)[order(pseudotime(cds.monpn))]),names(clusters(cds.monpn)))]
all_pseudotime <- pseudotime(cds.monpn)[order(pseudotime(cds.monpn))]

# part1 6,10,3,13,2,12
pt.matrix <- exprs(cds.monpn)[match(genes,rownames(rowData(cds.monpn))),order(pseudotime(cds.monpn))]

all_id <- colnames(pt.matrix)
id <- rownames(colData(cds.monpn))[colData(cds.monpn)$clusters %in% c(6,10,3,13,2,12)]
#Can also use "normalized_counts" instead of "exprs" to use various normalization methods, for example:
#normalized_counts(cds, norm_method = "log")

pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes
pt.matrix <- pt.matrix[,all_id %in% rownames(colData(cds.monpn))[colData(cds.monpn)$clusters %in% c(6,10,3,13,2,12)]]


df <- data.frame(Cluster = all_cluster[names(all_cluster) %in% id], 
                 Pseudotime = all_pseudotime[names(all_pseudotime)  %in% id])

col <-  c("#6ab662",
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
          "#fcafff")

col <- col[sort(c(6,10,3,13,2,12))]
ClusterCol <- col
names(ClusterCol) <- sort(c(6,10,3,13,2,12))
ClusterCol

max_pseudo <- ceiling(max(pseudotime(cds.monpn)))

ha <- HeatmapAnnotation(df = df,col = list(Cluster = ClusterCol, 
                                           Pseudotime = colorRamp2(seq(from=0,to=max_pseudo,length=11),rev(brewer.pal(11, "RdBu")))),
                        annotation_name_side = "left")


df1 <- data.frame(Cluster = all_cluster[names(all_cluster) %in% id])
df2 <- data.frame(Pseudotime = all_pseudotime[names(all_pseudotime)  %in% id])

ha1 <- HeatmapAnnotation(df = df1, col = list(Cluster = ClusterCol), 
                         annotation_name_side = "left")
ha2 <- HeatmapAnnotation(df = df2, col = list(Pseudotime = colorRamp2(seq(from=0,to=max_pseudo,length=11),rev(brewer.pal(11, "RdBu")))), 
                         annotation_name_side = "left")
#K means with 6 groups
htkm <- Heatmap(
  as.matrix(pt.matrix),
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  km = 6,
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  top_annotation               = ha
)

htkm
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
  cluster_rows                 = FALSE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  top_annotation               = ha1,
  bottom_annotation            = ha2
)

hthc

pdf("20220607_s4hV6Sub_Monocle3Neuron.heatmap1.uni-zscore.pdf",width=6,height = 2)
print(htkm)
print(hthc)
dev.off()

# part2 6,4,11,1,9,8,5,7
pt.matrix <- exprs(cds.monpn)[match(genes,rownames(rowData(cds.monpn))),order(pseudotime(cds.monpn))]
#
all_id <- colnames(pt.matrix)
id <- rownames(colData(cds.monpn))[colData(cds.monpn)$clusters %in% c(6,4,11,1,9,8,5,7)]
#Can also use "normalized_counts" instead of "exprs" to use various normalization methods, for example:
#normalized_counts(cds, norm_method = "log")

pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes
pt.matrix <- pt.matrix[,all_id %in% rownames(colData(cds.monpn))[colData(cds.monpn)$clusters %in% c(6,4,11,1,9,8,5,7)]]

df <- data.frame(Cluster = all_cluster[names(all_cluster) %in% id], 
                 Pseudotime = all_pseudotime[names(all_pseudotime)  %in% id])

col <-  c("#6ab662",
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
          "#fcafff")

col <- col[sort(c(6,4,11,1,9,8,5,7))]
ClusterCol <- col
names(ClusterCol) <- sort(c(6,4,11,1,9,8,5,7))
ClusterCol

max_pseudo <- ceiling(max(pseudotime(cds.monpn)))

ha <- HeatmapAnnotation(df = df,col = list(Cluster = ClusterCol, 
                                           Pseudotime = colorRamp2(seq(from=0,to=max_pseudo,length=11),rev(brewer.pal(11, "RdBu")))),
                        annotation_name_side = "left")


df1 <- data.frame(Cluster = all_cluster[names(all_cluster) %in% id])
df2 <- data.frame(Pseudotime = all_pseudotime[names(all_pseudotime)  %in% id])

ha1 <- HeatmapAnnotation(df = df1, col = list(Cluster = ClusterCol), 
                         annotation_name_side = "left")
ha2 <- HeatmapAnnotation(df = df2, col = list(Pseudotime = colorRamp2(seq(from=0,to=max_pseudo,length=11),rev(brewer.pal(11, "RdBu")))), 
                         annotation_name_side = "left")
#K means with 6 groups
htkm <- Heatmap(
  as.matrix(pt.matrix),
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  km = 6,
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  top_annotation               = ha
)

htkm
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
  cluster_rows                 = FALSE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  top_annotation               = ha1,
  bottom_annotation            = ha2
)

hthc

pdf("20220607_s4hV6Sub_Monocle3Neuron.heatmap2.uni-zscore.pdf.pdf",width=6,height = 2)
print(htkm)
print(hthc)
dev.off()

##### 20220614 PGsplit rainbow
genes <- c("HOPX",
            "SLC1A3",
            "BCAN",
            "CLU",
            "MT2A",
            "CST3",
            "SPC25",
            "ASCL1",
            "HES6",
            "MKI67",
            "ASPM",
            "CKS1B",
            "UBE2C",
            "UBE2T",
            "CENPF",
            "CCNB1",
            "CCNB2")

library(ComplexHeatmap)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(circlize)
library(monocle3)
library(pheatmap)

cds.pgsplit <- readRDS("../RDS/cds.pgsplit.RDS")
plot_cells(cds.pgsplit)
plot_cells(cds.pgsplit,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = T,
           graph_label_size=1.5)

project.pgsplit <- readRDS("../RDS/project.pgsplit.RDS")
col <- c("#94cdf6",
         "#fb9a99",
         "#a3ca3e",
         "#6f892a",
         "#ff7f01",
         "#77a5ca",
         "#fdc884",
         "#ccfb4d")
Seurat::DimPlot(project.pgsplit, reduction = "umap" ,label=F , cols = col) +coord_fixed()

names(project.pgsplit@active.ident)[project.pgsplit@active.ident %in% c(1,4,5)]
project.pgsplit@meta.data$active.ident <- project.pgsplit@active.ident

way1 <- names(clusters(cds.pgsplit))[clusters(cds.pgsplit) %in% c(1,4,5)]
way1_SC8 <- names(clusters(cds.pgsplit))[project.pgsplit@meta.data$Sample == "SC8" & clusters(cds.pgsplit) %in% c(1,4,5)]
way1_SC9 <- names(clusters(cds.pgsplit))[project.pgsplit@meta.data$Sample == "SC9" & clusters(cds.pgsplit) %in% c(1,4,5)]
way1_SC10 <- names(clusters(cds.pgsplit))[project.pgsplit@meta.data$Sample == "SC10" & clusters(cds.pgsplit) %in% c(1,4,5)]
way1_SC12 <- names(clusters(cds.pgsplit))[project.pgsplit@meta.data$Sample == "SC12" & clusters(cds.pgsplit) %in% c(1,4,5)]

way2 <- names(clusters(cds.pgsplit))[clusters(cds.pgsplit) %in% c(2,3,6,7,8,9,10)]
way2_SC8 <- names(clusters(cds.pgsplit))[project.pgsplit@meta.data$Sample == "SC8" & clusters(cds.pgsplit) %in% c(2,3,6,7,8,9,10)]
way2_SC9 <- names(clusters(cds.pgsplit))[project.pgsplit@meta.data$Sample == "SC9" & clusters(cds.pgsplit) %in% c(2,3,6,7,8,9,10)]
way2_SC10 <- names(clusters(cds.pgsplit))[project.pgsplit@meta.data$Sample == "SC10" & clusters(cds.pgsplit) %in% c(2,3,6,7,8,9,10)]
way2_SC12 <- names(clusters(cds.pgsplit))[project.pgsplit@meta.data$Sample == "SC12" & clusters(cds.pgsplit) %in% c(2,3,6,7,8,9,10)]

##### All sample
all_cluster <- colData(cds.pgsplit)$Seurat_res0.2[match(names(pseudotime(cds.pgsplit)[order(pseudotime(cds.pgsplit))]),names(colData(cds.pgsplit)$Seurat_res0.2))]
all_pseudotime <- pseudotime(cds.pgsplit)[order(pseudotime(cds.pgsplit))]

# part1 1,4,5
pt.matrix <- exprs(cds.pgsplit)[match(genes,rownames(rowData(cds.pgsplit))),order(pseudotime(cds.pgsplit))]

all_id <- colnames(pt.matrix)
id <- way1
#Can also use "normalized_counts" instead of "exprs" to use various normalization methods, for example:
#normalized_counts(cds, norm_method = "log")

pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes
pt.matrix <- pt.matrix[,all_id %in% id]


df <- data.frame(Cluster = all_cluster[names(all_cluster) %in% id], 
                 Pseudotime = all_pseudotime[names(all_pseudotime)  %in% id])

col <- c("#94cdf6",
         "#fb9a99",
         "#a3ca3e",
         "#6f892a",
         "#ff7f01",
         "#77a5ca",
         "#fdc884",
         "#ccfb4d")

col <- col[sort(unique(df$Cluster))]
ClusterCol <- col
names(ClusterCol) <- sort(unique(df$Cluster))
ClusterCol

max_pseudo <- ceiling(max(pseudotime(cds.pgsplit)))

ha <- HeatmapAnnotation(df = df,col = list(Cluster = ClusterCol, 
                                           Pseudotime = colorRamp2(seq(from=0,to=max_pseudo,length=11),rev(brewer.pal(11, "RdBu")))),
                        annotation_name_side = "left")


df1 <- data.frame(Cluster = all_cluster[names(all_cluster) %in% id])
df2 <- data.frame(Pseudotime = all_pseudotime[names(all_pseudotime)  %in% id])

ha1 <- HeatmapAnnotation(df = df1, col = list(Cluster = ClusterCol), 
                         annotation_name_side = "left")
ha2 <- HeatmapAnnotation(df = df2, col = list(Pseudotime = colorRamp2(seq(from=0,to=max_pseudo,length=11),rev(brewer.pal(11, "RdBu")))), 
                         annotation_name_side = "left")
#K means with 6 groups
htkm <- Heatmap(
  as.matrix(pt.matrix),
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  km = 6,
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  top_annotation               = ha
)

htkm
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
  cluster_rows                 = FALSE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  top_annotation               = ha1,
  bottom_annotation            = ha2
)

hthc

pdf("20220614_s4hV6Sub_Monocle3PGsplit.heatmap1.uni-zscore.pdf",width=6,height = 2)
print(htkm)
print(hthc)
dev.off()

# part2 2,3,6,7,8,9,10
pt.matrix <- exprs(cds.pgsplit)[match(genes,rownames(rowData(cds.pgsplit))),order(pseudotime(cds.pgsplit))]

all_id <- colnames(pt.matrix)
id <- way2
#Can also use "normalized_counts" instead of "exprs" to use various normalization methods, for example:
#normalized_counts(cds, norm_method = "log")

pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes
pt.matrix <- pt.matrix[,all_id %in% id]


df <- data.frame(Cluster = all_cluster[names(all_cluster) %in% id], 
                 Pseudotime = all_pseudotime[names(all_pseudotime)  %in% id])

col <- c("#94cdf6",
         "#fb9a99",
         "#a3ca3e",
         "#6f892a",
         "#ff7f01",
         "#77a5ca",
         "#fdc884",
         "#ccfb4d")

col <- col[sort(unique(df$Cluster))]
ClusterCol <- col
names(ClusterCol) <- sort(unique(df$Cluster))
ClusterCol

max_pseudo <- ceiling(max(pseudotime(cds.pgsplit)))

ha <- HeatmapAnnotation(df = df,col = list(Cluster = ClusterCol, 
                                           Pseudotime = colorRamp2(seq(from=0,to=max_pseudo,length=11),rev(brewer.pal(11, "RdBu")))),
                        annotation_name_side = "left")


df1 <- data.frame(Cluster = all_cluster[names(all_cluster) %in% id])
df2 <- data.frame(Pseudotime = all_pseudotime[names(all_pseudotime)  %in% id])

ha1 <- HeatmapAnnotation(df = df1, col = list(Cluster = ClusterCol), 
                         annotation_name_side = "left")
ha2 <- HeatmapAnnotation(df = df2, col = list(Pseudotime = colorRamp2(seq(from=0,to=max_pseudo,length=11),rev(brewer.pal(11, "RdBu")))), 
                         annotation_name_side = "left")
#K means with 6 groups
htkm <- Heatmap(
  as.matrix(pt.matrix),
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  km = 6,
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  top_annotation               = ha
)

htkm
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
  cluster_rows                 = FALSE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  top_annotation               = ha1,
  bottom_annotation            = ha2
)

hthc

pdf("20220614_s4hV6Sub_Monocle3PGsplit.heatmap2.uni-zscore.pdf",width=6,height = 2)
  print(htkm)
  print(hthc)
dev.off()


##### SC8
all_cluster <- colData(cds.pgsplit)$Seurat_res0.2[match(names(pseudotime(cds.pgsplit)[order(pseudotime(cds.pgsplit))]),names(colData(cds.pgsplit)$Seurat_res0.2))]
all_pseudotime <- pseudotime(cds.pgsplit)[order(pseudotime(cds.pgsplit))]

# part1 1,4,5
pt.matrix <- exprs(cds.pgsplit)[match(genes,rownames(rowData(cds.pgsplit))),order(pseudotime(cds.pgsplit))]

all_id <- colnames(pt.matrix)
id <- way1_SC8
#Can also use "normalized_counts" instead of "exprs" to use various normalization methods, for example:
#normalized_counts(cds, norm_method = "log")

pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes
pt.matrix <- pt.matrix[,all_id %in% id]


df <- data.frame(Cluster = all_cluster[names(all_cluster) %in% id], 
                 Pseudotime = all_pseudotime[names(all_pseudotime)  %in% id])

col <- c("#94cdf6",
         "#fb9a99",
         "#a3ca3e",
         "#6f892a",
         "#ff7f01",
         "#77a5ca",
         "#fdc884",
         "#ccfb4d")

col <- col[sort(unique(df$Cluster))]
ClusterCol <- col
names(ClusterCol) <- sort(unique(df$Cluster))
ClusterCol

max_pseudo <- ceiling(max(pseudotime(cds.pgsplit)))

ha <- HeatmapAnnotation(df = df,col = list(Cluster = ClusterCol, 
                                           Pseudotime = colorRamp2(seq(from=0,to=max_pseudo,length=11),rev(brewer.pal(11, "RdBu")))),
                        annotation_name_side = "left")


df1 <- data.frame(Cluster = all_cluster[names(all_cluster) %in% id])
df2 <- data.frame(Pseudotime = all_pseudotime[names(all_pseudotime)  %in% id])

ha1 <- HeatmapAnnotation(df = df1, col = list(Cluster = ClusterCol), 
                         annotation_name_side = "left")
ha2 <- HeatmapAnnotation(df = df2, col = list(Pseudotime = colorRamp2(seq(from=0,to=max_pseudo,length=11),rev(brewer.pal(11, "RdBu")))), 
                         annotation_name_side = "left")
#K means with 6 groups
htkm <- Heatmap(
  as.matrix(pt.matrix),
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  km = 6,
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  top_annotation               = ha
)

htkm
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
  cluster_rows                 = FALSE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  top_annotation               = ha1,
  bottom_annotation            = ha2
)

hthc

pdf("20220614_s4hV6Sub_Monocle3PGsplit.heatmap1.SC8.uni-zscore.pdf",width=6,height = 2)
print(htkm)
print(hthc)
dev.off()

# part2 2,3,6,7,8,9,10
pt.matrix <- exprs(cds.pgsplit)[match(genes,rownames(rowData(cds.pgsplit))),order(pseudotime(cds.pgsplit))]

all_id <- colnames(pt.matrix)
id <- way2_SC8
#Can also use "normalized_counts" instead of "exprs" to use various normalization methods, for example:
#normalized_counts(cds, norm_method = "log")

pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes
pt.matrix <- pt.matrix[,all_id %in% id]


df <- data.frame(Cluster = all_cluster[names(all_cluster) %in% id], 
                 Pseudotime = all_pseudotime[names(all_pseudotime)  %in% id])

col <- c("#94cdf6",
         "#fb9a99",
         "#a3ca3e",
         "#6f892a",
         "#ff7f01",
         "#77a5ca",
         "#fdc884",
         "#ccfb4d")

col <- col[sort(unique(df$Cluster))]
ClusterCol <- col
names(ClusterCol) <- sort(unique(df$Cluster))
ClusterCol

max_pseudo <- ceiling(max(pseudotime(cds.pgsplit)))

ha <- HeatmapAnnotation(df = df,col = list(Cluster = ClusterCol, 
                                           Pseudotime = colorRamp2(seq(from=0,to=max_pseudo,length=11),rev(brewer.pal(11, "RdBu")))),
                        annotation_name_side = "left")


df1 <- data.frame(Cluster = all_cluster[names(all_cluster) %in% id])
df2 <- data.frame(Pseudotime = all_pseudotime[names(all_pseudotime)  %in% id])

ha1 <- HeatmapAnnotation(df = df1, col = list(Cluster = ClusterCol), 
                         annotation_name_side = "left")
ha2 <- HeatmapAnnotation(df = df2, col = list(Pseudotime = colorRamp2(seq(from=0,to=max_pseudo,length=11),rev(brewer.pal(11, "RdBu")))), 
                         annotation_name_side = "left")
#K means with 6 groups
htkm <- Heatmap(
  as.matrix(pt.matrix),
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  km = 6,
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  top_annotation               = ha
)

htkm
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
  cluster_rows                 = FALSE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  top_annotation               = ha1,
  bottom_annotation            = ha2
)

hthc

pdf("20220614_s4hV6Sub_Monocle3PGsplit.heatmap2.SC8.uni-zscore.pdf",width=6,height = 2)
print(htkm)
print(hthc)
dev.off()

##### SC9
all_cluster <- colData(cds.pgsplit)$Seurat_res0.2[match(names(pseudotime(cds.pgsplit)[order(pseudotime(cds.pgsplit))]),names(colData(cds.pgsplit)$Seurat_res0.2))]
all_pseudotime <- pseudotime(cds.pgsplit)[order(pseudotime(cds.pgsplit))]

# part1 1,4,5
pt.matrix <- exprs(cds.pgsplit)[match(genes,rownames(rowData(cds.pgsplit))),order(pseudotime(cds.pgsplit))]

all_id <- colnames(pt.matrix)
id <- way1_SC9
#Can also use "normalized_counts" instead of "exprs" to use various normalization methods, for example:
#normalized_counts(cds, norm_method = "log")

pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes
pt.matrix <- pt.matrix[,all_id %in% id]


df <- data.frame(Cluster = all_cluster[names(all_cluster) %in% id], 
                 Pseudotime = all_pseudotime[names(all_pseudotime)  %in% id])

col <- c("#94cdf6",
         "#fb9a99",
         "#a3ca3e",
         "#6f892a",
         "#ff7f01",
         "#77a5ca",
         "#fdc884",
         "#ccfb4d")

col <- col[sort(unique(df$Cluster))]
ClusterCol <- col
names(ClusterCol) <- sort(unique(df$Cluster))
ClusterCol

max_pseudo <- ceiling(max(pseudotime(cds.pgsplit)))

ha <- HeatmapAnnotation(df = df,col = list(Cluster = ClusterCol, 
                                           Pseudotime = colorRamp2(seq(from=0,to=max_pseudo,length=11),rev(brewer.pal(11, "RdBu")))),
                        annotation_name_side = "left")


df1 <- data.frame(Cluster = all_cluster[names(all_cluster) %in% id])
df2 <- data.frame(Pseudotime = all_pseudotime[names(all_pseudotime)  %in% id])

ha1 <- HeatmapAnnotation(df = df1, col = list(Cluster = ClusterCol), 
                         annotation_name_side = "left")
ha2 <- HeatmapAnnotation(df = df2, col = list(Pseudotime = colorRamp2(seq(from=0,to=max_pseudo,length=11),rev(brewer.pal(11, "RdBu")))), 
                         annotation_name_side = "left")
#K means with 6 groups
htkm <- Heatmap(
  as.matrix(pt.matrix),
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  km = 6,
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  top_annotation               = ha
)

htkm
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
  cluster_rows                 = FALSE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  top_annotation               = ha1,
  bottom_annotation            = ha2
)

hthc

pdf("20220614_s4hV6Sub_Monocle3PGsplit.heatmap1.SC9.uni-zscore.pdf",width=6,height = 2)
print(htkm)
print(hthc)
dev.off()

# part2 2,3,6,7,8,9,10
pt.matrix <- exprs(cds.pgsplit)[match(genes,rownames(rowData(cds.pgsplit))),order(pseudotime(cds.pgsplit))]

all_id <- colnames(pt.matrix)
id <- way2_SC9
#Can also use "normalized_counts" instead of "exprs" to use various normalization methods, for example:
#normalized_counts(cds, norm_method = "log")

pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes
pt.matrix <- pt.matrix[,all_id %in% id]


df <- data.frame(Cluster = all_cluster[names(all_cluster) %in% id], 
                 Pseudotime = all_pseudotime[names(all_pseudotime)  %in% id])

col <- c("#94cdf6",
         "#fb9a99",
         "#a3ca3e",
         "#6f892a",
         "#ff7f01",
         "#77a5ca",
         "#fdc884",
         "#ccfb4d")

col <- col[sort(unique(df$Cluster))]
ClusterCol <- col
names(ClusterCol) <- sort(unique(df$Cluster))
ClusterCol

max_pseudo <- ceiling(max(pseudotime(cds.pgsplit)))

ha <- HeatmapAnnotation(df = df,col = list(Cluster = ClusterCol, 
                                           Pseudotime = colorRamp2(seq(from=0,to=max_pseudo,length=11),rev(brewer.pal(11, "RdBu")))),
                        annotation_name_side = "left")


df1 <- data.frame(Cluster = all_cluster[names(all_cluster) %in% id])
df2 <- data.frame(Pseudotime = all_pseudotime[names(all_pseudotime)  %in% id])

ha1 <- HeatmapAnnotation(df = df1, col = list(Cluster = ClusterCol), 
                         annotation_name_side = "left")
ha2 <- HeatmapAnnotation(df = df2, col = list(Pseudotime = colorRamp2(seq(from=0,to=max_pseudo,length=11),rev(brewer.pal(11, "RdBu")))), 
                         annotation_name_side = "left")
#K means with 6 groups
htkm <- Heatmap(
  as.matrix(pt.matrix),
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  km = 6,
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  top_annotation               = ha
)

htkm
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
  cluster_rows                 = FALSE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  top_annotation               = ha1,
  bottom_annotation            = ha2
)

hthc

pdf("20220614_s4hV6Sub_Monocle3PGsplit.heatmap2.SC9.uni-zscore.pdf",width=6,height = 2)
print(htkm)
print(hthc)
dev.off()

##### SC10
all_cluster <- colData(cds.pgsplit)$Seurat_res0.2[match(names(pseudotime(cds.pgsplit)[order(pseudotime(cds.pgsplit))]),names(colData(cds.pgsplit)$Seurat_res0.2))]
all_pseudotime <- pseudotime(cds.pgsplit)[order(pseudotime(cds.pgsplit))]

# part1 1,4,5
pt.matrix <- exprs(cds.pgsplit)[match(genes,rownames(rowData(cds.pgsplit))),order(pseudotime(cds.pgsplit))]

all_id <- colnames(pt.matrix)
id <- way1_SC10
#Can also use "normalized_counts" instead of "exprs" to use various normalization methods, for example:
#normalized_counts(cds, norm_method = "log")

pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes
pt.matrix <- pt.matrix[,all_id %in% id]


df <- data.frame(Cluster = all_cluster[names(all_cluster) %in% id], 
                 Pseudotime = all_pseudotime[names(all_pseudotime)  %in% id])

col <- c("#94cdf6",
         "#fb9a99",
         "#a3ca3e",
         "#6f892a",
         "#ff7f01",
         "#77a5ca",
         "#fdc884",
         "#ccfb4d")

col <- col[sort(unique(df$Cluster))]
ClusterCol <- col
names(ClusterCol) <- sort(unique(df$Cluster))
ClusterCol

max_pseudo <- ceiling(max(pseudotime(cds.pgsplit)))

ha <- HeatmapAnnotation(df = df,col = list(Cluster = ClusterCol, 
                                           Pseudotime = colorRamp2(seq(from=0,to=max_pseudo,length=11),rev(brewer.pal(11, "RdBu")))),
                        annotation_name_side = "left")


df1 <- data.frame(Cluster = all_cluster[names(all_cluster) %in% id])
df2 <- data.frame(Pseudotime = all_pseudotime[names(all_pseudotime)  %in% id])

ha1 <- HeatmapAnnotation(df = df1, col = list(Cluster = ClusterCol), 
                         annotation_name_side = "left")
ha2 <- HeatmapAnnotation(df = df2, col = list(Pseudotime = colorRamp2(seq(from=0,to=max_pseudo,length=11),rev(brewer.pal(11, "RdBu")))), 
                         annotation_name_side = "left")
#K means with 6 groups
htkm <- Heatmap(
  as.matrix(pt.matrix),
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  km = 6,
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  top_annotation               = ha
)

htkm
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
  cluster_rows                 = FALSE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  top_annotation               = ha1,
  bottom_annotation            = ha2
)

hthc

pdf("20220614_s4hV6Sub_Monocle3PGsplit.heatmap1.SC10.uni-zscore.pdf",width=6,height = 2)
print(htkm)
print(hthc)
dev.off()

# part2 2,3,6,7,8,9,10
pt.matrix <- exprs(cds.pgsplit)[match(genes,rownames(rowData(cds.pgsplit))),order(pseudotime(cds.pgsplit))]

all_id <- colnames(pt.matrix)
id <- way2_SC10
#Can also use "normalized_counts" instead of "exprs" to use various normalization methods, for example:
#normalized_counts(cds, norm_method = "log")

pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes
pt.matrix <- pt.matrix[,all_id %in% id]


df <- data.frame(Cluster = all_cluster[names(all_cluster) %in% id], 
                 Pseudotime = all_pseudotime[names(all_pseudotime)  %in% id])

col <- c("#94cdf6",
         "#fb9a99",
         "#a3ca3e",
         "#6f892a",
         "#ff7f01",
         "#77a5ca",
         "#fdc884",
         "#ccfb4d")

col <- col[sort(unique(df$Cluster))]
ClusterCol <- col
names(ClusterCol) <- sort(unique(df$Cluster))
ClusterCol

max_pseudo <- ceiling(max(pseudotime(cds.pgsplit)))

ha <- HeatmapAnnotation(df = df,col = list(Cluster = ClusterCol, 
                                           Pseudotime = colorRamp2(seq(from=0,to=max_pseudo,length=11),rev(brewer.pal(11, "RdBu")))),
                        annotation_name_side = "left")


df1 <- data.frame(Cluster = all_cluster[names(all_cluster) %in% id])
df2 <- data.frame(Pseudotime = all_pseudotime[names(all_pseudotime)  %in% id])

ha1 <- HeatmapAnnotation(df = df1, col = list(Cluster = ClusterCol), 
                         annotation_name_side = "left")
ha2 <- HeatmapAnnotation(df = df2, col = list(Pseudotime = colorRamp2(seq(from=0,to=max_pseudo,length=11),rev(brewer.pal(11, "RdBu")))), 
                         annotation_name_side = "left")
#K means with 6 groups
htkm <- Heatmap(
  as.matrix(pt.matrix),
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  km = 6,
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  top_annotation               = ha
)

htkm
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
  cluster_rows                 = FALSE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  top_annotation               = ha1,
  bottom_annotation            = ha2
)

hthc

pdf("20220614_s4hV6Sub_Monocle3PGsplit.heatmap2.SC10.uni-zscore.pdf",width=6,height = 2)
print(htkm)
print(hthc)
dev.off()

##### SC12
all_cluster <- colData(cds.pgsplit)$Seurat_res0.2[match(names(pseudotime(cds.pgsplit)[order(pseudotime(cds.pgsplit))]),names(colData(cds.pgsplit)$Seurat_res0.2))]
all_pseudotime <- pseudotime(cds.pgsplit)[order(pseudotime(cds.pgsplit))]

# part1 1,4,5
pt.matrix <- exprs(cds.pgsplit)[match(genes,rownames(rowData(cds.pgsplit))),order(pseudotime(cds.pgsplit))]

all_id <- colnames(pt.matrix)
id <- way1_SC12
#Can also use "normalized_counts" instead of "exprs" to use various normalization methods, for example:
#normalized_counts(cds, norm_method = "log")

pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes
pt.matrix <- pt.matrix[,all_id %in% id]


df <- data.frame(Cluster = all_cluster[names(all_cluster) %in% id], 
                 Pseudotime = all_pseudotime[names(all_pseudotime)  %in% id])

col <- c("#94cdf6",
         "#fb9a99",
         "#a3ca3e",
         "#6f892a",
         "#ff7f01",
         "#77a5ca",
         "#fdc884",
         "#ccfb4d")

col <- col[sort(unique(df$Cluster))]
ClusterCol <- col
names(ClusterCol) <- sort(unique(df$Cluster))
ClusterCol

max_pseudo <- ceiling(max(pseudotime(cds.pgsplit)))

ha <- HeatmapAnnotation(df = df,col = list(Cluster = ClusterCol, 
                                           Pseudotime = colorRamp2(seq(from=0,to=max_pseudo,length=11),rev(brewer.pal(11, "RdBu")))),
                        annotation_name_side = "left")


df1 <- data.frame(Cluster = all_cluster[names(all_cluster) %in% id])
df2 <- data.frame(Pseudotime = all_pseudotime[names(all_pseudotime)  %in% id])

ha1 <- HeatmapAnnotation(df = df1, col = list(Cluster = ClusterCol), 
                         annotation_name_side = "left")
ha2 <- HeatmapAnnotation(df = df2, col = list(Pseudotime = colorRamp2(seq(from=0,to=max_pseudo,length=11),rev(brewer.pal(11, "RdBu")))), 
                         annotation_name_side = "left")
#K means with 6 groups
htkm <- Heatmap(
  as.matrix(pt.matrix),
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  km = 6,
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  top_annotation               = ha
)

htkm
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
  cluster_rows                 = FALSE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  top_annotation               = ha1,
  bottom_annotation            = ha2
)

hthc

pdf("20220614_s4hV6Sub_Monocle3PGsplit.heatmap1.SC12.uni-zscore.pdf",width=6,height = 2)
print(htkm)
print(hthc)
dev.off()

# part2 2,3,6,7,8,9,10
pt.matrix <- exprs(cds.pgsplit)[match(genes,rownames(rowData(cds.pgsplit))),order(pseudotime(cds.pgsplit))]

all_id <- colnames(pt.matrix)
id <- way2_SC12
#Can also use "normalized_counts" instead of "exprs" to use various normalization methods, for example:
#normalized_counts(cds, norm_method = "log")

pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes
pt.matrix <- pt.matrix[,all_id %in% id]


df <- data.frame(Cluster = all_cluster[names(all_cluster) %in% id], 
                 Pseudotime = all_pseudotime[names(all_pseudotime)  %in% id])

col <- c("#94cdf6",
         "#fb9a99",
         "#a3ca3e",
         "#6f892a",
         "#ff7f01",
         "#77a5ca",
         "#fdc884",
         "#ccfb4d")

col <- col[sort(unique(df$Cluster))]
ClusterCol <- col
names(ClusterCol) <- sort(unique(df$Cluster))
ClusterCol

max_pseudo <- ceiling(max(pseudotime(cds.pgsplit)))

ha <- HeatmapAnnotation(df = df,col = list(Cluster = ClusterCol, 
                                           Pseudotime = colorRamp2(seq(from=0,to=max_pseudo,length=11),rev(brewer.pal(11, "RdBu")))),
                        annotation_name_side = "left")


df1 <- data.frame(Cluster = all_cluster[names(all_cluster) %in% id])
df2 <- data.frame(Pseudotime = all_pseudotime[names(all_pseudotime)  %in% id])

ha1 <- HeatmapAnnotation(df = df1, col = list(Cluster = ClusterCol), 
                         annotation_name_side = "left")
ha2 <- HeatmapAnnotation(df = df2, col = list(Pseudotime = colorRamp2(seq(from=0,to=max_pseudo,length=11),rev(brewer.pal(11, "RdBu")))), 
                         annotation_name_side = "left")
#K means with 6 groups
htkm <- Heatmap(
  as.matrix(pt.matrix),
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  km = 6,
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  top_annotation               = ha
)

htkm
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
  cluster_rows                 = FALSE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  top_annotation               = ha1,
  bottom_annotation            = ha2
)

hthc

pdf("20220614_s4hV6Sub_Monocle3PGsplit.heatmap2.SC12.uni-zscore.pdf",width=6,height = 2)
print(htkm)
print(hthc)
dev.off()

##### 20220714
library(ggplot2)
library(monocle3)

cds.pgsplit <- readRDS("../RDS/cds.pgsplit.RDS")
plot_cells(cds.pgsplit)

cds.pgsplit.t <- order_cells(cds.pgsplit, root_pr_nodes = c("Y_95","Y_126"))

plot_cells(cds.pgsplit.t, 
           color_cells_by = "pseudotime",
                               label_cell_groups=FALSE,
                              label_leaves=FALSE,
                              label_branch_points=FALSE,
                              graph_label_size=1.5)+ggtitle("Root Y_95 Y_126")

#### 20220906
library(ggplot2)
library(patchwork)
library(Seurat)

#marker1
marker <- read.table("20220960.2.txt", header = F, sep = "\t")

project.insplit <- readRDS("D:/Projects/wj_sc_project2/RDS/project.insplit.RDS")

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

plot.data <- project.insplit@assays$RNA@scale.data[rownames(project.insplit@assays$RNA@scale.data) %in% rev(marker$V1),]
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
plot.data.final$sample <-project.insplit@meta.data$Sample
dim(plot.data.final)

plot.data.final <- reshape2::melt(plot.data.final,id.vars=c("barcode","cluster","sample"))
head(plot.data.final)
colnames(plot.data.final) <- c("barcode","cluster","sample","symbol","expression")
plot.data.final$symbol <- factor(plot.data.final$symbol, levels = rev(marker$V1))
plot.data.final$cluster <- factor(plot.data.final$cluster, 
                                  levels = rev(c(20,16,17,8,21,2,6,22,18,1,28,12,19)))
plot.data.final$type[plot.data.final$cluster %in% c(20,16,17,8,21,2,6)] <- "dorsal"
plot.data.final$type[plot.data.final$cluster %in% c(18,1,28,12,19)] <- "ventral"
plot.data.final$type[plot.data.final$cluster %in% c(22)] <- "mix"

col <- c("#FA9A99", "#1F78B4", "#B1DF89")

p <-ggplot(plot.data.final)+
  geom_violin(aes(x=expression,y=cluster,fill = type),scale="width",width=0.75)+
  theme_bw()+ ylab("")+xlab("")+
  facet_grid(.~symbol) +
  theme(plot.margin=unit(x=c(0.05,0,0,0),units="inches"),strip.background = element_blank()) +
  theme(strip.text.x = element_text(angle = 270,face = "italic")) +
  theme(axis.line.y = element_blank(),axis.text.y=element_text(color = "black",angle = 270)) +
  theme(axis.line.x = element_blank(),axis.ticks = element_blank()) +
  theme(axis.text.x = element_blank(),) +
  theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid=element_blank()) +
  scale_fill_manual(values = col) + ggtitle("All")
# theme(legend.position = "none")
p

ggsave("20220906.IN.miniviolin1.pdf",p,width=6,height=6)

p1 <- ggplot(plot.data.final[plot.data.final$sample == "SC8",])+
  geom_violin(aes(x=expression,y=cluster,fill = type),scale="width",width=0.75)+
  theme_bw()+ ylab("")+xlab("")+
  facet_grid(.~symbol) +
  theme(plot.margin=unit(x=c(0.05,0,0,0),units="inches"),strip.background = element_blank()) +
  theme(strip.text.x = element_text(angle = 270,face = "italic")) +
  theme(axis.line.y = element_blank(),axis.text.y=element_text(color = "black",angle = 270)) +
  theme(axis.line.x = element_blank(),axis.ticks = element_blank()) +
  theme(axis.text.x = element_blank(),) +
  theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid=element_blank()) +
  scale_fill_manual(values = col) +ggtitle("SC8")
p1

ggsave("20220906.IN.miniviolin1.SC8.pdf",p1,width=6,height=6)

p2 <- ggplot(plot.data.final[plot.data.final$sample == "SC9",])+
  geom_violin(aes(x=expression,y=cluster,fill = type),scale="width",width=0.75)+
  theme_bw()+ ylab("")+xlab("")+
  facet_grid(.~symbol) +
  theme(plot.margin=unit(x=c(0.05,0,0,0),units="inches"),strip.background = element_blank()) +
  theme(strip.text.x = element_text(angle = 270,face = "italic")) +
  theme(axis.line.y = element_blank(),axis.text.y=element_text(color = "black",angle = 270)) +
  theme(axis.line.x = element_blank(),axis.ticks = element_blank()) +
  theme(axis.text.x = element_blank(),) +
  theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid=element_blank()) +
  scale_fill_manual(values = col) +ggtitle("SC9")
p2

ggsave("20220906.IN.miniviolin1.SC9.pdf",p2,width=6,height=6)

p3 <- ggplot(plot.data.final[plot.data.final$sample == "SC10",])+
  geom_violin(aes(x=expression,y=cluster,fill = type),scale="width",width=0.75)+
  theme_bw()+ ylab("")+xlab("")+
  facet_grid(.~symbol) +
  theme(plot.margin=unit(x=c(0.05,0,0,0),units="inches"),strip.background = element_blank()) +
  theme(strip.text.x = element_text(angle = 270,face = "italic")) +
  theme(axis.line.y = element_blank(),axis.text.y=element_text(color = "black",angle = 270)) +
  theme(axis.line.x = element_blank(),axis.ticks = element_blank()) +
  theme(axis.text.x = element_blank(),) +
  theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid=element_blank()) +
  scale_fill_manual(values = col) +ggtitle("SC10")
p3

ggsave("20220906.IN.miniviolin1.SC10.pdf",p3,width=6,height=6)

p4 <- ggplot(plot.data.final[plot.data.final$sample == "SC12",])+
  geom_violin(aes(x=expression,y=cluster,fill = type),scale="width",width=0.75)+
  theme_bw()+ ylab("")+xlab("")+
  facet_grid(.~symbol) +
  theme(plot.margin=unit(x=c(0.05,0,0,0),units="inches"),strip.background = element_blank()) +
  theme(strip.text.x = element_text(angle = 270,face = "italic")) +
  theme(axis.line.y = element_blank(),axis.text.y=element_text(color = "black",angle = 270)) +
  theme(axis.line.x = element_blank(),axis.ticks = element_blank()) +
  theme(axis.text.x = element_blank(),) +
  theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid=element_blank()) +
  scale_fill_manual(values = col) +ggtitle("SC12")
p4

ggsave("20220906.IN.miniviolin1.SC12.pdf",p4,width=6,height=6)

p.all <- p/p1/p2/p3/p4 +plot_layout(ncol = 1, nrow = 5)

ggsave("20220906.IN.miniviolin1.sample.pdf",p.all,width=6,height=30)

#marker2
marker <- read.table("20220960.3.txt", header = F, sep = "\t")

project.insplit <- readRDS("D:/Projects/wj_sc_project2/RDS/project.insplit.RDS")

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

plot.data <- project.insplit@assays$RNA@scale.data[rownames(project.insplit@assays$RNA@scale.data) %in% rev(marker$V1),]
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
plot.data.final$sample <-project.insplit@meta.data$Sample
dim(plot.data.final)

plot.data.final <- reshape2::melt(plot.data.final,id.vars=c("barcode","cluster","sample"))
head(plot.data.final)
colnames(plot.data.final) <- c("barcode","cluster","sample","symbol","expression")
plot.data.final$symbol <- factor(plot.data.final$symbol, levels = rev(marker$V1))
plot.data.final$cluster <- factor(plot.data.final$cluster, 
                                  levels = rev(c(20,16,17,8,21,2,6,22,18,1,28,12,19)))
plot.data.final$type[plot.data.final$cluster %in% c(20,16,17,8,21,2,6)] <- "dorsal"
plot.data.final$type[plot.data.final$cluster %in% c(18,1,28,12,19)] <- "ventral"
plot.data.final$type[plot.data.final$cluster %in% c(22)] <- "mix"

col <- c("#FA9A99", "#1F78B4", "#B1DF89")

p <-ggplot(plot.data.final)+
  geom_violin(aes(x=expression,y=cluster,fill = type),scale="width",width=0.75)+
  theme_bw()+ ylab("")+xlab("")+
  facet_grid(.~symbol) +
  theme(plot.margin=unit(x=c(0.05,0,0,0),units="inches"),strip.background = element_blank()) +
  theme(strip.text.x = element_text(angle = 270,face = "italic")) +
  theme(axis.line.y = element_blank(),axis.text.y=element_text(color = "black",angle = 270)) +
  theme(axis.line.x = element_blank(),axis.ticks = element_blank()) +
  theme(axis.text.x = element_blank(),) +
  theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid=element_blank()) +
  scale_fill_manual(values = col) + ggtitle("All")
# theme(legend.position = "none")
p

ggsave("20220906.IN.miniviolin2.pdf",p,width=6,height=6)

p1 <- ggplot(plot.data.final[plot.data.final$sample == "SC8",])+
  geom_violin(aes(x=expression,y=cluster,fill = type),scale="width",width=0.75)+
  theme_bw()+ ylab("")+xlab("")+
  facet_grid(.~symbol) +
  theme(plot.margin=unit(x=c(0.05,0,0,0),units="inches"),strip.background = element_blank()) +
  theme(strip.text.x = element_text(angle = 270,face = "italic")) +
  theme(axis.line.y = element_blank(),axis.text.y=element_text(color = "black",angle = 270)) +
  theme(axis.line.x = element_blank(),axis.ticks = element_blank()) +
  theme(axis.text.x = element_blank(),) +
  theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid=element_blank()) +
  scale_fill_manual(values = col) +ggtitle("SC8")
p1

ggsave("20220906.IN.miniviolin2.SC8.pdf",p1,width=6,height=6)

p2 <- ggplot(plot.data.final[plot.data.final$sample == "SC9",])+
  geom_violin(aes(x=expression,y=cluster,fill = type),scale="width",width=0.75)+
  theme_bw()+ ylab("")+xlab("")+
  facet_grid(.~symbol) +
  theme(plot.margin=unit(x=c(0.05,0,0,0),units="inches"),strip.background = element_blank()) +
  theme(strip.text.x = element_text(angle = 270,face = "italic")) +
  theme(axis.line.y = element_blank(),axis.text.y=element_text(color = "black",angle = 270)) +
  theme(axis.line.x = element_blank(),axis.ticks = element_blank()) +
  theme(axis.text.x = element_blank(),) +
  theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid=element_blank()) +
  scale_fill_manual(values = col) +ggtitle("SC9")
p2

ggsave("20220906.IN.miniviolin2.SC9.pdf",p2,width=6,height=6)

p3 <- ggplot(plot.data.final[plot.data.final$sample == "SC10",])+
  geom_violin(aes(x=expression,y=cluster,fill = type),scale="width",width=0.75)+
  theme_bw()+ ylab("")+xlab("")+
  facet_grid(.~symbol) +
  theme(plot.margin=unit(x=c(0.05,0,0,0),units="inches"),strip.background = element_blank()) +
  theme(strip.text.x = element_text(angle = 270,face = "italic")) +
  theme(axis.line.y = element_blank(),axis.text.y=element_text(color = "black",angle = 270)) +
  theme(axis.line.x = element_blank(),axis.ticks = element_blank()) +
  theme(axis.text.x = element_blank(),) +
  theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid=element_blank()) +
  scale_fill_manual(values = col) +ggtitle("SC10")
p3

ggsave("20220906.IN.miniviolin2.SC10.pdf",p3,width=6,height=6)

p4 <- ggplot(plot.data.final[plot.data.final$sample == "SC12",])+
  geom_violin(aes(x=expression,y=cluster,fill = type),scale="width",width=0.75)+
  theme_bw()+ ylab("")+xlab("")+
  facet_grid(.~symbol) +
  theme(plot.margin=unit(x=c(0.05,0,0,0),units="inches"),strip.background = element_blank()) +
  theme(strip.text.x = element_text(angle = 270,face = "italic")) +
  theme(axis.line.y = element_blank(),axis.text.y=element_text(color = "black",angle = 270)) +
  theme(axis.line.x = element_blank(),axis.ticks = element_blank()) +
  theme(axis.text.x = element_blank(),) +
  theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid=element_blank()) +
  scale_fill_manual(values = col) +ggtitle("SC12")
p4

ggsave("20220906.IN.miniviolin2.SC12.pdf",p4,width=6,height=6)

p.all <- p/p1/p2/p3/p4 +plot_layout(ncol = 1, nrow = 5)

ggsave("20220906.IN.miniviolin2.sample.pdf",p.all,width=6,height=30)

#### 20220906 pIN split sample violin
library(Seurat)
library(ggplot2)
library(patchwork)

marker <- read.table("20220960.1.txt", header = F, sep = "\t")

project.insplit <- readRDS("D:/Projects/wj_sc_project2/RDS/project.insplit.RDS")

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

plot.data <- project.insplit@assays$RNA@scale.data[rownames(project.insplit@assays$RNA@scale.data) %in% rev(marker$V1),]
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
plot.data.final$sample <-project.insplit@meta.data$Sample
dim(plot.data.final)

plot.data.final <- reshape2::melt(plot.data.final,id.vars=c("barcode","cluster","sample"))
head(plot.data.final)
colnames(plot.data.final) <- c("barcode","cluster","sample","symbol","expression")
plot.data.final$symbol <- factor(plot.data.final$symbol, levels = marker$V1)
plot.data.final$cluster <- factor(plot.data.final$cluster, 
                                  levels = rev(c(20,16,17,8,21,2,6,22,18,1,28,12,19)))
plot.data.final$type[plot.data.final$cluster %in% c(20,16,17,8,21,2,6)] <- "dorsal"
plot.data.final$type[plot.data.final$cluster %in% c(18,1,28,12,19)] <- "ventral"
plot.data.final$type[plot.data.final$cluster %in% c(22)] <- "mix"

col <- c("#FA9A99", "#1F78B4", "#B1DF89")

p <-ggplot(plot.data.final)+
  geom_violin(aes(x=expression,y=cluster,fill = type),scale="width",width=0.75)+
  theme_bw()+ ylab("")+xlab("")+
  facet_grid(.~symbol) +
  theme(plot.margin=unit(x=c(0.05,0,0,0),units="inches"),strip.background = element_blank()) +
  theme(strip.text.x = element_text(angle = 270,face = "italic")) +
  theme(axis.line.y = element_blank(),axis.text.y=element_text(color = "black",angle = 270)) +
  theme(axis.line.x = element_blank(),axis.ticks = element_blank()) +
  theme(axis.text.x = element_blank(),) +
  theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid=element_blank()) +
  scale_fill_manual(values = col) + ggtitle("All")
# theme(legend.position = "none")
p

ggsave("20220906.INsplit.violin.All.pdf",p,width=12,height=6)

p1 <- ggplot(plot.data.final[plot.data.final$sample == "SC8",])+
  geom_violin(aes(x=expression,y=cluster,fill = type),scale="width",width=0.75)+
  theme_bw()+ ylab("")+xlab("")+
  facet_grid(.~symbol) +
  theme(plot.margin=unit(x=c(0.05,0,0,0),units="inches"),strip.background = element_blank()) +
  theme(strip.text.x = element_text(angle = 270,face = "italic")) +
  theme(axis.line.y = element_blank(),axis.text.y=element_text(color = "black",angle = 270)) +
  theme(axis.line.x = element_blank(),axis.ticks = element_blank()) +
  theme(axis.text.x = element_blank(),) +
  theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid=element_blank()) +
  scale_fill_manual(values = col) +ggtitle("SC8")
p1

ggsave("20220906.INsplit.violin.SC8.pdf",p1,width=12,height=6)

p2 <- ggplot(plot.data.final[plot.data.final$sample == "SC9",])+
  geom_violin(aes(x=expression,y=cluster,fill = type),scale="width",width=0.75)+
  theme_bw()+ ylab("")+xlab("")+
  facet_grid(.~symbol) +
  theme(plot.margin=unit(x=c(0.05,0,0,0),units="inches"),strip.background = element_blank()) +
  theme(strip.text.x = element_text(angle = 270,face = "italic")) +
  theme(axis.line.y = element_blank(),axis.text.y=element_text(color = "black",angle = 270)) +
  theme(axis.line.x = element_blank(),axis.ticks = element_blank()) +
  theme(axis.text.x = element_blank(),) +
  theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid=element_blank()) +
  scale_fill_manual(values = col) +ggtitle("SC9")
p2

ggsave("20220906.INsplit.violin.SC9.pdf",p2,width=12,height=6)

p3 <- ggplot(plot.data.final[plot.data.final$sample == "SC10",])+
  geom_violin(aes(x=expression,y=cluster,fill = type),scale="width",width=0.75)+
  theme_bw()+ ylab("")+xlab("")+
  facet_grid(.~symbol) +
  theme(plot.margin=unit(x=c(0.05,0,0,0),units="inches"),strip.background = element_blank()) +
  theme(strip.text.x = element_text(angle = 270,face = "italic")) +
  theme(axis.line.y = element_blank(),axis.text.y=element_text(color = "black",angle = 270)) +
  theme(axis.line.x = element_blank(),axis.ticks = element_blank()) +
  theme(axis.text.x = element_blank(),) +
  theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid=element_blank()) +
  scale_fill_manual(values = col) +ggtitle("SC10")
p3

ggsave("20220906.INsplit.violin.SC10.pdf",p3,width=12,height=6)

p4 <- ggplot(plot.data.final[plot.data.final$sample == "SC12",])+
  geom_violin(aes(x=expression,y=cluster,fill = type),scale="width",width=0.75)+
  theme_bw()+ ylab("")+xlab("")+
  facet_grid(.~symbol) +
  theme(plot.margin=unit(x=c(0.05,0,0,0),units="inches"),strip.background = element_blank()) +
  theme(strip.text.x = element_text(angle = 270,face = "italic")) +
  theme(axis.line.y = element_blank(),axis.text.y=element_text(color = "black",angle = 270)) +
  theme(axis.line.x = element_blank(),axis.ticks = element_blank()) +
  theme(axis.text.x = element_blank(),) +
  theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid=element_blank()) +
  scale_fill_manual(values = col) +ggtitle("SC12")
p4

ggsave("20220906.INsplit.violin.SC12.pdf",p4,width=12,height=6)

p.all <- p/p1/p2/p3/p4 +plot_layout(ncol = 1, nrow = 5)

ggsave("20220906.INsplit.violin.sample.pdf",p.all,width=12,height=30)
