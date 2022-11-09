setwd("D:/Projects/wj_sc_project2/20210803_s4hV6Sub_Monocle3")
setwd("/data2/chenzixi/20210419_WJ_SC/20210803_s4hV6Sub_Monocle3")

##### import libs
library(Seurat)
library(reshape2)
library(gplots)
library(hash)
library(data.table)
library(grid)
library(RColorBrewer)
library(ComplexHeatmap)
library(scales)
library(magrittr)
library(dplyr)
library(harmony)
library(ggplot2)
library(clustree)
library(monocle3)
library(viridis)
library(patchwork)
library(plyr)
#memory.limit(size=100000)

cols=c(brewer.pal(12,"Paired"),brewer.pal(8,"Dark2"),
       brewer.pal(9,"Set1"),brewer.pal(8,"Set2"),
       brewer.pal(12,"Set3"),brewer.pal(8,"Pastel2"),
       brewer.pal(9,"Pastel1"),brewer.pal(8,"Accent"))
col=unique(cols)

#####
project.monsub <- subset(project, idents = c(5,9,12,13,23))

##### import Seurat data to Monocle
data <- GetAssayData(object = project.monsub, assay ="RNA", slot = "counts")
celldata <- as.data.frame(project.monsub@meta.data)
genedata <- as.data.frame(x = row.names(data), row.names = row.names(data))
colnames(genedata) <- "gene_short_name"
cds <- new_cell_data_set(data, cell_metadata = celldata, gene_metadata = genedata)

##### Estimate size factors and dispersions
cds <- estimate_size_factors(cds)
cds <- preprocess_cds(cds)
cds <- align_cds(cds, alignment_group = "Sample")
cds <- reduce_dimension(cds,reduction_method = "UMAP")

# ##### use Seurat umap
# cds.embed <- cds@int_colData$reducedDims$UMAP
# int.embed <- Embeddings(project, reduction = "umap")
# int.embed <- int.embed[rownames(cds.embed),]
# cds@int_colData$reducedDims$UMAP <- int.embed

cds <- cluster_cells(cds, verbose = T)

r1 <- plot_cells(cds, reduction_method = "UMAP", graph_label_size=10, show_trajectory_graph = F,
                 color_cells_by = "seurat_clusters", 
                 label_cell_groups = T, label_groups_by_cluster = T,
                 label_branch_points = F, label_roots = F, label_leaves = F,
                 label_principal_points = F) + ggtitle("Seurat_cluster")+
  scale_color_manual(values = col)
r1$theme$legend.position <- "right"

r2 <- plot_cells(cds, reduction_method = "UMAP", graph_label_size=10, show_trajectory_graph = F,
                 color_cells_by = "cluster", 
                 label_cell_groups = F, label_groups_by_cluster = T,
                 label_branch_points = F, label_roots = F, label_leaves = F,
                 label_principal_points = F) + ggtitle("Monocle_cluster")+
  scale_color_manual(values = col)

pdf("20210830_s4hV6Sub_Monocle3.pdf",width = 18,height = 6)
print(r1 + r2)
dev.off()

### learn_graph
cds <- learn_graph(cds)

pdf("20210830_s4hV6Sub_Monocle3.graph.pdf")

plot_cells(cds, reduction_method = "UMAP", graph_label_size=5, show_trajectory_graph = F,
           color_cells_by = "cluster", 
           label_cell_groups = T, label_groups_by_cluster = T,
           label_branch_points = F, label_roots = F, label_leaves = F,
           label_principal_points = F) + ggtitle("No_Label")+
  scale_color_manual(values = col) + coord_fixed()

plot_cells(cds, reduction_method = "UMAP", graph_label_size=5, show_trajectory_graph = F,
           color_cells_by = "partition", 
           label_cell_groups = T, label_groups_by_cluster = T,
           label_branch_points = F, label_roots = F, label_leaves = F,
           label_principal_points = F) + ggtitle("Partition")+
  scale_color_manual(values = col) + coord_fixed()

plot_cells(cds, reduction_method = "UMAP", graph_label_size=3, 
           color_cells_by = "cluster", 
           label_cell_groups = T, label_groups_by_cluster = F,
           label_branch_points = F, label_roots = T, label_leaves = F,
           label_principal_points = F) + ggtitle("Label_roots")+
  scale_color_manual(values = col) + coord_fixed()

plot_cells(cds, reduction_method = "UMAP", graph_label_size=3, 
           color_cells_by = "cluster", 
           label_cell_groups = T, label_groups_by_cluster = F,
           label_branch_points = F, label_roots = F, label_leaves = T,
           label_principal_points = F) + ggtitle("Label_leaves")+
  scale_color_manual(values = col) + coord_fixed()

plot_cells(cds, reduction_method = "UMAP", graph_label_size=3, 
           color_cells_by = "cluster", 
           label_cell_groups = T, label_groups_by_cluster = F,
           label_branch_points = T, label_roots = F, label_leaves = F,
           label_principal_points = F) + ggtitle("Label_branch_points")+
  scale_color_manual(values = col) + coord_fixed()

plot_cells(cds, reduction_method = "UMAP", graph_label_size=2, 
           color_cells_by = "cluster", 
           label_cell_groups = F, label_groups_by_cluster = T,
           label_branch_points = F, label_roots = F, label_leaves = F,
           label_principal_points = T ) +   ggtitle("label_principal_points")+
  scale_color_manual(values = col) + coord_fixed()

dev.off()

# cds <- order_cells(cds, root_pr_nodes= c("Y_18","Y_160","Y_175"))
# 
# pse1 <- plot_cells(cds,
#                    color_cells_by = "pseudotime",
#                    label_cell_groups=FALSE,
#                    label_leaves=FALSE,
#                    label_branch_points=FALSE,
#                    graph_label_size=1.5)+ggtitle("Root Y_18 Y_160 Y_175")

# ggsave("20210725.s4hV6Sub.MonSub.MN.Monocle3.pseudotime1.pdf",pse1)
### plot marker gene
marker <- read.table("marker1.txt",header = F)[,1]
p <- plot_cells(cds, genes = marker, 
                show_trajectory_graph=FALSE, 
                label_cell_groups=FALSE)+
  scale_color_viridis(option = "C")+ coord_fixed()
ggsave("20210803_s4hV6Sub_Monocle3.marker1.pdf",p,height=32,width = 32)

marker <- read.table("marker2.txt",header = F)[,1]
p<-plot_cells(cds, genes = marker, 
              show_trajectory_graph=FALSE, 
              label_cell_groups=FALSE)+
  scale_color_viridis(option = "C")+ coord_fixed()
ggsave("20210803_s4hV6Sub_Monocle3.marker2.pdf",p,height=24,width = 32)

marker <- read.table("marker3.txt",header = F)[,1]
p<-plot_cells(cds, genes = unique(marker), 
              show_trajectory_graph=FALSE, 
              label_cell_groups=FALSE)+
  scale_color_viridis(option = "C")+ coord_fixed()
ggsave("20210803_s4hV6Sub_Monocle3.marker3.pdf",p,height=48,width = 56,limitsize = FALSE)

marker <- read.table("marker4.txt",header = F)[,1]
p<-plot_cells(cds, genes = unique(marker), 
              show_trajectory_graph=FALSE, 
              label_cell_groups=FALSE)+
  scale_color_viridis(option = "C")+ coord_fixed()
ggsave("20210803_s4hV6Sub_Monocle3.marker4.pdf",p,height=32,width = 32)

##### 
# IN
# 5,9,12,13,23,  1,3,7
# MN
# 5,9,12,13,23, 11,21
# neuron 
# 1.3.7.11.21
#####

#### IN
#####
project.monsub <- subset(project, idents = c(1,3,7,5,9,12,13,23))

##### import Seurat data to Monocle
data <- GetAssayData(object = project.monsub, assay ="RNA", slot = "counts")
celldata <- as.data.frame(project.monsub@meta.data)
genedata <- as.data.frame(x = row.names(data), row.names = row.names(data))
colnames(genedata) <- "gene_short_name"
cds <- new_cell_data_set(data, cell_metadata = celldata, gene_metadata = genedata)

##### Estimate size factors and dispersions
cds <- estimate_size_factors(cds)
cds <- preprocess_cds(cds)
cds <- align_cds(cds, alignment_group = "Sample")
cds <- reduce_dimension(cds,reduction_method = "UMAP")

# ##### use Seurat umap
# cds.embed <- cds@int_colData$reducedDims$UMAP
# int.embed <- Embeddings(project, reduction = "umap")
# int.embed <- int.embed[rownames(cds.embed),]
# cds@int_colData$reducedDims$UMAP <- int.embed

cds <- cluster_cells(cds, verbose = T)

r1 <- plot_cells(cds, reduction_method = "UMAP", graph_label_size=10, show_trajectory_graph = F,
                 color_cells_by = "seurat_clusters", 
                 label_cell_groups = T, label_groups_by_cluster = T,
                 label_branch_points = F, label_roots = F, label_leaves = F,
                 label_principal_points = F) + ggtitle("Seurat_cluster")+
  scale_color_manual(values = col)
r1$theme$legend.position <- "right"

r2 <- plot_cells(cds, reduction_method = "UMAP", graph_label_size=10, show_trajectory_graph = F,
                 color_cells_by = "cluster", 
                 label_cell_groups = F, label_groups_by_cluster = T,
                 label_branch_points = F, label_roots = F, label_leaves = F,
                 label_principal_points = F) + ggtitle("Monocle_cluster")+
  scale_color_manual(values = col)

pdf("20210803_s4hV6Sub_Monocle3.IN.pdf",width = 18,height = 6)
print(r1 + r2)
dev.off()

### learn_graph
cds <- learn_graph(cds)

pdf("20210803_s4hV6Sub_Monocle3.IN.graph.pdf")

plot_cells(cds, reduction_method = "UMAP", graph_label_size=5, show_trajectory_graph = F,
           color_cells_by = "cluster", 
           label_cell_groups = T, label_groups_by_cluster = T,
           label_branch_points = F, label_roots = F, label_leaves = F,
           label_principal_points = F) + ggtitle("No_Label")+
  scale_color_manual(values = col) + coord_fixed()

plot_cells(cds, reduction_method = "UMAP", graph_label_size=5, show_trajectory_graph = F,
           color_cells_by = "partition", 
           label_cell_groups = T, label_groups_by_cluster = T,
           label_branch_points = F, label_roots = F, label_leaves = F,
           label_principal_points = F) + ggtitle("Partition")+
  scale_color_manual(values = col) + coord_fixed()

plot_cells(cds, reduction_method = "UMAP", graph_label_size=3, 
           color_cells_by = "cluster", 
           label_cell_groups = T, label_groups_by_cluster = F,
           label_branch_points = F, label_roots = T, label_leaves = F,
           label_principal_points = F) + ggtitle("Label_roots")+
  scale_color_manual(values = col) + coord_fixed()

plot_cells(cds, reduction_method = "UMAP", graph_label_size=3, 
           color_cells_by = "cluster", 
           label_cell_groups = T, label_groups_by_cluster = F,
           label_branch_points = F, label_roots = F, label_leaves = T,
           label_principal_points = F) + ggtitle("Label_leaves")+
  scale_color_manual(values = col) + coord_fixed()

plot_cells(cds, reduction_method = "UMAP", graph_label_size=3, 
           color_cells_by = "cluster", 
           label_cell_groups = T, label_groups_by_cluster = F,
           label_branch_points = T, label_roots = F, label_leaves = F,
           label_principal_points = F) + ggtitle("Label_branch_points")+
  scale_color_manual(values = col) + coord_fixed()

plot_cells(cds, reduction_method = "UMAP", graph_label_size=2, 
           color_cells_by = "cluster", 
           label_cell_groups = F, label_groups_by_cluster = T,
           label_branch_points = F, label_roots = F, label_leaves = F,
           label_principal_points = T ) +   ggtitle("label_principal_points")+
  scale_color_manual(values = col) + coord_fixed()

dev.off()

# cds <- order_cells(cds, root_pr_nodes= c("Y_18","Y_160","Y_175"))
# 
# pse1 <- plot_cells(cds,
#                    color_cells_by = "pseudotime",
#                    label_cell_groups=FALSE,
#                    label_leaves=FALSE,
#                    label_branch_points=FALSE,
#                    graph_label_size=1.5)+ggtitle("Root Y_18 Y_160 Y_175")

# ggsave("20210725.s4hV6Sub.MonSub.MN.Monocle3.pseudotime1.pdf",pse1)
### plot marker gene
marker <- read.table("marker1.txt",header = F)[,1]
p <- plot_cells(cds, genes = marker, 
                show_trajectory_graph=FALSE, 
                label_cell_groups=FALSE)+
  scale_color_viridis(option = "C")+ coord_fixed()
ggsave("20210803_s4hV6Sub_Monocle3.IN.marker1.pdf",p,height=32,width = 32)

marker <- read.table("marker2.txt",header = F)[,1]
p<-plot_cells(cds, genes = marker, 
              show_trajectory_graph=FALSE, 
              label_cell_groups=FALSE)+
  scale_color_viridis(option = "C")+ coord_fixed()
ggsave("20210803_s4hV6Sub_Monocle3.IN.marker2.pdf",p,height=24,width = 32)

marker <- read.table("marker3.txt",header = F)[,1]
p<-plot_cells(cds, genes = unique(marker), 
              show_trajectory_graph=FALSE, 
              label_cell_groups=FALSE)+
  scale_color_viridis(option = "C")+ coord_fixed()
ggsave("20210803_s4hV6Sub_Monocle3.IN.marker3.pdf",p,height=48,width = 56,limitsize = FALSE)

marker <- read.table("marker4.txt",header = F)[,1]
p<-plot_cells(cds, genes = unique(marker), 
              show_trajectory_graph=FALSE, 
              label_cell_groups=FALSE)+
  scale_color_viridis(option = "C")+ coord_fixed()
ggsave("20210803_s4hV6Sub_Monocle3.IN.marker4.pdf",p,height=32,width = 32)

##### MN
#####
project.monsub <- subset(project, idents = c(11,21,5,9,12,13,23))

##### import Seurat data to Monocle
data <- GetAssayData(object = project.monsub, assay ="RNA", slot = "counts")
celldata <- as.data.frame(project.monsub@meta.data)
genedata <- as.data.frame(x = row.names(data), row.names = row.names(data))
colnames(genedata) <- "gene_short_name"
cds <- new_cell_data_set(data, cell_metadata = celldata, gene_metadata = genedata)

##### Estimate size factors and dispersions
cds <- estimate_size_factors(cds)
cds <- preprocess_cds(cds)
cds <- align_cds(cds, alignment_group = "Sample")
cds <- reduce_dimension(cds,reduction_method = "UMAP")

# ##### use Seurat umap
# cds.embed <- cds@int_colData$reducedDims$UMAP
# int.embed <- Embeddings(project, reduction = "umap")
# int.embed <- int.embed[rownames(cds.embed),]
# cds@int_colData$reducedDims$UMAP <- int.embed

cds <- cluster_cells(cds, verbose = T)

r1 <- plot_cells(cds, reduction_method = "UMAP", graph_label_size=10, show_trajectory_graph = F,
                 color_cells_by = "seurat_clusters", 
                 label_cell_groups = T, label_groups_by_cluster = T,
                 label_branch_points = F, label_roots = F, label_leaves = F,
                 label_principal_points = F) + ggtitle("Seurat_cluster")+
  scale_color_manual(values = col)
r1$theme$legend.position <- "right"

r2 <- plot_cells(cds, reduction_method = "UMAP", graph_label_size=10, show_trajectory_graph = F,
                 color_cells_by = "cluster", 
                 label_cell_groups = F, label_groups_by_cluster = T,
                 label_branch_points = F, label_roots = F, label_leaves = F,
                 label_principal_points = F) + ggtitle("Monocle_cluster")+
  scale_color_manual(values = col)

pdf("20210803_s4hV6Sub_Monocle3.MN.pdf",width = 18,height = 6)
print(r1 + r2)
dev.off()

### learn_graph
cds <- learn_graph(cds)

pdf("20210803_s4hV6Sub_Monocle3.MN.graph.pdf")

plot_cells(cds, reduction_method = "UMAP", graph_label_size=5, show_trajectory_graph = F,
           color_cells_by = "cluster", 
           label_cell_groups = T, label_groups_by_cluster = T,
           label_branch_points = F, label_roots = F, label_leaves = F,
           label_principal_points = F) + ggtitle("No_Label")+
  scale_color_manual(values = col) + coord_fixed()

plot_cells(cds, reduction_method = "UMAP", graph_label_size=5, show_trajectory_graph = F,
           color_cells_by = "partition", 
           label_cell_groups = T, label_groups_by_cluster = T,
           label_branch_points = F, label_roots = F, label_leaves = F,
           label_principal_points = F) + ggtitle("Partition")+
  scale_color_manual(values = col) + coord_fixed()

plot_cells(cds, reduction_method = "UMAP", graph_label_size=3, 
           color_cells_by = "cluster", 
           label_cell_groups = T, label_groups_by_cluster = F,
           label_branch_points = F, label_roots = T, label_leaves = F,
           label_principal_points = F) + ggtitle("Label_roots")+
  scale_color_manual(values = col) + coord_fixed()

plot_cells(cds, reduction_method = "UMAP", graph_label_size=3, 
           color_cells_by = "cluster", 
           label_cell_groups = T, label_groups_by_cluster = F,
           label_branch_points = F, label_roots = F, label_leaves = T,
           label_principal_points = F) + ggtitle("Label_leaves")+
  scale_color_manual(values = col) + coord_fixed()

plot_cells(cds, reduction_method = "UMAP", graph_label_size=3, 
           color_cells_by = "cluster", 
           label_cell_groups = T, label_groups_by_cluster = F,
           label_branch_points = T, label_roots = F, label_leaves = F,
           label_principal_points = F) + ggtitle("Label_branch_points")+
  scale_color_manual(values = col) + coord_fixed()

plot_cells(cds, reduction_method = "UMAP", graph_label_size=2, 
           color_cells_by = "cluster", 
           label_cell_groups = F, label_groups_by_cluster = T,
           label_branch_points = F, label_roots = F, label_leaves = F,
           label_principal_points = T ) +   ggtitle("label_principal_points")+
  scale_color_manual(values = col) + coord_fixed()

dev.off()

# cds <- order_cells(cds, root_pr_nodes= c("Y_18","Y_160","Y_175"))
# 
# pse1 <- plot_cells(cds,
#                    color_cells_by = "pseudotime",
#                    label_cell_groups=FALSE,
#                    label_leaves=FALSE,
#                    label_branch_points=FALSE,
#                    graph_label_size=1.5)+ggtitle("Root Y_18 Y_160 Y_175")

# ggsave("20210725.s4hV6Sub.MonSub.MN.Monocle3.pseudotime1.pdf",pse1)
### plot marker gene
marker <- read.table("marker1.txt",header = F)[,1]
p <- plot_cells(cds, genes = marker, 
                show_trajectory_graph=FALSE, 
                label_cell_groups=FALSE)+
  scale_color_viridis(option = "C")+ coord_fixed()
ggsave("20210803_s4hV6Sub_Monocle3.MN.marker1.pdf",p,height=32,width = 32)

marker <- read.table("marker2.txt",header = F)[,1]
p<-plot_cells(cds, genes = marker, 
              show_trajectory_graph=FALSE, 
              label_cell_groups=FALSE)+
  scale_color_viridis(option = "C")+ coord_fixed()
ggsave("20210803_s4hV6Sub_Monocle3.MN.marker2.pdf",p,height=24,width = 32)

marker <- read.table("marker3.txt",header = F)[,1]
p<-plot_cells(cds, genes = unique(marker), 
              show_trajectory_graph=FALSE, 
              label_cell_groups=FALSE)+
  scale_color_viridis(option = "C")+ coord_fixed()
ggsave("20210803_s4hV6Sub_Monocle3.MN.marker3.pdf",p,height=48,width = 56,limitsize = FALSE)

marker <- read.table("marker4.txt",header = F)[,1]
p<-plot_cells(cds, genes = unique(marker), 
              show_trajectory_graph=FALSE, 
              label_cell_groups=FALSE)+
  scale_color_viridis(option = "C")+ coord_fixed()
ggsave("20210803_s4hV6Sub_Monocle3.MN.marker4.pdf",p,height=32,width = 32)

#### Neuron 1.3.7.11.21
#####
project.monsub <- subset(project, idents = c(1,3,7,11,21))

##### import Seurat data to Monocle
data <- GetAssayData(object = project.monsub, assay ="RNA", slot = "counts")
celldata <- as.data.frame(project.monsub@meta.data)
genedata <- as.data.frame(x = row.names(data), row.names = row.names(data))
colnames(genedata) <- "gene_short_name"
cds <- new_cell_data_set(data, cell_metadata = celldata, gene_metadata = genedata)

##### Estimate size factors and dispersions
cds <- estimate_size_factors(cds)
cds <- preprocess_cds(cds)
cds <- align_cds(cds, alignment_group = "Sample")
cds <- reduce_dimension(cds,reduction_method = "UMAP")

# ##### use Seurat umap
# cds.embed <- cds@int_colData$reducedDims$UMAP
# int.embed <- Embeddings(project, reduction = "umap")
# int.embed <- int.embed[rownames(cds.embed),]
# cds@int_colData$reducedDims$UMAP <- int.embed

cds <- cluster_cells(cds, verbose = T)

r1 <- plot_cells(cds, reduction_method = "UMAP", graph_label_size=10, show_trajectory_graph = F,
                 color_cells_by = "seurat_clusters", 
                 label_cell_groups = T, label_groups_by_cluster = T,
                 label_branch_points = F, label_roots = F, label_leaves = F,
                 label_principal_points = F) + ggtitle("Seurat_cluster")+
  scale_color_manual(values = col)
r1$theme$legend.position <- "right"

r2 <- plot_cells(cds, reduction_method = "UMAP", graph_label_size=10, show_trajectory_graph = F,
                 color_cells_by = "cluster", 
                 label_cell_groups = F, label_groups_by_cluster = T,
                 label_branch_points = F, label_roots = F, label_leaves = F,
                 label_principal_points = F) + ggtitle("Monocle_cluster")+
  scale_color_manual(values = col)

pdf("20210803_s4hV6Sub_Monocle3.Neuron.pdf",width = 18,height = 6)
print(r1 + r2)
dev.off()

### learn_graph
cds <- learn_graph(cds)

pdf("20210803_s4hV6Sub_Monocle3.Neuron.graph.pdf")

plot_cells(cds, reduction_method = "UMAP", graph_label_size=5, show_trajectory_graph = F,
           color_cells_by = "cluster", 
           label_cell_groups = T, label_groups_by_cluster = T,
           label_branch_points = F, label_roots = F, label_leaves = F,
           label_principal_points = F) + ggtitle("No_Label")+
  scale_color_manual(values = col) + coord_fixed()

plot_cells(cds, reduction_method = "UMAP", graph_label_size=5, show_trajectory_graph = F,
           color_cells_by = "partition", 
           label_cell_groups = T, label_groups_by_cluster = T,
           label_branch_points = F, label_roots = F, label_leaves = F,
           label_principal_points = F) + ggtitle("Partition")+
  scale_color_manual(values = col) + coord_fixed()

plot_cells(cds, reduction_method = "UMAP", graph_label_size=3, 
           color_cells_by = "cluster", 
           label_cell_groups = T, label_groups_by_cluster = F,
           label_branch_points = F, label_roots = T, label_leaves = F,
           label_principal_points = F) + ggtitle("Label_roots")+
  scale_color_manual(values = col) + coord_fixed()

plot_cells(cds, reduction_method = "UMAP", graph_label_size=3, 
           color_cells_by = "cluster", 
           label_cell_groups = T, label_groups_by_cluster = F,
           label_branch_points = F, label_roots = F, label_leaves = T,
           label_principal_points = F) + ggtitle("Label_leaves")+
  scale_color_manual(values = col) + coord_fixed()

plot_cells(cds, reduction_method = "UMAP", graph_label_size=3, 
           color_cells_by = "cluster", 
           label_cell_groups = T, label_groups_by_cluster = F,
           label_branch_points = T, label_roots = F, label_leaves = F,
           label_principal_points = F) + ggtitle("Label_branch_points")+
  scale_color_manual(values = col) + coord_fixed()

plot_cells(cds, reduction_method = "UMAP", graph_label_size=2, 
           color_cells_by = "cluster", 
           label_cell_groups = F, label_groups_by_cluster = T,
           label_branch_points = F, label_roots = F, label_leaves = F,
           label_principal_points = T ) +   ggtitle("label_principal_points")+
  scale_color_manual(values = col) + coord_fixed()

dev.off()

# cds <- order_cells(cds, root_pr_nodes= c("Y_18","Y_160","Y_175"))
# 
# pse1 <- plot_cells(cds,
#                    color_cells_by = "pseudotime",
#                    label_cell_groups=FALSE,
#                    label_leaves=FALSE,
#                    label_branch_points=FALSE,
#                    graph_label_size=1.5)+ggtitle("Root Y_18 Y_160 Y_175")

# ggsave("20210725.s4hV6Sub.MonSub.MN.Monocle3.pseudotime1.pdf",pse1)
### plot marker gene
marker <- read.table("marker1.txt",header = F)[,1]
p <- plot_cells(cds, genes = marker, 
                show_trajectory_graph=FALSE, 
                label_cell_groups=FALSE)+
  scale_color_viridis(option = "C")+ coord_fixed()
ggsave("20210803_s4hV6Sub_Monocle3.Neuron.marker1.pdf",p,height=32,width = 32)

marker <- read.table("marker2.txt",header = F)[,1]
p<-plot_cells(cds, genes = marker, 
              show_trajectory_graph=FALSE, 
              label_cell_groups=FALSE)+
  scale_color_viridis(option = "C")+ coord_fixed()
ggsave("20210803_s4hV6Sub_Monocle3.Neuron.marker2.pdf",p,height=24,width = 32)

marker <- read.table("marker3.txt",header = F)[,1]
p<-plot_cells(cds, genes = unique(marker), 
              show_trajectory_graph=FALSE, 
              label_cell_groups=FALSE)+
  scale_color_viridis(option = "C")+ coord_fixed()
ggsave("20210803_s4hV6Sub_Monocle3.Neuron.marker3.pdf",p,height=48,width = 56,limitsize = FALSE)

marker <- read.table("marker4.txt",header = F)[,1]
p<-plot_cells(cds, genes = unique(marker), 
              show_trajectory_graph=FALSE, 
              label_cell_groups=FALSE)+
  scale_color_viridis(option = "C")+ coord_fixed()
ggsave("20210803_s4hV6Sub_Monocle3.Neuron.marker4.pdf",p,height=32,width = 32)


#### 20210804 1,2,3,4,5,7,9,11,12,13,14,15,16,18,19,21,23
#####
project.monsub <- subset(project, idents = c(1,2,3,4,5,7,9,11,12,13,14,15,16,18,19,21,23))

##### import Seurat data to Monocle
data <- GetAssayData(object = project.monsub, assay ="RNA", slot = "counts")
celldata <- as.data.frame(project.monsub@meta.data)
genedata <- as.data.frame(x = row.names(data), row.names = row.names(data))
colnames(genedata) <- "gene_short_name"
cds <- new_cell_data_set(data, cell_metadata = celldata, gene_metadata = genedata)

##### Estimate size factors and dispersions
cds <- estimate_size_factors(cds)
cds <- preprocess_cds(cds)
cds <- align_cds(cds, alignment_group = "Sample")
cds <- reduce_dimension(cds,reduction_method = "UMAP")

# ##### use Seurat umap
# cds.embed <- cds@int_colData$reducedDims$UMAP
# int.embed <- Embeddings(project, reduction = "umap")
# int.embed <- int.embed[rownames(cds.embed),]
# cds@int_colData$reducedDims$UMAP <- int.embed

cds <- cluster_cells(cds, verbose = T)

r1 <- plot_cells(cds, reduction_method = "UMAP", graph_label_size=10, show_trajectory_graph = F,
                 color_cells_by = "seurat_clusters", 
                 label_cell_groups = T, label_groups_by_cluster = T,
                 label_branch_points = F, label_roots = F, label_leaves = F,
                 label_principal_points = F) + ggtitle("Seurat_cluster")+
  scale_color_manual(values = col)
r1$theme$legend.position <- "right"

r2 <- plot_cells(cds, reduction_method = "UMAP", graph_label_size=10, show_trajectory_graph = F,
                 color_cells_by = "cluster", 
                 label_cell_groups = F, label_groups_by_cluster = T,
                 label_branch_points = F, label_roots = F, label_leaves = F,
                 label_principal_points = F) + ggtitle("Monocle_cluster")+
  scale_color_manual(values = col)

pdf("20210804_s4hV6Sub_Monocle3.pdf",width = 18,height = 6)
print(r1 + r2)
dev.off()

### learn_graph
cds <- learn_graph(cds)

pdf("20210804_s4hV6Sub_Monocle3.graph.pdf")

plot_cells(cds, reduction_method = "UMAP", graph_label_size=5, show_trajectory_graph = F,
           color_cells_by = "cluster", 
           label_cell_groups = T, label_groups_by_cluster = T,
           label_branch_points = F, label_roots = F, label_leaves = F,
           label_principal_points = F) + ggtitle("No_Label")+
  scale_color_manual(values = col) + coord_fixed()

plot_cells(cds, reduction_method = "UMAP", graph_label_size=5, show_trajectory_graph = F,
           color_cells_by = "partition", 
           label_cell_groups = T, label_groups_by_cluster = T,
           label_branch_points = F, label_roots = F, label_leaves = F,
           label_principal_points = F) + ggtitle("Partition")+
  scale_color_manual(values = col) + coord_fixed()

plot_cells(cds, reduction_method = "UMAP", graph_label_size=3, 
           color_cells_by = "cluster", 
           label_cell_groups = T, label_groups_by_cluster = F,
           label_branch_points = F, label_roots = T, label_leaves = F,
           label_principal_points = F) + ggtitle("Label_roots")+
  scale_color_manual(values = col) + coord_fixed()

plot_cells(cds, reduction_method = "UMAP", graph_label_size=3, 
           color_cells_by = "cluster", 
           label_cell_groups = T, label_groups_by_cluster = F,
           label_branch_points = F, label_roots = F, label_leaves = T,
           label_principal_points = F) + ggtitle("Label_leaves")+
  scale_color_manual(values = col) + coord_fixed()

plot_cells(cds, reduction_method = "UMAP", graph_label_size=3, 
           color_cells_by = "cluster", 
           label_cell_groups = T, label_groups_by_cluster = F,
           label_branch_points = T, label_roots = F, label_leaves = F,
           label_principal_points = F) + ggtitle("Label_branch_points")+
  scale_color_manual(values = col) + coord_fixed()

plot_cells(cds, reduction_method = "UMAP", graph_label_size=2, 
           color_cells_by = "cluster", 
           label_cell_groups = F, label_groups_by_cluster = T,
           label_branch_points = F, label_roots = F, label_leaves = F,
           label_principal_points = T ) +   ggtitle("label_principal_points")+
  scale_color_manual(values = col) + coord_fixed()

dev.off()

# cds <- order_cells(cds, root_pr_nodes= c("Y_18","Y_160","Y_175"))
# 
# pse1 <- plot_cells(cds,
#                    color_cells_by = "pseudotime",
#                    label_cell_groups=FALSE,
#                    label_leaves=FALSE,
#                    label_branch_points=FALSE,
#                    graph_label_size=1.5)+ggtitle("Root Y_18 Y_160 Y_175")

# ggsave("20210725.s4hV6Sub.MonSub.MN.Monocle3.pseudotime1.pdf",pse1)
### plot marker gene
marker <- read.table("marker1.txt",header = F)[,1]
p <- plot_cells(cds, genes = marker, 
                show_trajectory_graph=FALSE, 
                label_cell_groups=FALSE)+
  scale_color_viridis(option = "C")+ coord_fixed()
ggsave("20210804_s4hV6Sub_Monocle3.marker1.pdf",p,height=32,width = 32)

marker <- read.table("marker2.txt",header = F)[,1]
p<-plot_cells(cds, genes = marker, 
              show_trajectory_graph=FALSE, 
              label_cell_groups=FALSE)+
  scale_color_viridis(option = "C")+ coord_fixed()
ggsave("20210804_s4hV6Sub_Monocle3.marker2.pdf",p,height=24,width = 32)

marker <- read.table("marker3.txt",header = F)[,1]
p<-plot_cells(cds, genes = unique(marker), 
              show_trajectory_graph=FALSE, 
              label_cell_groups=FALSE)+
  scale_color_viridis(option = "C")+ coord_fixed()
ggsave("20210804_s4hV6Sub_Monocle3.marker3.pdf",p,height=48,width = 56,limitsize = FALSE)

marker <- read.table("marker4.txt",header = F)[,1]
p<-plot_cells(cds, genes = unique(marker), 
              show_trajectory_graph=FALSE, 
              label_cell_groups=FALSE)+
  scale_color_viridis(option = "C")+ coord_fixed()
ggsave("20210804_s4hV6Sub_Monocle3.marker4.pdf",p,height=32,width = 32)

p1 <- plot_cells(cds, reduction_method = "UMAP", graph_label_size=10, show_trajectory_graph = F,
                 color_cells_by = "seurat_clusters", 
                 label_cell_groups = T, label_groups_by_cluster = T,
                 label_branch_points = F, label_roots = F, label_leaves = F,
                 label_principal_points = F) + ggtitle("Seurat_cluster")+
  scale_color_manual(values = col)+ facet_wrap(~Sample)
p1$theme$legend.position <- "right"

p2<-plot_cells(cds, reduction_method = "UMAP", graph_label_size=10, show_trajectory_graph = F,
                 color_cells_by = "cluster", 
                 label_cell_groups = F, label_groups_by_cluster = T,
                 label_branch_points = F, label_roots = F, label_leaves = F,
                 label_principal_points = F) + ggtitle("Monocle_cluster")+
  scale_color_manual(values = col)+ facet_wrap(~Sample)

pdf("20210803_s4hV6Sub_Monocle3.splitSample.pdf",width = 12,height = 12)
print(p1)
print(p2)
dev.off()

#### 20210804 1,2,3,4,5,7,9,11,12,13,14,15,17,18,21,23
#####
project.monsub <- subset(project, idents = c(1,2,3,4,5,7,9,11,12,13,14,15,17,18,21,23))

##### import Seurat data to Monocle
data <- GetAssayData(object = project.monsub, assay ="RNA", slot = "counts")
celldata <- as.data.frame(project.monsub@meta.data)
genedata <- as.data.frame(x = row.names(data), row.names = row.names(data))
colnames(genedata) <- "gene_short_name"
cds <- new_cell_data_set(data, cell_metadata = celldata, gene_metadata = genedata)

##### Estimate size factors and dispersions
cds <- estimate_size_factors(cds)
cds <- preprocess_cds(cds)
cds <- align_cds(cds, alignment_group = "Sample")
cds <- reduce_dimension(cds,reduction_method = "UMAP")

# ##### use Seurat umap
# cds.embed <- cds@int_colData$reducedDims$UMAP
# int.embed <- Embeddings(project, reduction = "umap")
# int.embed <- int.embed[rownames(cds.embed),]
# cds@int_colData$reducedDims$UMAP <- int.embed

#cds_temp <- cluster_cells(cds, verbose = T, resolution = 3e-04, random_seed = 1024)
cds <- cluster_cells(cds, verbose = T, resolution = 3e-04, random_seed = 1024)

r1 <- plot_cells(cds, reduction_method = "UMAP", graph_label_size=10, show_trajectory_graph = F,
                 color_cells_by = "seurat_clusters", 
                 label_cell_groups = T, label_groups_by_cluster = T,
                 label_branch_points = F, label_roots = F, label_leaves = F,
                 label_principal_points = F) + ggtitle("Seurat_cluster")+
  scale_color_manual(values = col)
r1$theme$legend.position <- "right"

r2 <- plot_cells(cds, reduction_method = "UMAP", graph_label_size=10, show_trajectory_graph = F,
                 color_cells_by = "cluster", 
                 label_cell_groups = F, label_groups_by_cluster = T,
                 label_branch_points = F, label_roots = F, label_leaves = F,
                 label_principal_points = F) + ggtitle("Monocle_cluster")+
  scale_color_manual(values = col)

pdf("20210804_s4hV6Sub2_Monocle3.pdf",width = 18,height = 6)
print(r1 + r2)
dev.off()

### learn_graph
cds <- learn_graph(cds)

pdf("20210804_s4hV6Sub2_Monocle3.graph.pdf")

plot_cells(cds, reduction_method = "UMAP", graph_label_size=5, show_trajectory_graph = F,
           color_cells_by = "cluster", 
           label_cell_groups = T, label_groups_by_cluster = T,
           label_branch_points = F, label_roots = F, label_leaves = F,
           label_principal_points = F) + ggtitle("No_Label")+
  scale_color_manual(values = col) + coord_fixed()

plot_cells(cds, reduction_method = "UMAP", graph_label_size=5, show_trajectory_graph = F,
           color_cells_by = "partition", 
           label_cell_groups = T, label_groups_by_cluster = T,
           label_branch_points = F, label_roots = F, label_leaves = F,
           label_principal_points = F) + ggtitle("Partition")+
  scale_color_manual(values = col) + coord_fixed()

plot_cells(cds, reduction_method = "UMAP", graph_label_size=3, 
           color_cells_by = "cluster", 
           label_cell_groups = T, label_groups_by_cluster = F,
           label_branch_points = F, label_roots = T, label_leaves = F,
           label_principal_points = F) + ggtitle("Label_roots")+
  scale_color_manual(values = col) + coord_fixed()

plot_cells(cds, reduction_method = "UMAP", graph_label_size=3, 
           color_cells_by = "cluster", 
           label_cell_groups = T, label_groups_by_cluster = F,
           label_branch_points = F, label_roots = F, label_leaves = T,
           label_principal_points = F) + ggtitle("Label_leaves")+
  scale_color_manual(values = col) + coord_fixed()

plot_cells(cds, reduction_method = "UMAP", graph_label_size=3, 
           color_cells_by = "cluster", 
           label_cell_groups = T, label_groups_by_cluster = F,
           label_branch_points = T, label_roots = F, label_leaves = F,
           label_principal_points = F) + ggtitle("Label_branch_points")+
  scale_color_manual(values = col) + coord_fixed()

plot_cells(cds, reduction_method = "UMAP", graph_label_size=2, 
           color_cells_by = "cluster", 
           label_cell_groups = F, label_groups_by_cluster = T,
           label_branch_points = F, label_roots = F, label_leaves = F,
           label_principal_points = T ) +   ggtitle("label_principal_points")+
  scale_color_manual(values = col) + coord_fixed()

dev.off()

# cds <- order_cells(cds, root_pr_nodes= c("Y_18","Y_160","Y_175"))
# 
# pse1 <- plot_cells(cds,
#                    color_cells_by = "pseudotime",
#                    label_cell_groups=FALSE,
#                    label_leaves=FALSE,
#                    label_branch_points=FALSE,
#                    graph_label_size=1.5)+ggtitle("Root Y_18 Y_160 Y_175")

# ggsave("20210725.s4hV6Sub.MonSub.MN.Monocle3.pseudotime1.pdf",pse1)
### plot marker gene
marker <- read.table("marker1.txt",header = F)[,1]
p <- plot_cells(cds, genes = marker, 
                show_trajectory_graph=FALSE, 
                label_cell_groups=FALSE)+
  scale_color_viridis(option = "C")+ coord_fixed()
ggsave("20210804_s4hV6Sub2_Monocle3.marker1.pdf",p,height=32,width = 32)

marker <- read.table("marker2.txt",header = F)[,1]
p<-plot_cells(cds, genes = marker, 
              show_trajectory_graph=FALSE, 
              label_cell_groups=FALSE)+
  scale_color_viridis(option = "C")+ coord_fixed()
ggsave("20210804_s4hV6Sub2_Monocle3.marker2.pdf",p,height=24,width = 32)

marker <- read.table("marker3.txt",header = F)[,1]
p<-plot_cells(cds, genes = unique(marker), 
              show_trajectory_graph=FALSE, 
              label_cell_groups=FALSE)+
  scale_color_viridis(option = "C")+ coord_fixed()
ggsave("20210804_s4hV6Sub2_Monocle3.marker3.pdf",p,height=48,width = 56,limitsize = FALSE)

marker <- read.table("marker4.txt",header = F)[,1]
p<-plot_cells(cds, genes = unique(marker), 
              show_trajectory_graph=FALSE, 
              label_cell_groups=FALSE)+
  scale_color_viridis(option = "C")+ coord_fixed()
ggsave("20210804_s4hV6Sub2_Monocle3.marker4.pdf",p,height=32,width = 32)

marker <- read.table("marker5.txt",header = F)[,1]
p<-plot_cells(cds, genes = unique(marker), 
              show_trajectory_graph=FALSE, 
              label_cell_groups=FALSE)+
  scale_color_viridis(option = "C")+ coord_fixed()
ggsave("20210804_s4hV6Sub2_Monocle3.marker5.pdf",p,height=32,width = 32)


p1 <- plot_cells(cds, reduction_method = "UMAP", graph_label_size=10, show_trajectory_graph = F,
                 color_cells_by = "seurat_clusters", 
                 label_cell_groups = T, label_groups_by_cluster = T,
                 label_branch_points = F, label_roots = F, label_leaves = F,
                 label_principal_points = F) + ggtitle("Seurat_cluster")+
  scale_color_manual(values = col)+ facet_wrap(~Sample)
p1$theme$legend.position <- "right"

p2<-plot_cells(cds, reduction_method = "UMAP", graph_label_size=10, show_trajectory_graph = F,
               color_cells_by = "cluster", 
               label_cell_groups = F, label_groups_by_cluster = T,
               label_branch_points = F, label_roots = F, label_leaves = F,
               label_principal_points = F) + ggtitle("Monocle_cluster")+
  scale_color_manual(values = col)+ facet_wrap(~Sample)

pdf("20210803_s4hV6Sub2_Monocle3.splitSample.pdf",width = 12,height = 12)
print(p1)
print(p2)
dev.off()

save.image("20210804_s4hV6Sub2_Monocle3.RData")

##### highlight 
hl_col <- rep("white",length(col))
hl_col[c(1,2)] <- col[c(1,2)]

plot_cells(cds, reduction_method = "UMAP", graph_label_size=10, show_trajectory_graph = F,
           color_cells_by = "cluster", 
           label_cell_groups = F, label_groups_by_cluster = T,
           label_branch_points = F, label_roots = F, label_leaves = F,
           label_principal_points = F) + ggtitle("Monocle_cluster")+
  scale_color_manual(values = hl_col)

##### order
cds <- order_cells(cds, root_pr_nodes = c("Y_269","Y_1154"))

pse1 <- plot_cells(cds,
                   color_cells_by = "pseudotime", trajectory_graph_segment_size = 0.25,
                   label_cell_groups=FALSE,
                   label_leaves=FALSE, 
                   label_branch_points=FALSE,
                   graph_label_size=1.5)+ggtitle("Root Y_269 Y_1154")
pse1
ggsave("20210803_s4hV6Sub2_Monocle3.pseudotime1.pdf",pse1)

### subset for plot
data1 <- cds@clusters$UMAP$clusters
data2 <- monocle3::pseudotime(cds)
data3 <- cds@int_colData@listData$reducedDims$UMAP
data <- as.data.frame(cbind(data1,data2,data3))
colnames(data) <- c("cluster","pseudotime","umap1","umap2")
data$cluster <- factor(as.vector.factor(data$cluster))

sub1 <- subset(data, cluster %in% c(13,14,29,4,5,24,16,22,37,7,9,30,31,20,17,29,1,23,32))
dim(sub1)

p1 <- ggplot(sub1) + geom_point(aes(umap1,umap2,color=pseudotime),size=0.1)+ theme_classic()+
  scale_color_viridis(option = "C") + coord_fixed()

p2 <- ggplot(sub1) + geom_point(aes(umap1,umap2,color=cluster),size=0.1)+ theme_classic()+
   coord_fixed()

p1+p2

### monocle violin
cds_subset <- cds[marker[c(1:10)],]
plot_genes_violin(cds_subset,group_cells_by = "seurat_clusters")

### monocle diff
#24,16,22,37； 13,19； 20,30 ；23,32； 10,36； 2,35； 8,39；
colData(cds)$monocle_cluster <- cds@clusters$UMAP$clusters
colData(cds)$monocle_cluster[colData(cds)$monocle_cluster %in% c(24,16,22,37)] <- 16
colData(cds)$monocle_cluster[colData(cds)$monocle_cluster %in% c(13,19)] <- 13
colData(cds)$monocle_cluster[colData(cds)$monocle_cluster %in% c(20,30)] <- 20
colData(cds)$monocle_cluster[colData(cds)$monocle_cluster %in% c(23,32)] <- 23
colData(cds)$monocle_cluster[colData(cds)$monocle_cluster %in% c(10,36)] <- 10
colData(cds)$monocle_cluster[colData(cds)$monocle_cluster %in% c(2,35)] <- 2
colData(cds)$monocle_cluster[colData(cds)$monocle_cluster %in% c(8,39)] <- 8
colData(cds)$monocle_cluster <- factor(colData(cds)$monocle_cluster)

gene_fits <- fit_models(cds, model_formula_str = "~monocle_cluster")
fit_coefs <- coefficient_table(gene_fits)
monocle_cluster_terms <- fit_coefs %>% filter(term == "monocle_cluster")
monocle_cluster_terms %>% filter (q_value < 0.05) %>% elect(gene_short_name, term, q_value, estimate)

save.image("20210804_s4hV6Sub2_Monocle3.RData")


##### 20210807
moncluster <- cds@clusters$UMAP$clusters

# subset monocle cluster 19,13,14,17,1,23,32
sub_id1 <- names(cds@clusters$UMAP$clusters[cds@clusters$UMAP$clusters %in% c(19,13,14,17,1,23,32)])
length(sub_id1)

sub_cluster1 <- moncluster[moncluster %in% c(19,13,14,17,1,23,32)]
length(sub_cluster1)

project.monmn <- subset(project, cells = sub_id1)
project.monmn

##### import Seurat data to Monocle
data <- GetAssayData(object = project.monmn, assay ="RNA", slot = "counts")
celldata <- as.data.frame(project.monmn@meta.data)
genedata <- as.data.frame(x = row.names(data), row.names = row.names(data))
colnames(genedata) <- "gene_short_name"
cds <- new_cell_data_set(data, cell_metadata = celldata, gene_metadata = genedata)

##### Estimate size factors and dispersions
cds <- estimate_size_factors(cds)
cds <- preprocess_cds(cds)
cds <- align_cds(cds, alignment_group = "Sample")
cds <- reduce_dimension(cds,reduction_method = "UMAP")

head(sub_cluster1)
tail(sub_cluster1)
colData(cds)$moncuster <- sub_cluster1

# ##### use Seurat umap
# cds.embed <- cds@int_colData$reducedDims$UMAP
# int.embed <- Embeddings(project, reduction = "umap")
# int.embed <- int.embed[rownames(cds.embed),]
# cds@int_colData$reducedDims$UMAP <- int.embed

cds <- cluster_cells(cds, verbose = T)

r1 <- plot_cells(cds, reduction_method = "UMAP", graph_label_size=10, show_trajectory_graph = F,
                 color_cells_by = "seurat_clusters", 
                 label_cell_groups = T, label_groups_by_cluster = T,
                 label_branch_points = F, label_roots = F, label_leaves = F,
                 label_principal_points = F) + ggtitle("Seurat_cluster")+
  scale_color_manual(values = col)
r1$theme$legend.position <- "right"

r2 <- plot_cells(cds, reduction_method = "UMAP", graph_label_size=10, show_trajectory_graph = F,
                 color_cells_by = "cluster", 
                 label_cell_groups = F, label_groups_by_cluster = T,
                 label_branch_points = F, label_roots = F, label_leaves = F,
                 label_principal_points = F) + ggtitle("Monocle_cluster")+
  scale_color_manual(values = col)

r3 <- plot_cells(cds, reduction_method = "UMAP", graph_label_size=10, show_trajectory_graph = F,
                 color_cells_by = "moncuster", 
                 label_cell_groups = F, label_groups_by_cluster = T,
                 label_branch_points = F, label_roots = F, label_leaves = F,
                 label_principal_points = F) + ggtitle("moncuster_cluster")+
  scale_color_manual(values = col)



pdf("20210807_s4hV6Sub_Monocle3_monMN.pdf",width = 27,height = 6)
print(r1 + r2 +r3)
dev.off()

### learn_graph
cds <- learn_graph(cds)

pdf("20210807_s4hV6Sub_Monocle3_monMN.graph.pdf")

plot_cells(cds, reduction_method = "UMAP", graph_label_size=5, show_trajectory_graph = F,
           color_cells_by = "cluster", 
           label_cell_groups = T, label_groups_by_cluster = T,
           label_branch_points = F, label_roots = F, label_leaves = F,
           label_principal_points = F) + ggtitle("No_Label")+
  scale_color_manual(values = col) + coord_fixed()

plot_cells(cds, reduction_method = "UMAP", graph_label_size=5, show_trajectory_graph = F,
           color_cells_by = "partition", 
           label_cell_groups = T, label_groups_by_cluster = T,
           label_branch_points = F, label_roots = F, label_leaves = F,
           label_principal_points = F) + ggtitle("Partition")+
  scale_color_manual(values = col) + coord_fixed()

plot_cells(cds, reduction_method = "UMAP", graph_label_size=3, 
           color_cells_by = "cluster", 
           label_cell_groups = T, label_groups_by_cluster = F,
           label_branch_points = F, label_roots = T, label_leaves = F,
           label_principal_points = F) + ggtitle("Label_roots")+
  scale_color_manual(values = col) + coord_fixed()

plot_cells(cds, reduction_method = "UMAP", graph_label_size=3, 
           color_cells_by = "cluster", 
           label_cell_groups = T, label_groups_by_cluster = F,
           label_branch_points = F, label_roots = F, label_leaves = T,
           label_principal_points = F) + ggtitle("Label_leaves")+
  scale_color_manual(values = col) + coord_fixed()

plot_cells(cds, reduction_method = "UMAP", graph_label_size=3, 
           color_cells_by = "cluster", 
           label_cell_groups = T, label_groups_by_cluster = F,
           label_branch_points = T, label_roots = F, label_leaves = F,
           label_principal_points = F) + ggtitle("Label_branch_points")+
  scale_color_manual(values = col) + coord_fixed()

plot_cells(cds, reduction_method = "UMAP", graph_label_size=2, 
           color_cells_by = "cluster", 
           label_cell_groups = F, label_groups_by_cluster = T,
           label_branch_points = F, label_roots = F, label_leaves = F,
           label_principal_points = T ) +   ggtitle("label_principal_points")+
  scale_color_manual(values = col) + coord_fixed()

dev.off()

cds <- order_cells(cds, root_pr_nodes= c("Y_25"))
pse1 <- plot_cells(cds,trajectory_graph_segment_size = 0.25,
                   color_cells_by = "pseudotime",
                   label_cell_groups=FALSE,
                   label_leaves=FALSE,
                   label_branch_points=FALSE,
                   graph_label_size=1.5)+ggtitle("Root Y_25")

cds <- order_cells(cds, root_pr_nodes= c("Y_25","Y_49"))
pse2 <- plot_cells(cds, trajectory_graph_segment_size = 0.25,
                   color_cells_by = "pseudotime",
                   label_cell_groups=FALSE,
                   label_leaves=FALSE,
                   label_branch_points=FALSE,
                   graph_label_size=1.5)+ggtitle("Root Y_25 Y_49")

pse <-  pse1+pse2
ggsave("20210807_s4hV6Sub_Monocle3_monPN.pseudotime.pdf",pse,width = 18,height = 6)
### plot marker gene
marker <- read.table("marker1.txt",header = F)[,1]
p <- plot_cells(cds, genes = marker, 
                show_trajectory_graph=FALSE, 
                label_cell_groups=FALSE)+
  scale_color_viridis(option = "C")+ coord_fixed()
ggsave("20210807_s4hV6Sub_Monocle3_monMN.marker1.pdf",p,height=32,width = 32)

marker <- read.table("marker2.txt",header = F)[,1]
p<-plot_cells(cds, genes = marker, 
              show_trajectory_graph=FALSE, 
              label_cell_groups=FALSE)+
  scale_color_viridis(option = "C")+ coord_fixed()
ggsave("20210807_s4hV6Sub_Monocle3_monMN.marker2.pdf",p,height=24,width = 32)

marker <- read.table("marker3.txt",header = F)[,1]
p<-plot_cells(cds, genes = unique(marker), 
              show_trajectory_graph=FALSE, 
              label_cell_groups=FALSE)+
  scale_color_viridis(option = "C")+ coord_fixed()
ggsave("20210807_s4hV6Sub_Monocle3_monMN.marker3.pdf",p,height=48,width = 56,limitsize = FALSE)

marker <- read.table("marker4.txt",header = F)[,1]
p<-plot_cells(cds, genes = unique(marker), 
              show_trajectory_graph=FALSE, 
              label_cell_groups=FALSE)+
  scale_color_viridis(option = "C")+ coord_fixed()
ggsave("20210807_s4hV6Sub_Monocle3_monMN.marker4.pdf",p,height=32,width = 32)


# subset monocle cluster 4,29,19,13,14,17,1,23,32,15,12,18,11,8,3,2,35,6,34,10,36
moncluster <- cds@clusters$UMAP$clusters
sub_id1 <- names(cds@clusters$UMAP$clusters[cds@clusters$UMAP$clusters 
                                            %in% c(4,29,19,13,14,17,1,23,32,15,12,18,11,8,3,2,35,6,34,10,36)])
length(sub_id1)

sub_cluster1 <- moncluster[moncluster %in% c(4,29,19,13,14,17,1,23,32,15,12,18,11,8,3,2,35,6,34,10,36)]
length(sub_cluster1)

project.monpn <- subset(project, cells = sub_id1)
project.monpn

##### import Seurat data to Monocle
data <- GetAssayData(object = project.monpn, assay ="RNA", slot = "counts")
celldata <- as.data.frame(project.monpn@meta.data)
genedata <- as.data.frame(x = row.names(data), row.names = row.names(data))
colnames(genedata) <- "gene_short_name"
cds.monpn <- new_cell_data_set(data, cell_metadata = celldata, gene_metadata = genedata)

##### Estimate size factors and dispersions
cds.monpn <- estimate_size_factors(cds.monpn)
cds.monpn <- preprocess_cds(cds.monpn)
cds.monpn <- align_cds(cds.monpn, alignment_group = "Sample")
cds.monpn <- reduce_dimension(cds.monpn,reduction_method = "UMAP")

head(sub_cluster1)
tail(sub_cluster1)
colData(cds.monpn)$moncuster <- sub_cluster1

# ##### use Seurat umap
# cds.embed <- cds@int_colData$reducedDims$UMAP
# int.embed <- Embeddings(project, reduction = "umap")
# int.embed <- int.embed[rownames(cds.embed),]
# cds@int_colData$reducedDims$UMAP <- int.embed

cds.monpn <- cluster_cells(cds.monpn, verbose = T)

r1 <- plot_cells(cds.monpn, reduction_method = "UMAP", graph_label_size=10, show_trajectory_graph = F,
                 color_cells_by = "seurat_clusters", 
                 label_cell_groups = T, label_groups_by_cluster = T,
                 label_branch_points = F, label_roots = F, label_leaves = F,
                 label_principal_points = F) + ggtitle("Seurat_cluster")+
  scale_color_manual(values = col)
r1$theme$legend.position <- "right"

r2 <- plot_cells(cds.monpn, reduction_method = "UMAP", graph_label_size=10, show_trajectory_graph = F,
                 color_cells_by = "cluster", 
                 label_cell_groups = F, label_groups_by_cluster = T,
                 label_branch_points = F, label_roots = F, label_leaves = F,
                 label_principal_points = F) + ggtitle("Monocle_cluster")+
  scale_color_manual(values = col)

r3 <- plot_cells(cds.monpn, reduction_method = "UMAP", graph_label_size=10, show_trajectory_graph = F,
                 color_cells_by = "moncuster", 
                 label_cell_groups = T, label_groups_by_cluster = T,
                 label_branch_points = F, label_roots = F, label_leaves = F,
                 label_principal_points = F) + ggtitle("moncuster_cluster")+
  scale_color_manual(values = col)
r3$theme$legend.position <- "right"


pdf("20210807_s4hV6Sub_Monocle3_monPN.pdf",width = 27,height = 6)
print(r1 + r2 +r3)
dev.off()

### learn_graph
cds.monpn <- learn_graph(cds.monpn)

pdf("20210807_s4hV6Sub_Monocle3_monPN.graph.pdf")

plot_cells(cds.monpn, reduction_method = "UMAP", graph_label_size=5, show_trajectory_graph = F,
           color_cells_by = "cluster", 
           label_cell_groups = T, label_groups_by_cluster = T,
           label_branch_points = F, label_roots = F, label_leaves = F,
           label_principal_points = F) + ggtitle("No_Label")+
  scale_color_manual(values = col) + coord_fixed()

plot_cells(cds.monpn, reduction_method = "UMAP", graph_label_size=5, show_trajectory_graph = F,
           color_cells_by = "partition", 
           label_cell_groups = T, label_groups_by_cluster = T,
           label_branch_points = F, label_roots = F, label_leaves = F,
           label_principal_points = F) + ggtitle("Partition")+
  scale_color_manual(values = col) + coord_fixed()

plot_cells(cds.monpn, reduction_method = "UMAP", graph_label_size=3, 
           color_cells_by = "cluster", 
           label_cell_groups = T, label_groups_by_cluster = F,
           label_branch_points = F, label_roots = T, label_leaves = F,
           label_principal_points = F) + ggtitle("Label_roots")+
  scale_color_manual(values = col) + coord_fixed()

plot_cells(cds.monpn, reduction_method = "UMAP", graph_label_size=3, 
           color_cells_by = "cluster", 
           label_cell_groups = T, label_groups_by_cluster = F,
           label_branch_points = F, label_roots = F, label_leaves = T,
           label_principal_points = F) + ggtitle("Label_leaves")+
  scale_color_manual(values = col) + coord_fixed()

plot_cells(cds.monpn, reduction_method = "UMAP", graph_label_size=3, 
           color_cells_by = "cluster", 
           label_cell_groups = T, label_groups_by_cluster = F,
           label_branch_points = T, label_roots = F, label_leaves = F,
           label_principal_points = F) + ggtitle("Label_branch_points")+
  scale_color_manual(values = col) + coord_fixed()

plot_cells(cds.monpn, reduction_method = "UMAP", graph_label_size=2, 
           color_cells_by = "cluster", 
           label_cell_groups = F, label_groups_by_cluster = T,
           label_branch_points = F, label_roots = F, label_leaves = F,
           label_principal_points = T ) +   ggtitle("label_principal_points")+
  scale_color_manual(values = col) + coord_fixed()

dev.off()

cds.monpn <- order_cells(cds.monpn, root_pr_nodes= c("Y_25"))
pse1 <- plot_cells(cds.monpn,trajectory_graph_segment_size = 0.25,
                   color_cells_by = "pseudotime",
                   label_cell_groups=FALSE,
                   label_leaves=FALSE,
                   label_branch_points=FALSE,
                   graph_label_size=1.5)+ggtitle("Root Y_25")

cds.monpn <- order_cells(cds.monpn, root_pr_nodes= c("Y_25","Y_49"))
pse2 <- plot_cells(cds.monpn, trajectory_graph_segment_size = 0.25,
                   color_cells_by = "pseudotime",
                   label_cell_groups=FALSE,
                   label_leaves=FALSE,
                   label_branch_points=FALSE,
                   graph_label_size=1.5)+ggtitle("Root Y_25 Y_49")

pse <-  pse1+pse2
ggsave("20210807_s4hV6Sub_Monocle3_monPN.pseudotime.pdf",pse,width = 18,height = 6)
### plot marker gene
marker <- read.table("marker1.txt",header = F)[,1]
p <- plot_cells(cds, genes = marker, 
                show_trajectory_graph=FALSE, 
                label_cell_groups=FALSE)+
  scale_color_viridis(option = "C")+ coord_fixed()
ggsave("20210807_s4hV6Sub_Monocle3_monPN.marker1.pdf",p,height=32,width = 32)

marker <- read.table("marker2.txt",header = F)[,1]
p<-plot_cells(cds, genes = marker, 
              show_trajectory_graph=FALSE, 
              label_cell_groups=FALSE)+
  scale_color_viridis(option = "C")+ coord_fixed()
ggsave("20210807_s4hV6Sub_Monocle3_monPN.marker2.pdf",p,height=24,width = 32)

marker <- read.table("marker3.txt",header = F)[,1]
p<-plot_cells(cds, genes = unique(marker), 
              show_trajectory_graph=FALSE, 
              label_cell_groups=FALSE)+
  scale_color_viridis(option = "C")+ coord_fixed()
ggsave("20210807_s4hV6Sub_Monocle3_monPN.marker3.pdf",p,height=48,width = 56,limitsize = FALSE)

marker <- read.table("marker4.txt",header = F)[,1]
p<-plot_cells(cds, genes = unique(marker), 
              show_trajectory_graph=FALSE, 
              label_cell_groups=FALSE)+
  scale_color_viridis(option = "C")+ coord_fixed()
ggsave("20210807_s4hV6Sub_Monocle3_monPN.marker4.pdf",p,height=32,width = 32)

save.image("20210810.s4hV6Sub_Monocle3_monPN.RData")
