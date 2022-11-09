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
  col                          = colorRampPalette(c("#0f86a9", "white", "#ed8b10"))(200),
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
  col                          = colorRampPalette(c("#0f86a9", "white", "#ed8b10"))(200),
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
  col                          = colorRampPalette(c("#0f86a9", "white", "#ed8b10"))(200),
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
  col                          = colorRampPalette(c("#0f86a9", "white", "#ed8b10"))(200),
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
  col                          = colorRampPalette(c("#0f86a9", "white", "#ed8b10"))(200),
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
  col                          = colorRampPalette(c("#0f86a9", "white", "#ed8b10"))(200),
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
  col                          = colorRampPalette(c("#0f86a9", "white", "#ed8b10"))(200),
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
  col                          = colorRampPalette(c("#0f86a9", "white", "#ed8b10"))(200),
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
  col                          = colorRampPalette(c("#0f86a9", "white", "#ed8b10"))(200),
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
  col                          = colorRampPalette(c("#0f86a9", "white", "#ed8b10"))(200),
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
  col                          = colorRampPalette(c("#0f86a9", "white", "#ed8b10"))(200),
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
  col                          = colorRampPalette(c("#0f86a9", "white", "#ed8b10"))(200),
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
  col                          = colorRampPalette(c("#0f86a9", "white", "#ed8b10"))(200),
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
  col                          = colorRampPalette(c("#0f86a9", "white", "#ed8b10"))(200),
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
  col                          = colorRampPalette(c("#0f86a9", "white", "#ed8b10"))(200),
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
  col                          = colorRampPalette(c("#0f86a9", "white", "#ed8b10"))(200),
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
  col                          = colorRampPalette(c("#0f86a9", "white", "#ed8b10"))(200),
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
  col                          = colorRampPalette(c("#0f86a9", "white", "#ed8b10"))(200),
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
  col                          = colorRampPalette(c("#0f86a9", "white", "#ed8b10"))(200),
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
  col                          = colorRampPalette(c("#0f86a9", "white", "#ed8b10"))(200),
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
  col                          = colorRampPalette(c("#0f86a9", "white", "#ed8b10"))(200),
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
  col                          = colorRampPalette(c("#0f86a9", "white", "#ed8b10"))(200),
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
  col                          = colorRampPalette(c("#0f86a9", "white", "#ed8b10"))(200),
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
  col                          = colorRampPalette(c("#0f86a9", "white", "#ed8b10"))(200),
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
