load("20210924.s4hV6Sub.monpn.pgsplit.RData")
library(Seurat)
library(pheatmap)
##### 20211026 PGsplit marker
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

# marker1
marker <- read.table("20211026_PGsplit_marker1.txt", header = F)

plotdata <- cluster_mean.pgsplit[,colnames(cluster_mean.pgsplit) %in% marker$V1]
plotdata <- t(plotdata)

colnames(plotdata) <- c(1:8)
plotdata <- plotdata[match(marker$V1,rownames(plotdata)),match(c(1:8),colnames(plotdata))]

pdf("20211026.s4hV6Sub.PGsplit.marker1.heatmap.pdf", height=10)
	pheatmap(plotdata, cluster_cols = F, cluster_rows = T, scale = "row")
dev.off()


# marker2
marker <- read.table("20211026_PGsplit_marker2.txt", header = F)

plotdata <- cluster_mean.pgsplit[,colnames(cluster_mean.pgsplit) %in% marker$V1]
plotdata <- t(plotdata)

colnames(plotdata) <- c(1:8)
plotdata <- plotdata[match(marker$V1,rownames(plotdata)),match(c(1:8),colnames(plotdata))]

pdf("20211026.s4hV6Sub.PGsplit.marker2.heatmap.pdf", height=20)
	pheatmap(plotdata, cluster_cols = F, cluster_rows = T, scale = "row")
dev.off()

##### 20211130
marker <- read.table("20211130.IN.marker1.txt",sep = "\t")

pdf("20211130.s4hV6Sub.PGsplit.marker1.pdf")
  DotPlot(project.pgsplit,features = unique(marker$V1)) + coord_flip()
dev.off()

plot.data <- project.pgsplit@assays$RNA@scale.data[rownames(project.pgsplit@assays$RNA@scale.data) %in% marker$V1,]
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

color <- col[c(1:8)][c(3,4,5,7,8,2,1,6)]

p<-ggplot(plot.data.final)+
  geom_violin(aes(x=expression,y=cluster,fill = cluster),scale="width")+
  theme_bw()+ ylab("")+xlab("")+
  facet_grid(.~symbol) +
  theme(plot.margin=unit(x=c(0.05,0,0,0),units="inches"),strip.background = element_blank()) +
  theme(strip.text.x = element_text(angle = 90)) +
  theme(axis.line.y = element_line(color = "black"),axis.title.y=element_text(size=8)) +
  theme(axis.line.x = element_blank(),axis.ticks.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid=element_blank()) +
  scale_fill_manual(values = color)
p

ggsave("20211130.PG.miniviolin.pdf",p,width = 8)


plot.data.final2 <- plot.data.final[!(plot.data.final$symbol %in% c("DLX2","DLX5","DLX6-AS1","NKX2-1","LHX6")),]
p<-ggplot(plot.data.final2)+
  geom_violin(aes(x=expression,y=cluster,fill = cluster),scale="width")+
  theme_bw()+ ylab("")+xlab("")+
  facet_grid(sample~symbol) +
  theme(plot.margin=unit(x=c(0.05,0,0,0),units="inches"),strip.background = element_blank()) +
  theme(strip.text.x = element_text(angle = 90)) +
  theme(axis.line.y = element_line(color = "black"),axis.title.y=element_text(size=8)) +
  theme(axis.line.x = element_blank(),axis.ticks.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid=element_blank()) +
  scale_fill_manual(values = color)
p
ggsave("20211130.PG.miniviolin.sample.pdf",p,width = 12,height = 8)

##### 20211201
# marker1
marker <- read.table("20211201.PG.marker1.txt",sep = "\t")

plot.data <- project.pgsplit@assays$RNA@scale.data[rownames(project.pgsplit@assays$RNA@scale.data) %in% marker$V1,]
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

color <- col[c(1:8)][c(3,4,5,7,8,2,1,6)]

p<-ggplot(plot.data.final)+
  geom_violin(aes(x=expression,y=cluster,fill = cluster),scale="width")+
  theme_bw()+ ylab("")+xlab("")+
  facet_grid(sample~symbol) +
  theme(plot.margin=unit(x=c(0.05,0,0,0),units="inches"),strip.background = element_blank()) +
  theme(strip.text.x = element_text(angle = 90)) +
  theme(axis.line.y = element_line(color = "black"),axis.title.y=element_text(size=8)) +
  theme(axis.line.x = element_blank(),axis.ticks.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid=element_blank()) +
  scale_fill_manual(values = color)
p

ggsave("20211201.PG.miniviolin1.sample.pdf",p,width = 12,height = 8)

p<-ggplot(plot.data.final)+
  geom_violin(aes(x=expression,y=cluster,fill = cluster),scale="width")+
  theme_bw()+ ylab("")+xlab("")+
  facet_grid(.~symbol) +
  theme(plot.margin=unit(x=c(0.05,0,0,0),units="inches"),strip.background = element_blank()) +
  theme(strip.text.x = element_text(angle = 90)) +
  theme(axis.line.y = element_line(color = "black"),axis.title.y=element_text(size=8)) +
  theme(axis.line.x = element_blank(),axis.ticks.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid=element_blank()) +
  scale_fill_manual(values = color)
p

ggsave("20211201.PG.miniviolin1.pdf",p,width = 12,height = 8)

# marker2
marker <- read.table("20211201.PG.marker2.txt",sep = "\t")

plot.data <- project.pgsplit@assays$RNA@scale.data[rownames(project.pgsplit@assays$RNA@scale.data) %in% marker$V1,]
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

color <- col[c(1:8)][c(3,4,5,7,8,2,1,6)]

p<-ggplot(plot.data.final)+
  geom_violin(aes(x=expression,y=cluster,fill = cluster),scale="width")+
  theme_bw()+ ylab("")+xlab("")+
  facet_grid(sample~symbol) +
  theme(plot.margin=unit(x=c(0.05,0,0,0),units="inches"),strip.background = element_blank()) +
  theme(strip.text.x = element_text(angle = 90)) +
  theme(axis.line.y = element_line(color = "black"),axis.title.y=element_text(size=8)) +
  theme(axis.line.x = element_blank(),axis.ticks.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid=element_blank()) +
  scale_fill_manual(values = color)
p

ggsave("20211201.PG.miniviolin2.sample.pdf",p,width = 20,height = 8)

p<-ggplot(plot.data.final)+
  geom_violin(aes(x=expression,y=cluster,fill = cluster),scale="width")+
  theme_bw()+ ylab("")+xlab("")+
  facet_grid(.~symbol) +
  theme(plot.margin=unit(x=c(0.05,0,0,0),units="inches"),strip.background = element_blank()) +
  theme(strip.text.x = element_text(angle = 90)) +
  theme(axis.line.y = element_line(color = "black"),axis.title.y=element_text(size=8)) +
  theme(axis.line.x = element_blank(),axis.ticks.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid=element_blank()) +
  scale_fill_manual(values = color)
p

ggsave("20211201.PG.miniviolin2.pdf",p,width = 20,height = 8)

# marker3
marker <- read.table("20211201.PG.marker3.txt",sep = "\t")

plot.data <- project.pgsplit@assays$RNA@scale.data[rownames(project.pgsplit@assays$RNA@scale.data) %in% marker$V1,]
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
plot.data.final <- plot.data.final[plot.data.final$sample %in% c("SC8","SC9","SC12"),]

color <- col[c(1:8)][c(3,4,5,7,8,2,1,6)]

p<-ggplot(plot.data.final)+
  geom_violin(aes(x=expression,y=cluster,fill = cluster),scale="width")+
  theme_bw()+ ylab("")+xlab("")+
  facet_grid(sample~symbol) +
  theme(plot.margin=unit(x=c(0.05,0,0,0),units="inches"),strip.background = element_blank()) +
  theme(strip.text.x = element_text(angle = 90)) +
  theme(axis.line.y = element_line(color = "black"),axis.title.y=element_text(size=8)) +
  theme(axis.line.x = element_blank(),axis.ticks.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid=element_blank()) +
  scale_fill_manual(values = color)
p

ggsave("20211201.PG.miniviolin3.sample.pdf",p,width = 12,height = 8)

p<-ggplot(plot.data.final)+
  geom_violin(aes(x=expression,y=cluster,fill = cluster),scale="width")+
  theme_bw()+ ylab("")+xlab("")+
  facet_grid(.~symbol) +
  theme(plot.margin=unit(x=c(0.05,0,0,0),units="inches"),strip.background = element_blank()) +
  theme(strip.text.x = element_text(angle = 90)) +
  theme(axis.line.y = element_line(color = "black"),axis.title.y=element_text(size=8)) +
  theme(axis.line.x = element_blank(),axis.ticks.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid=element_blank()) +
  scale_fill_manual(values = color)
p

ggsave("20211201.PG.miniviolin3.pdf",p,width = 12,height = 8)


##### 20211209
marker <- read.table("20211209.marker1.txt",sep = "\t")

plot.data <- project.pgsplit@assays$RNA@scale.data[rownames(project.pgsplit@assays$RNA@scale.data) %in% marker$V1,]
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


plot.data.final1 <- plot.data.final[plot.data.final$sample %in% c("SC8","SC12"),]
plot.data.final1 <- plot.data.final1[plot.data.final1$cluster %in% c(3,4,5,7,8),]

color <- col[c(1:8)][c(3,4,5,7,8,2,1,6)]

p<-ggplot(plot.data.final1)+
  geom_violin(aes(x=expression,y=cluster,fill = cluster),scale="width")+
  theme_bw()+ ylab("")+xlab("")+
  facet_grid(sample~symbol) +
  theme(plot.margin=unit(x=c(0.05,0,0,0),units="inches"),strip.background = element_blank()) +
  theme(strip.text.x = element_text(angle = 90)) +
  theme(axis.line.y = element_line(color = "black"),axis.title.y=element_text(size=8)) +
  theme(axis.line.x = element_blank(),axis.ticks.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid=element_blank()) +
  scale_fill_manual(values = color)
p

ggsave("20211209.PG.miniviolin1.sample.pdf",p,width = 6,height = 2.5)


plot.data.final2 <- plot.data.final[plot.data.final$cluster %in% c(3,4,5,7,8),]

p<-ggplot(plot.data.final2)+
  geom_violin(aes(x=expression,y=cluster,fill = cluster),scale="width")+
  theme_bw()+ ylab("")+xlab("")+
  facet_grid(sample~symbol) +
  theme(plot.margin=unit(x=c(0.05,0,0,0),units="inches"),strip.background = element_blank()) +
  theme(strip.text.x = element_text(angle = 90)) +
  theme(axis.line.y = element_line(color = "black"),axis.title.y=element_text(size=8)) +
  theme(axis.line.x = element_blank(),axis.ticks.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid=element_blank()) +
  scale_fill_manual(values = color)
p

ggsave("20211209.PG.miniviolin1.pdf",p,width = 6,height = 4)

##### 20220107
load("/data2/chenzixi/20210419_WJ_SC/20210803_s4hV6Sub_Monocle3/20211111.s4hV6.all.Rdata")
# marker1
marker <- read.table("20220107.PG.marker1.txt",sep = "\t")

plot.data <- project.pgsplit@assays$RNA@scale.data[rownames(project.pgsplit@assays$RNA@scale.data) %in% marker$V1,]
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


plot.data.final1 <- plot.data.final[plot.data.final$sample %in% c("SC8","SC12"),]
plot.data.final1 <- plot.data.final1[plot.data.final1$cluster %in% c(3,4,5,7,8),]

color <- col[c(1:8)][c(3,4,5,7,8,2,1,6)]

p1 <- ggplot(plot.data.final1)+
  geom_violin(aes(x=expression,y=cluster,fill = cluster),scale="width")+
  theme_bw()+ ylab("")+xlab("")+
  facet_grid(sample~symbol) +
  theme(plot.margin=unit(x=c(0.05,0,0,0),units="inches"),strip.background = element_blank()) +
  theme(strip.text.x = element_text(angle = 90)) +
  theme(axis.line.y = element_line(color = "black"),axis.title.y=element_text(size=8)) +
  theme(axis.line.x = element_blank(),axis.ticks.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid=element_blank()) +
  scale_fill_manual(values = color)
p1

ggsave("20220107.PG.miniviolin1.sample.pdf",p1,width = 20,height = 8 )

plot.data.final2 <- plot.data.final[plot.data.final$cluster %in% c(3,4,5,7,8),]

p2 <-ggplot(plot.data.final2)+
  geom_violin(aes(x=expression,y=cluster,fill = cluster),scale="width")+
  theme_bw()+ ylab("")+xlab("")+
  facet_grid(sample~symbol) +
  theme(plot.margin=unit(x=c(0.05,0,0,0),units="inches"),strip.background = element_blank()) +
  theme(strip.text.x = element_text(angle = 90)) +
  theme(axis.line.y = element_line(color = "black"),axis.title.y=element_text(size=8)) +
  theme(axis.line.x = element_blank(),axis.ticks.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid=element_blank()) +
  scale_fill_manual(values = color)
p2

ggsave("20220107.PG.miniviolin1.pdf",p2,width = 20,height = 8 )

# marker2
marker <- read.table("20220107.PG.marker2.txt",sep = "\t")

plot.data <- project.pgsplit@assays$RNA@scale.data[rownames(project.pgsplit@assays$RNA@scale.data) %in% marker$V1,]
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


plot.data.final1 <- plot.data.final[plot.data.final$sample %in% c("SC8","SC12"),]
plot.data.final1 <- plot.data.final1[plot.data.final1$cluster %in% c(3,4,5,7,8),]

color <- col[c(1:8)][c(3,4,5,7,8,2,1,6)]

p1 <- ggplot(plot.data.final1)+
  geom_violin(aes(x=expression,y=cluster,fill = cluster),scale="width")+
  theme_bw()+ ylab("")+xlab("")+
  facet_grid(sample~symbol) +
  theme(plot.margin=unit(x=c(0.05,0,0,0),units="inches"),strip.background = element_blank()) +
  theme(strip.text.x = element_text(angle = 90)) +
  theme(axis.line.y = element_line(color = "black"),axis.title.y=element_text(size=8)) +
  theme(axis.line.x = element_blank(),axis.ticks.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid=element_blank()) +
  scale_fill_manual(values = color)
#p1

ggsave("20220107.PG.miniviolin2.sample.pdf",p1,width = 30,height = 8 )

plot.data.final2 <- plot.data.final[plot.data.final$cluster %in% c(3,4,5,7,8),]

p2 <-ggplot(plot.data.final2)+
  geom_violin(aes(x=expression,y=cluster,fill = cluster),scale="width")+
  theme_bw()+ ylab("")+xlab("")+
  facet_grid(sample~symbol) +
  theme(plot.margin=unit(x=c(0.05,0,0,0),units="inches"),strip.background = element_blank()) +
  theme(strip.text.x = element_text(angle = 90)) +
  theme(axis.line.y = element_line(color = "black"),axis.title.y=element_text(size=8)) +
  theme(axis.line.x = element_blank(),axis.ticks.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid=element_blank()) +
  scale_fill_manual(values = color)
#p2

ggsave("20220107.PG.miniviolin2.pdf",p2,width = 30,height = 8 )

### feature plot pMN
marker <- read.table("20220107.PG.pMN.marker1.txt", header = F, sep = "\t")

pdf("20220107.PG.pMN.s4hV6Sub.PGsplit.FeaturePlot.pdf", width = 24, height = ceiling(length(marker$V1)))
Seurat::FeaturePlot(project.pgsplit, reduction = "umap",
                    features = as.vector(unique(marker$V1)),
                    keep.scale = "all", ncol = 4,
                    cols = c("lightgrey", "red"))
dev.off()

### feature plot pIN
marker <- read.table("20220107.PG.pIN.marker1.txt", header = F, sep = "\t")

pdf("20220107.PG.pIN.s4hV6Sub.PGsplit.FeaturePlot.pdf", width = 24, height = ceiling(length(marker$V1)))
Seurat::FeaturePlot(project.pgsplit, reduction = "umap",
                    features = as.vector(unique(marker$V1)),
                    keep.scale = "all", ncol = 4,
                    cols = c("lightgrey", "red"))
dev.off()

##### 20220109
load("/data2/chenzixi/20210419_WJ_SC/20210803_s4hV6Sub_Monocle3/20211111.s4hV6.all.Rdata")
# marker1
marker <- read.table("20220107.PG.marker1.txt",sep = "\t")

plot.data <- project.pgsplit@assays$RNA@scale.data[rownames(project.pgsplit@assays$RNA@scale.data) %in% marker$V1,]
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

p2 <-ggplot(plot.data.final2)+
  geom_violin(aes(x=expression,y=cluster,fill = cluster),scale="width")+
  theme_bw()+ ylab("")+xlab("")+
  facet_grid(.~symbol) +
  theme(plot.margin=unit(x=c(0.05,0,0,0),units="inches"),strip.background = element_blank()) +
  theme(strip.text.x = element_text(angle = 90)) +
  theme(axis.line.y = element_line(color = "black"),axis.title.y=element_text(size=8)) +
  theme(axis.line.x = element_blank(),axis.ticks.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid=element_blank()) +
  scale_fill_manual(values = color)
p2

ggsave("20220109.PG.miniviolin1.pdf",p2,width = 40,height = 8 )

# marker2
marker <- read.table("20220107.PG.marker2.txt",sep = "\t")

plot.data <- project.pgsplit@assays$RNA@scale.data[rownames(project.pgsplit@assays$RNA@scale.data) %in% marker$V1,]
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

p2 <-ggplot(plot.data.final2)+
  geom_violin(aes(x=expression,y=cluster,fill = cluster),scale="width")+
  theme_bw()+ ylab("")+xlab("")+
  facet_grid(.~symbol) +
  theme(plot.margin=unit(x=c(0.05,0,0,0),units="inches"),strip.background = element_blank()) +
  theme(strip.text.x = element_text(angle = 90)) +
  theme(axis.line.y = element_line(color = "black"),axis.title.y=element_text(size=8)) +
  theme(axis.line.x = element_blank(),axis.ticks.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid=element_blank()) +
  scale_fill_manual(values = color)
#p2

ggsave("20220109.PG.miniviolin2.pdf",p2,width = 60,height = 8 ,limitsize = F)

##### 20220109 PGsplit
library(monocle3)

marker <- read.table("20220107.PG.pMN.marker1.txt", header = F, sep = "\t")
p<-plot_cells(cds.pgsplit, genes=marker$V1, 
           show_trajectory_graph=FALSE, 
           label_cell_groups=FALSE) + 
  scale_color_gradient2(low="lightgrey",mid = "grey", high="red", breaks = c(0,0.5,1,1.5,2))


ggsave("20220109.PG.pMN.feature.pdf",p,width = 35,height = 21 )

marker <- read.table("20220107.PG.pIN.marker1.txt", header = F, sep = "\t")
p<-plot_cells(cds.pgsplit, genes=marker$V1, 
              show_trajectory_graph=FALSE, 
              label_cell_groups=FALSE) + 
  scale_color_gradient2(low="lightgrey",mid = "grey", high="red", breaks = c(0,0.5,1,1.5,2))
ggsave("20220109.PG.pIN.feature.pdf",p,width = 42,height = 21 )

##### 20220110 violin
marker <- read.table("20220110.marker1.txt",sep = "\t")

plot.data <- project.pgsplit@assays$RNA@scale.data[rownames(project.pgsplit@assays$RNA@scale.data) %in% marker$V1,]
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

p2 <-ggplot(plot.data.final2)+
  geom_violin(aes(x=expression,y=cluster,fill = cluster),scale="width")+
  theme_bw()+ ylab("")+xlab("")+
  facet_grid(.~symbol) +
  theme(plot.margin=unit(x=c(0.05,0,0,0),units="inches"),strip.background = element_blank()) +
  theme(strip.text.x = element_text(angle = 90)) +
  theme(axis.line.y = element_line(color = "black"),axis.title.y=element_text(size=8)) +
  theme(axis.line.x = element_blank(),axis.ticks.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid=element_blank()) +
  scale_fill_manual(values = color)
p2

ggsave("20220110.PG.miniviolin1.pdf",p2,width = 40,height = 8 )

# marker2
marker <- read.table("20220110.marker2.txt",sep = "\t")

plot.data <- project.pgsplit@assays$RNA@scale.data[rownames(project.pgsplit@assays$RNA@scale.data) %in% marker$V1,]
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

p2 <-ggplot(plot.data.final2)+
  geom_violin(aes(x=expression,y=cluster,fill = cluster),scale="width")+
  theme_bw()+ ylab("")+xlab("")+
  facet_grid(.~symbol) +
  theme(plot.margin=unit(x=c(0.05,0,0,0),units="inches"),strip.background = element_blank()) +
  theme(strip.text.x = element_text(angle = 90)) +
  theme(axis.line.y = element_line(color = "black"),axis.title.y=element_text(size=8)) +
  theme(axis.line.x = element_blank(),axis.ticks.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid=element_blank()) +
  scale_fill_manual(values = color)
p2

ggsave("20220110.PG.miniviolin2.pdf",p2,width = 20,height = 8 )

##### 20220111 violin
marker <- read.table("20220111.marker1.txt",sep = "\t", header = F)

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
  geom_violin(aes(x=expression,y=cluster,fill = cluster),scale="width")+
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

ggsave("20220111.PG.miniviolin1.pdf",p2,width = 24,height = 6, units = "cm" )

##### 20220111 enrichment
library(clusterProfiler)
library(org.Hs.eg.db)

data <- read.table("20220111.GO.1.txt", header = T, sep = "\t")
data$Type <- "pMN"
data$Type[data$Group %in% c("GO6","GO7","GO8","GO9","GO10")] <- "pIN"
list1 <- split(data$Gene, data$Group)
list2 <- split(data$Gene, data$Type)

go.gp <- compareCluster(list1,"enrichGO", OrgDb='org.Hs.eg.db', ont = "BP", keyType = "SYMBOL")
go.tp <- compareCluster(list2,"enrichGO", OrgDb='org.Hs.eg.db', ont = "BP", keyType = "SYMBOL")
go.mn <- enrichGO(data$Gene[1:128], OrgDb='org.Hs.eg.db', ont = "BP", keyType = "SYMBOL")
go.in <- enrichGO(data$Gene[129:431], OrgDb='org.Hs.eg.db', ont = "BP", keyType = "SYMBOL")

p1 <- clusterProfiler::dotplot(go.gp, showCategory = 100)
p1
p2 <- dotplot(go.tp)

p1 + p2 

write.table(go.gp, "20220111.PG.GO.xls", sep = "\t", quote = F, row.names = F)

dotplot(go.mn)
dotplot(go.in)

##### 20220112
library(RColorBrewer)
cols=c(brewer.pal(12,"Paired"),brewer.pal(8,"Dark2"),
       brewer.pal(9,"Set1"),brewer.pal(8,"Set2"),
       brewer.pal(12,"Set3"),brewer.pal(8,"Pastel2"),
       brewer.pal(9,"Pastel1"),brewer.pal(8,"Accent"))
col=unique(cols)
color <- col[c(1:8)][c(1,6,3,4,5,2,7,8)]

marker <- read.table("20220112.marker1.txt", header = F, sep = "\t")
cds.plot <- cds.pgsplit[rowData(cds.pgsplit)$gene_short_name %in% marker$V1,]
p<-monocle3::plot_genes_in_pseudotime(cds.plot, color_cells_by = "Seurat_res0.2", ncol = 4, cell_size = 1) + 
  scale_color_manual(values = color)
p
ggsave("20220112.s4hV6Sub.PGsplit.curve.pdf",p, width = 16, height = ceiling(length(marker$V1)/4)*2)

# marker2
marker <- read.table("20220112.marker2.txt", header = F, sep = "\t")

plotdata <- cluster_mean.pgsplit[,colnames(cluster_mean.pgsplit) %in% marker$V1]
plotdata <- t(plotdata)

colnames(plotdata) <- c(1:8)
plotdata <- plotdata[match(marker$V1,rownames(plotdata)),match(c(5,7,4,3,8,2,6,1),colnames(plotdata))]

pdf("20220112.s4hV6Sub.PGsplit.marker2.heatmap.pdf", height=10)
  pheatmap::pheatmap(plotdata, cluster_cols = F, cluster_rows = F, scale = "row", 
                                          angle_col = "0")
dev.off()

##### 20220113 pMN GO
marker <- read.table("20220113.pMN.GO.txt", header = F)
pMN.go <- enrichGO(marker$V1, OrgDb='org.Hs.eg.db', ont = "BP", keyType = "SYMBOL")
write.table(pMN.go, "20220113.pMN.go.xls", sep = "\t", quote = F, row.names = F)

marker <- read.table("20220113.pIN.GO.txt", header = F)
pIN.go <- enrichGO(marker$V1, OrgDb='org.Hs.eg.db', ont = "BP", keyType = "SYMBOL")
write.table(pIN.go, "20220113.pIN.go.xls", sep = "\t", quote = F, row.names = F)

##### 20220113 feature
marker <- read.table("20220113.PG.pMN.marker1.txt", header = F, sep = "\t")
p<-plot_cells(cds.pgsplit, genes=marker$V1, 
              show_trajectory_graph=FALSE, 
              label_cell_groups=FALSE) + 
  scale_color_gradient2(low="lightgrey",mid = "grey", high="red", breaks = c(0,0.5,1,1.5,2))
ggsave("20220113.PG.pMN.feature.pdf",p,width = 42,height = 21 )

marker <- read.table("20220113.PG.pIN.marker1.txt", header = F, sep = "\t")
p<-plot_cells(cds.pgsplit, genes=marker$V1, 
              show_trajectory_graph=FALSE, 
              label_cell_groups=FALSE) + 
  scale_color_gradient2(low="lightgrey",mid = "grey", high="red", breaks = c(0,0.5,1,1.5,2))
ggsave("20220113.PG.pIN.feature.pdf",p,width = 42,height = 21 )