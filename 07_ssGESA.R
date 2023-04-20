###########################Do ssGESA 
library(tibble)
library(magrittr)
library(corto)
library(GSVA)
library(limma)
rm(list=ls())

data<-read.table(file='GSCs.tpm.150',sep='\t',header=T,row.names=1)
minRowFPKM=rowMeans(data)>2 # filter by mean
minNumFPKM=rowSums(data>0)>10 #Screen the number of samples whose expression level is not 0
data=data[minRowFPKM & minNumFPKM,] #联合一下
data=log2(data+1)
dat <- as.matrix(data)
#ssgsea1<- gsva(dat, l,method='ssgsea')
Z=ssgsea(dat,l,scale = TRUE)
b=z2p(Z)
-log10(b)
write.table(-log10(b), "ssGESA-logp.out", col.names = TRUE, row.names = TRUE, sep="\t", quote = FALSE)
write.table(b, "ssGESA_p.out", col.names = TRUE, row.names = TRUE, sep="\t", quote = FALSE)
write.table(Z, "ssGESA_Z.out", col.names = TRUE, row.names = TRUE, sep="\t", quote = FALSE)
##########pLOT##################
library(pheatmap)
library("RColorBrewer")

pdf("Z.pdf",width=10,height=3)
hmcol <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
hmcol <- hmcol[length(hmcol):1]
pheatmap(ssgsea1,cellwidth=13,cellheight=15
         col=hmcol,scale="row",
         show_colnames = T,
         cluster_rows = F,
         cluster_cols = T,
         fontsize=15)
dev.off()
pdf("Z_sort.pdf")
a=read.table("ssGESA_Z_sort.out",header=T,row.names=1)
breaksList = seq(-4, 4, by = 0.1)
cols=colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(breaksList))
pheatmap(a,color= cols, cellwidth=13,cellheight=15,breaks = breaksList, angle_col = 45,cluster_rows = FALSE,  cluster_cols = FALSE,)
dev.off()
pdf("-log_sort.pdf")
a=read.table("ssGESA-logp_sort.out",header=T,row.names=1)
breaksList = seq(0, 5, by = 0.1)
cols=colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(breaksList))
pheatmap(a,color= cols, cellwidth=13,cellheight=15,breaks = breaksList, angle_col = 45,cluster_rows = FALSE,  cluster_cols = FALSE,)
dev.off()

#########################plot PCA


library("FactoMineR")
library("factoextra")
myfpkm<-read.table("merge.allgene.tpm",header=T,comment.char="",sep = "\t",check.names=FALSE,row.names=1)
minRowFPKM=rowMeans(myfpkm)>2   # filter by mean
minNumFPKM=rowSums(myfpkm>0)>3#Screen the number of samples whose expression level is not 0
myfpkm=myfpkm[minRowFPKM & minNumFPKM,] #联合
logmyfpkm=log2(myfpkm+1)
pdf('PCA.pdf')
 res.pca <- PCA(t(logmyfpkm), graph = FALSE)
fviz_pca_ind(res.pca, col.ind = "cos2", 
              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
              repel = TRUE # Avoid text overlapping (slow if many points)
              )
dev.off()
