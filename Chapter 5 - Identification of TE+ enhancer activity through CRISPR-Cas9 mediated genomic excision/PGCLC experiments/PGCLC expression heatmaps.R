mat=read.delim("~/R/Bioinformatic for Thesis/PGC_heatmap2.txt")
head(mat)
val.mat=as.matrix(mat[,2:ncol(mat)])
row.names(val.mat)=mat[,1]
library(gplots)
heatmap.2(val.mat,trace = "none",dendrogram = "none",Colv = F,
          Rowv = F,scale = "column",col=bluered,key=F,
          margins =c(14,14))

heatmap.2(val.mat,trace = "none",dendrogram = "none",Colv = F,Rowv = F,
         col=bluered,key=T,density.info = "none",
          margins =c(14,14))

sample.mat=read.delim("~/R/Bioinformatic for Thesis/PGC_Tdrd12_heatmap3.txt")
samp.val=as.matrix(sample.mat[,2:ncol(sample.mat)])
row.names(samp.val)=sample.mat$Sample

head(val.mat)

heatmap.2(samp.val,trace = "none",dendrogram = "none",Colv = F,Rowv = F,
          col=bluered,key=T,density.info = "none",
          margins =c(14,14))
