##set working directory
setwd('~/todd/Data/Homer_heatmaps')

library(ComplexHeatmap,lib= '~/Rpackages/')
library(circlize,lib= '~/Rpackages/')


##Set bed file for heatmaps

bed="~/todd/Data/Output/ROI_coord.txt"
 

##set directories for homerTag directories

tag.dir = '~/chris/Homer_Tags/'
tag.list = c("ES_ATAC","ES_H3K27ac","ES_H3K4me1","ES_Nanog","ES_Oct4","ES_Sox2","TS_ATAC","TS_Cdx2","TS_Elf5","TS_Eomes","TS_H3K27ac","TS_H3K4me1")

##Set colours you want to plot each tag in

col.list =c("red","red","red","red","red","red",
            "blue","blue","blue","blue","blue","blue")

##Set which classes you want to look at and which ones are relevant for each cell type

class.list<-c("ESC_REDE","ESC_noclass","ESC_NonEnh","ESC_NEDE","TSC_REDE","TSC_noclass","TSC_NonEnh","TSC_NEDE")

order.dir<-'~/todd/Data/Homer_heatmaps/'
tags.for.order<-c("ES_H3K27ac","ES_H3K27ac","ES_H3K27ac","ES_H3K27ac","ES_H3K27ac","ES_H3K27ac",
                  "TS_H3K27ac","TS_H3K27ac","TS_H3K27ac","TS_H3K27ac","TS_H3K27ac","TS_H3K27ac")
                  
group.1<-c("ESC_REDE","ESC_noclass","ESC_NonEnh","TSC_REDE","TSC_noclass","TSC_NonEnh")
group.2<-c("ESC_NEDE","TSC_NEDE")


##getting order for heatmaps for each class
lookup<-read.delim(bed)[,c(4,9)]


####Gets the order for each class based on the tags you have selected#### 

for(i in 1:length(class.list)){
  lookup.class<-as.character(lookup[grep(lookup[,2],pattern = class.list[i]),1])
  mat<-read.delim(paste(order.dir,tags.for.order[i],"_hist.txt",sep=""),h=T)
  mat[,1]<-as.character(mat[,1])
  mat<-mat[mat[,1]%in% lookup.class,]
  total = rowSums(mat[2:ncol(mat)])
  peak.order = mat[order(total,decreasing=T),1]
  df<-data.frame(peak.order,seq(1:length(peak.order)))
  #writes to a file for later use
  write.table(df,paste(class.list[i],"_heatmap_order.txt",sep=""),sep="\t",quote = F,col.names = F,row.names = F)}


####Plot complex heatmap ####

##Miguel's heatmap function
draw.heatmap = function(peak.class,homer.mat,order.IDs,name,width=1.5,height=5,col) {
  
  ##peak IDs
  peaks = read.delim(order.IDs,header=F,as.is=T)[,1]
  ##open heatmap and add 0 value to missing lines
  hm = read.delim(homer.mat)
  no.data = peaks[!(peaks %in% hm[,1])]
  ##seperate the names and matrix values of heatmap to plot
  hm.nam = c(as.character(hm[,1]),as.character(no.data))
  if(length(no.data>0)){hm.val = rbind(as.matrix(hm[,-1]),matrix(rep(0,length(no.data)*(ncol(hm)-1)),nrow=length(no.data)))}
     if(length(no.data)==0){hm.val=as.matrix(hm[,-1])}
  ##set saturation point
  
    if(peak.class %in% group.1){
    group.1.id<-lookup[grep(x=lookup[,2],pattern=paste(group.1,collapse="|")),1]
    sat = quantile(hm.val[match(group.1.id,hm.nam)],0.99,na.rm=T)
    }
    if(peak.class %in% group.2){
      group.2.id<-lookup[grep(x=lookup[,2],pattern=paste(group.2,collapse="|")),1]
      sat = quantile(hm.val[match(group.2.id,hm.nam)],0.95,na.rm=T)
      if(sat==0){sat=quantile(hm.val[match(group.2.id,hm.nam)],0.99,na.rm=T)}
  }
  ##Subset heatmap into the specific group
  peak.order =read.delim(order.IDs,h=F)[,1]
  mat = hm.val[match(peak.order,hm.nam),]
  ##make heatmap object
  colramp = colorRamp2(c(min(mat,na.rm=T),sat),c('white',col))
  hm.plot= Heatmap(mat,cluster_rows=F,cluster_columns=F,
                    show_row_names=F,show_column_names=F,show_heatmap_legend=F,col=colramp)
  jpeg(paste(name,'heatmap.jpg',sep='_'),w=width,h=height,units='in',res=600)
  draw(hm.plot)
  dev.off()
  }

for(m in 1:length(class.list)){
    #pdf(paste(class.list[m],"_heatmap_plots.pdf",sep=""),width = 3, height = 5)
    for(i in 1:length(tag.list)){
    draw.heatmap(peak.class=class.list[m],
                 order.IDs=paste(class.list[m],"_heatmap_order.txt",sep=""),
                 homer.mat=paste(tag.list[i],"_hist.txt",sep=""),
                 name=paste(class.list[m],tag.list[i],sep="_"),
                 col=col.list[i])
    }
    }
    
