setwd("~/R/Thesis Scripts/Results Chapter 2/Prox vs Chic vs Goth/")

##TSS info for baitmaps in the HiC data
bait=read.delim("~/R/Thesis Scripts/Data/PCHi-C/bait_RefSeq_TSS_look.txt")

#ES Chicago Data
es.chic=read.delim("~/R/Thesis Scripts/Data/PCHi-C/Chicago_Inter_lookup.txt",h=F)
chic=read.delim("~/R/Thesis Scripts/Data/PCHi-C/ESC_Chicago.txt")[,c(4,8)]
#ES Gothic Data
es.goth=read.delim("~/R/Thesis Scripts/Data/PCHi-C/Gothic_Inter_lookup.txt",h=F)
hic='/Users/hmx579/Dropbox/Todd_data/Data/PCHi-C/gothic_interactions_all_xy'
#nearest gene Data
ROI=read.delim("~/R/Thesis Scripts/Data/Output/ROI_coord.txt")
#Refseq lookup file
refseq.lookup="~/R/Thesis Scripts/Data/Annotations/RefSeq_mRNA_lookup.txt"

##Set proximity model cutoff
cutoff=100000

##load Gothic data
hic.data<- local({
  load(hic)
  stopifnot(length(ls())==1)
  environment()[[ls()]]
}) 

es=hic.data$ESC
es$bait.id=paste(es$chr.bait,es$start.bait,es$end.bait,sep="-")
es$int.id<-paste(es$chr.target,es$start.target,es$end.target,sep="-")
es=unique(es[,11:12])



get.hic.distance=function(input.df){
  int1=strsplit(as.character(input.df[,1]),split = "-")
  int2=strsplit(as.character(input.df[,2]),split = "-")
  ori=as.numeric(unlist(lapply(int1, function(x){x[2]})))<as.numeric(unlist(lapply(int2, function(x){x[2]})))
  dist=c()
  for(n in 1:nrow(input.df)){
    if(ori[n]==T){dist[n]=as.numeric(unlist(lapply(int2, function(x){x[2]}))[n])-as.numeric(unlist(lapply(int1, function(x){x[3]}))[n])}
    if(ori[n]==F){dist[n]=as.numeric(unlist(lapply(int1, function(x){x[2]}))[n])-as.numeric(unlist(lapply(int2, function(x){x[3]}))[n])}
  }
  return(dist)
}
get.nearest.gene=function(input.group,chr.col,mid.col){
  refseq= read.delim(refseq.lookup)
  nearest.gene=c()
  distance=c()
  for(i in 1:nrow(input.group)){
    chr.sub=refseq[refseq$chrom %in% input.group[i,chr.col],]
    nearest=chr.sub[which.min(abs(input.group[i,mid.col]-chr.sub$TSS)),]
    distance[i]=abs(input.group[i,mid.col]-nearest$TSS)
    if(nearest$mRNA_name!="null"){nearest.gene[i]=as.character(nearest$mRNA_name)}else{nearest.gene[i]=as.character(nearest$refseq_name)}
  }
  df=data.frame(input.group,distance,nearest.gene)
  
  return(df)
}
make.same.length.vectors=function(input.list){
  len=max(unlist(lapply(input.list, function(x){length(x)})))
  new.list=(lapply(input.list, function(x){append(x,rep("",len-length(x)))}))
  return(new.list)
}

class.list=c("REDE","NEDE")

for(i in 1:length(class.list)){
###Get Chicago for class
class=es.chic[grep(es.chic[,2],pattern = class.list[i]),]
es.class=class[grep(class[,2],pattern = "ESC"),]
colnames(bait)=c("bait.id","Gene")
colnames(es.class)=c("int.id","ES_class")
chic.mer=merge(es.class,chic,by="int.id")
chic.bait.mer=merge(chic.mer,bait,by="bait.id")


###Get Gothic for class
class=es.goth[grep(es.goth[,2],pattern = class.list[i]),]
es.class=class[grep(class[,2],pattern = "ESC"),]
colnames(es.class)=c("int.id","ES_class")
got.mer=merge(es.class,es,by="int.id")
got.bait.mer=merge(got.mer,bait,by="bait.id")
head(got.bait.mer)


###Get nearest gene for class
class.roi=ROI[ROI$class==paste("ESC_",class.list[i],sep=""),]
class.roi$mid=rowMeans(class.roi[,7:8])
nearest.roi=get.nearest.gene(class.roi,1,10)
nearest.roi=nearest.roi[nearest.roi$distance<cutoff,]


chic.distance=get.hic.distance(chic.bait.mer[,1:2])
goth.distance=get.hic.distance(got.bait.mer[,1:2])
nearest.distance=nearest.roi$distance
dist.list=list(nearest.distance,chic.distance,goth.distance)

dist.list=make.same.length.vectors(dist.list)
df=do.call(cbind,dist.list)
colnames(df)=c("Nearest.dist","Chicago.dist","Gothic.dist")
colnames(df)=paste(class.list[i],colnames(df))
write.table(df,paste(class.list[i],"distances.txt",sep="_"),sep = "\t",quote = F,col.names = T,row.names = F)
}
