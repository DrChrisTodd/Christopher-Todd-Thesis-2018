#### linking HiC data with Gene expression

##Set work directory
#setwd('~/todd')
setwd('~/R/Thesis Scripts/Results Chapter 2/Gothic Anlysis/')

##Set R object containing HiC data
#hic = './Data/PCHiC/gothic_interactions_all_xy'
hic ='/Users/hmx579/Dropbox/Todd_data/Data/PCHi-C/gothic_interactions_all_xy'

bait.file="~/R/Thesis Scripts/Data/PCHi-C/bait_RNAseq_look.txt"
gene.file<-"/Users/hmx579/Dropbox/Todd_data/Data/RNA-seq/TSC_RNA_FPKM.txt.gz"




hic.data<- local({
  load(hic)
  stopifnot(length(ls())==1)
  environment()[[ls()]]
})

es<-hic.data$ESC[hic.data$ES.pp==F,]
es$name.bait<-unlist(lapply(strsplit(es$name.bait,split = "-|,"),function(x){x[[1]]}))
es$target.ID = paste(es$chr.target,es$start.target,es$end.target,sep="-")
es$unique<-hic.data$ES.specific[hic.data$ES.pp==F]
es$bait.ID= paste(es$chr.bait,es$start.bait,es$end.bait,sep="-")


ts<-hic.data$TSC[hic.data$TS.pp==F,]
ts$name.bait<-unlist(lapply(strsplit(ts$name.bait,split = "-|,"),function(x){x[[1]]}))
ts$target.ID = paste(ts$chr.target,ts$start.target,ts$end.target,sep="-")
ts$unique<-hic.data$TS.specific[hic.data$TS.pp==F]
ts$bait.ID= paste(ts$chr.bait,ts$start.bait,ts$end.bait,sep="-")

#agg.look = "./Data/Output/Gothic_Inter_Agg_lookup.txt"
agg.look='Gothic_Inter_Agg_lookup.txt'
agg = read.delim(agg.look)


for(n in 1:2){
  if(n==1){type=es}else{type=ts}
type.agg<-merge(type,agg,by="target.ID",all.x=T)

type.agg<-type.agg[,c(13:21)]
type.agg[is.na(type.agg)]<-0
get.agg=function(df,ID.name){
for(i in 1:(ncol(df)-1)){
  class.agg<-aggregate(df[,(1+i)],list(df[,1]),sum)
  colnames(class.agg)<-c(ID.name,colnames(df)[(1+i)])
  if(i==1){id.agg=class.agg}else{id.agg<-merge(id.agg,class.agg,by=ID.name)}
}
return(id.agg)}
bait.agg=get.agg(type.agg,"Bait_ID")


if(n==1){write.table(bait.agg,"ES_gothic_bait_agg_look.txt",sep="\t",quote = F,col.names = T,row.names = F)}
if(n==2){write.table(bait.agg,"TS_gothic_bait_agg_look.txt",sep="\t",quote = F,col.names = T,row.names = F)}

gene.id<-read.delim(gene.file)[,c(1,6:9)]
colnames(gene.id)[1]<-"Gene_ID"
bait.look=read.delim(bait.file,as.is = T,stringsAsFactors = F,h=F)
colnames(bait.look)=c("Bait_ID","Gene_ID")
gene.mer=merge(gene.id,bait.look,by="Gene_ID")


bait.mer<-merge(gene.mer,bait.agg,by="Bait_ID",all.x=T)
bait.mer[is.na(bait.mer)]=0
gene.info=unique(bait.mer[,c(2:6)])
int.info=bait.mer[,c(2,7:ncol(bait.mer))]
gene.agg=get.agg(int.info,"Gene_ID")
gene.mer=merge(gene.agg,gene.info,by="Gene_ID")

##now I've linked baits to genes perform second aggregate function 

gene.mer$ES_mean<-(gene.mer$ES_E14+gene.mer$ES_J1)/2
gene.mer$TS_mean<-(gene.mer$TS_EGFP+gene.mer$TS_Rs26)/2
gene.mer$Rel_exp<-(gene.mer$ES_mean-gene.mer$TS_mean)

if(n==1){
  gene.mer$ES_REDE_only=gene.mer$ESC_REDE>0&gene.mer$ESC_NEDE==0
  gene.mer$ES_NEDE_only=gene.mer$ESC_REDE==0&gene.mer$ESC_NEDE>0
  gene.mer$ES_nonEnhTE_only=gene.mer$ESC_REDE==0&gene.mer$ESC_NEDE==0&gene.mer$ESC_NonEnh>0
  gene.mer$No_ES_Enh_only=gene.mer$ESC_REDE==0&gene.mer$ESC_NEDE==0
  gene.mer$All=rep(T)
}
if(n==2){
  gene.mer$TS_REDE_only=gene.mer$TSC_REDE>0&gene.mer$TSC_NEDE==0
  gene.mer$TS_NEDE_only=gene.mer$TSC_REDE==0&gene.mer$TSC_NEDE>0
  gene.mer$TS_nonEnhTE_only=gene.mer$TSC_REDE==0&gene.mer$TSC_NEDE==0&gene.mer$TSC_NonEnh>0
  gene.mer$No_TS_Enh_only=gene.mer$TSC_REDE==0&gene.mer$TSC_NEDE==0
  gene.mer$All=rep(T)
}
if(n==1){write.table(gene.mer,"ES_gothic_gene_lookup.txt",sep="\t",quote = F,col.names = T,row.names = F)}
if(n==2){write.table(gene.mer,"TS_gothic_gene_lookup.txt",sep="\t",quote = F,col.names = T,row.names = F)}

for(i in 1:5){
  len=max(colSums(gene.mer[,17:21]))
  sub=gene.mer[gene.mer[,(16+i)]==T,]
  es_mean<-sub$ES_mean
  ts_mean<-sub$TS_mean
  min.exp=es_mean>(-1)|ts_mean>(-1)
  Rel_exp<-sub$Rel_exp
  Rel_exp[!min.exp]=""
  es_mean<-append(es_mean,rep(" ",(len-(nrow(sub)))))
  ts_mean<-append(ts_mean,rep(" ",(len-(nrow(sub)))))
  Rel_exp<-append(Rel_exp,rep(" ",(len-(nrow(sub)))))
  class.name<-substr(colnames(gene.mer)[(16+i)],1,nchar(colnames(gene.mer)[(16+i)])-4)
  df<-data.frame(es_mean,ts_mean,Rel_exp)
  colnames(df)<-c(paste(class.name,"es_mean",sep=""),paste(class.name,"ts_mean",sep=""),paste(class.name,"Rel_exp",sep=""))
  if(i==1){graph.pad.df<-df}else{graph.pad.df<-cbind(graph.pad.df,df)}
}
if(n==1){write.table(graph.pad.df,"ES_gothic_gene_exp.txt",sep="\t",quote = F,col.names = T,row.names = F)}
if(n==2){write.table(graph.pad.df,"TS_gothic_gene_exp.txt",sep="\t",quote = F,col.names = T,row.names = F)}


}

