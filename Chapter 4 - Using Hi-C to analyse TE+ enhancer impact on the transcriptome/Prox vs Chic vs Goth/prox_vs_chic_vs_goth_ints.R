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



class.list=c("REDE","NEDE")

for(i in 1:length(class.list)){
#subset chicago lookups to class of interest
class=es.chic[grep(es.chic[,2],pattern = class.list[i]),]
es.class=class[grep(class[,2],pattern = "ESC"),]
colnames(bait)=c("bait.id","Gene")
colnames(es.class)=c("int.id","ES_class")
chic.mer=merge(es.class,chic,by="int.id")


chic.bait.mer=merge(chic.mer,bait,by="bait.id")
chic.re.gene=unique(chic.bait.mer[,3:4])
chic.re.gene$Gene=as.character(chic.re.gene$Gene)
chic.re.agg=aggregate(chic.re.gene$Gene,list(chic.re.gene$ES_class),paste)
colnames(chic.re.agg)=c("ES_class_ID","Chic_Genes")


#subset gothic lookups to class of interest
class=es.goth[grep(es.goth[,2],pattern = class.list[i]),]
es.class=class[grep(class[,2],pattern = "ESC"),]
colnames(es.class)=c("int.id","ES_class")


got.mer=merge(es.class,es,by="int.id")
got.bait.mer=merge(got.mer,bait,by="bait.id")
got.re.gene=unique(got.bait.mer[,3:4])
got.re.agg=aggregate(got.re.gene$Gene,list(got.re.gene$ES_class),paste)
colnames(got.re.agg)=c("ES_class_ID","Goth_Genes")


##get nearest df from "Nearest gene analysis" script
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
class.roi=ROI[ROI$class==paste("ESC_",class.list[i],sep=""),]
class.roi$mid=rowMeans(class.roi[,7:8])
nearest.roi=get.nearest.gene(class.roi,1,10)
nearest.roi=nearest.roi[nearest.roi$distance<cutoff,]

nearest=nearest.roi[,c(4,12)]
colnames(nearest)=c("ES_class_ID","Nearest_Genes")


tog.mer=merge(nearest,chic.re.agg,by="ES_class_ID",all=T)
tog.mer=merge(tog.mer,got.re.agg,by="ES_class_ID",all=T)


nearest.genes=unique(as.character(tog.mer$Nearest_Genes))
goth.genes=unique(unlist(tog.mer$Goth_Genes))
goth.genes=goth.genes[!is.na(goth.genes)]
chic.genes=unique(unlist(tog.mer$Chic_Genes))
chic.genes=chic.genes[!is.na(chic.genes)]
c.g=c()
for(i in 1:nrow(tog.mer)){
if(!is.na(tog.mer$Chic_Genes[i])){
  c.g[i]=paste(as.character(unlist(tog.mer$Chic_Genes[i])),collapse = "|")}
  else{c.g[i]=NA}
}

tog.mer$Chic_Genes=c.g
g.g=c()
for(i in 1:nrow(tog.mer)){
  if(!is.na(tog.mer$Goth_Genes[i])){
    g.g[i]=paste(as.character(unlist(tog.mer$Goth_Genes[i])),collapse = "|")}
  else{g.g[i]=NA}
}
tog.mer$Goth_Genes=g.g
write.table(tog.mer,paste(class.list[i],"prox_vs_chic_vs_goth.txt",sep="_"),sep = "\t",col.names = T,row.names = F,quote = F)

#writing files in form where can give to venn diagram website
write.table(nearest.genes,paste(class.list[i],"nearest_genes,txt",sep="_"),sep="\t",col.names = F,row.names = F,quote = F)
write.table(goth.genes,paste(class.list[i],"gothic_genes.txt",sep="_"),sep="\t",col.names = F,row.names = F,quote = F)
write.table(chic.genes,paste(class.list[i],"chicago_genes.txt",sep="_"),sep="\t",col.names = F,row.names = F,quote = F)

}

