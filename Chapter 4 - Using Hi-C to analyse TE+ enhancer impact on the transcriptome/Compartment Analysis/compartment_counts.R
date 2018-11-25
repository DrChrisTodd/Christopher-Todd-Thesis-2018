
setwd("~/R/Bioinformatic for Thesis/active_compartments/")

look=read.delim("Compartment_lookup.txt",h=F)
head(look)

base=read.delim("A_B_compartments.txt",h=F)

base$size=base$V3-base$V2

summary(base$size[base$V4=="A"])
summary(base$size[base$V4=="B"])
sum(base$size[base$V4=="A"])
sum(base$size[base$V4=="B"])

chrm.look=read.delim("../Nearest Gene Analysis/chromosome_lengths.txt")
head(chrm.look)
len=c()
for(i in 1:nrow(chrm.look)){
  sub=as.character(chrm.look$Total.length..bp.)[i]
  x=unlist(strsplit(sub,split = ","))
  y=as.numeric(paste(x,collapse = ""))
  len[i]=y
}


perA=(sum(base$size[base$V4=="A"])/sum(len))*100
perB=(sum(base$size[base$V4=="B"])/sum(len))*100
pernon=(100-sum(perA,perB))

###ROI shuffled for randomised controls generated using bedshuffle tool
#shuffleBed -i ROI_bed.txt -g ./Data/Annotations/mm10 -chrom -seed 927442958 >ROI_shuffle.txt

get.class.from.name=function(x){
  tissue=substr(x,1,3)
  class=unlist(lapply(strsplit(as.character(x),split = "_"),function(x){x[2]}))
  class=paste(tissue,class,sep="_")
  return(class)
}

shuffle=read.delim("ROI_shuffle_compartment_lookup.txt",h=F)

shuffle$ID_class=shuffle$V4
shuffle$Shuff_comp=shuffle$V8
look$ID_class=look$V4
look$Comp=look$V8
colnames(look)
look=look[,9:10]
shuffle=shuffle[,9:10]

ROI=read.delim("../ROI_coord.txt")
head(ROI)
mer=merge(ROI,look,by="ID_class",all.x = T)
mer=merge(mer,shuffle,by = "ID_class",all.x = T)
head(mer)
summary(mer$Comp)

df.row=list()
classes=as.character(unique(mer$class))
for(i in 1:length(classes)){
  sub=mer[mer$class==classes[i],10:11]
  comp.none=sum(is.na(sub$Comp))
  shuff.comp.none=sum(is.na(sub$Shuff_comp))
  comp=as.character(sub$Comp[!is.na(sub$Comp)])
  shuff=as.character(sub$Shuff_comp[!is.na(sub$Shuff_comp)])
  comp.A=sum(comp=="A")
  comp.B=sum(comp=="B")
  shuff.A=sum(shuff=="A")
  shuff.B=sum(shuff=="B")
  df.row[[i]]=data.frame(class=classes[i],num.none=comp.none,num.A=comp.A,num.B=comp.B,num.shuff.none=shuff.comp.none,num.shuff.A=shuff.A,num.shuff.B=shuff.B)
}
df=do.call(rbind,df.row)
write.table(df,"ROI_compartment_counts.txt",sep="\t",quote=F,col.names = T,row.names = F)

####look at what these distributions look like when only consider elements with interactions
##i.e. are only those in active compartments are interacting


hic ='/Users/hmx579/Dropbox/Todd_data/Data/PCHi-C/gothic_interactions_all_xy'
hic.data<- local({
  load(hic)
  stopifnot(length(ls())==1)
  environment()[[ls()]]
})

es<-hic.data$ESC[hic.data$ES.pp==F,]
es.targets = paste(es$chr.target,es$start.target,es$end.target,sep="-")
look=look[,c(4,8)]

got.int=read.delim("../Gothic_Inter_lookup.txt",h=F)
head(got.int)
es.ints=got.int[got.int$V1%in%es.targets,]
head(look)
head(es.ints)
colnames(es.ints)=c("target.ID","ID_class")
colnames(look)=c("ID_class","Comp")
mer=merge(es.ints,look,by="ID_class",all.x=T)
head(mer)
mer$class=get.class.from.name(mer$ID_class)

df.row=list()
classes=as.character(unique(mer$class))
for(i in 1:length(classes)){
  sub=mer[mer$class==classes[i],3:4]
  comp.none=sum(is.na(sub$Comp))
  comp=as.character(sub$Comp[!is.na(sub$Comp)])
  comp.A=sum(comp=="A")
  comp.B=sum(comp=="B")
  df.row[[i]]=data.frame(class=classes[i],num.none=comp.none,num.A=comp.A,num.B=comp.B)
}
df=do.call(rbind,df.row)
head(df)
write.table(df,"ES_ints_ROI_compartment_counts.txt",sep="\t",quote=F,col.names = T,row.names = F)
