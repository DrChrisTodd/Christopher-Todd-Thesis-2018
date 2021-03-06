###Getting candidate REs which are non Enh TEs
setwd('~/todd')
##will use the ATAC data and mappability score
##nonEnh elements will be those with little ATAC but at least a certain level of mappability

TE.look=read.delim("Data/Output/TE_coordinates_lookup.txt")

##perform homer ATAC heatmap analysis for TE_coord dataset and then return values so can pick lowest scoring


##Set path to homer
homer = '~/Homer/bin/'
tag.dir = '~/chris/Homer_Tags/'
es.mark = 'ES_ATAC'
ts.mark = 'TS_ATAC'

##import functions
source('/data/Blizard-BrancoLab/Chris/Todd_et_al/Functions/annotatePeaks.R')

es.REs=c("RLTR13D6","RLTR9")
ts.REs=c("RLTR13D5","RLTR13B")
map.lim = 0.5

#define non-enh to enh ratio
ratio=2

###Annotate Peaks function from Homer Package

#putting TE coords in bed format
bed=TE.look[,c(1,9:10,8,7,4)]

##home annotatepeaks function from homer
temp.bed.file<-tempfile()
write.table(bed,file=temp.bed.file,quote=F,sep='\t',col.names=F,row.names=F)
es.hist=annotatepeaks(bed.file = temp.bed.file,tag=paste(tag.dir,es.mark,sep=""),ref.genome = ' mm10 ',path.to.homer = homer)
ts.hist=annotatepeaks(bed.file = temp.bed.file,tag=paste(tag.dir,ts.mark,sep=""),ref.genome = ' mm10 ',path.to.homer = homer)
unlink(temp.bed.file)

#summarising each TE into single score
es.hist=data.frame(ID_class=es.hist[,1],score=rowSums(es.hist[,2:ncol(es.hist)]))
ts.hist=data.frame(ID_class=ts.hist[,1],score=rowSums(ts.hist[,2:ncol(ts.hist)]))



#getting ES relevant nonEnh TEs
for(i in 1:length(es.REs)){
  rep.sub=TE.look[grep(TE.look$repName,pattern=es.REs[i]),]
  enh.sub=rep.sub[grep(rep.sub$ID_class,pattern="REDE"),]
  non.REDE=rep.sub[!(grepl(rep.sub$ID_class,pattern="REDE")),]
  map.non=non.REDE[non.REDE$mappability>map.lim,]
  
  mer=merge(map.non,es.hist,by="ID_class",all.x = T)
  mer[is.na(mer)]=0
  mer.ord=mer[order(mer$score,decreasing = F),]
  nonenhgroup=mer.ord[1:(nrow(enh.sub)*ratio),]
  if(i==1){es.df=nonenhgroup}else{es.df=rbind(es.df,nonenhgroup)}
  }
#getting TS relevant nonEnh TEs
for(i in 1:length(ts.REs)){
  rep.sub=TE.look[grep(TE.look$repName,pattern=ts.REs[i]),]
  enh.sub=rep.sub[grep(rep.sub$ID_class,pattern="REDE"),]
  non.REDE=rep.sub[!(grepl(rep.sub$ID_class,pattern="REDE")),]
  map.non=non.REDE[non.REDE$mappability>map.lim,]
  
  mer=merge(map.non,ts.hist,by="ID_class",all.x=T)
  mer[is.na(mer)]=0
  mer.ord=mer[order(mer$score,decreasing = F),]
  nonenhgroup=mer.ord[1:(nrow(enh.sub)*ratio),]
  if(i==1){ts.df=nonenhgroup}else{ts.df=rbind(ts.df,nonenhgroup)}
}

tog.df=rbind(es.df,ts.df)
tog.df=tog.df[,c(2:8,1,9:10)]
tog.df$ID_class=paste(unlist(lapply(strsplit(as.character(tog.df$ID_class),split = "_"),function(x){x[1]})),"NonEnh",sep="_")

write.table(tog.df,"Data/Output/NonEnhTE_candidates.txt",sep="\t",quote = F,col.names = T,row.names = F)
