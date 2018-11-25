#### linking RE elements with HiC data

##Set work directory
setwd('~/todd')

##regions of interest that want to look up in the HiC datasets
ROI= read.delim("Data/Output/ROI_coord.txt")
##Defining each class to be treated differently
classes<-unique(as.character(ROI[,9]))

##HiC datasets
hic='Data/PCHiC/gothic_interactions_all_xy'

##import functions
source('Functions/intersectBed.R')
##path to bedtools
bedtools = '~/bedtools2/bin/'


#####Script######

hic.data<- local({
  load(hic)
  stopifnot(length(ls())==1)
  environment()[[ls()]]
}) 
es<-hic.data$ESC[hic.data$ES.pp==F,]
es$target.ID<-paste(es$chr.target,es$start.target,es$end.target,sep="-")
es$bait.ID<-paste(es$chr.bait,es$start.bait,es$end.bait,sep="-")
es$int.ID<-paste(es$bait.ID,es$target.ID,sep=":")
ts<-hic.data$TSC[hic.data$TS.pp==F,]
ts$target.ID<-paste(ts$chr.target,ts$start.target,ts$end.target,sep="-")
ts$bait.ID<-paste(ts$chr.bait,ts$start.bait,ts$end.bait,sep="-")
ts$int.ID<-paste(ts$bait.ID,ts$target.ID,sep=":")


##Generating a lookup for each interaction region in hic dataset
#join region of interest bed coordinates together
roi<-ROI[,c(1,7:8,4)]
hic.target<-unique(rbind(es[,c(7:9,11)],ts[,c(7:9,11)]))
target.look = intersectBed(hic.target,roi,opt.string='-wa -wb -F 0.5',path.to.bedtools=bedtools)
target.look<-unique(target.look[,c(4,8)])
write.table(target.look,"Data/Output/Gothic_Inter_lookup.txt",sep="\t",quote = F,row.names = F,col.names = F)


##getting aggregate values for each inter region for each of the 6 groups
##es TE enh, ts TE enh, es TE non enh, ts TE non enh, ES nonTE enh, ts nonTE enh

for(i in 1:length(classes)){
  class.id<-as.character(ROI[ROI[,9]==classes[i],4])
  class.log<-target.look[,2] %in% class.id
  target.look$class<-class.log*1
  agg<-aggregate(target.look$class,list(target.look[,1]),sum)
  colnames(agg)=c("target.ID",classes[i])
  if(i==1){int.mer=agg}else{int.mer=merge(int.mer,agg,by="target.ID",all=T)}
}

write.table(int.mer,"Data/Output/Gothic_Inter_Agg_lookup.txt",sep="\t",quote = F, col.names = T,row.names = F)

