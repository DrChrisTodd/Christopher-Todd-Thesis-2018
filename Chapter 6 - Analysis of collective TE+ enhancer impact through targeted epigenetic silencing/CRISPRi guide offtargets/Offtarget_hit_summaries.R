setwd("~/CRISPRi_guides")

##path to bedtools
bedtools = '~/bedtools2/bin/'


##import functions
source('~/todd/Functions/intersectBed.R')

Repmask=read.delim("./RepeatMasker.txt",h=F)



files=list.files(pattern = "Cas_OFF")
for(n in 1:length(files)){
file=read.delim(files[n])
guide=unlist(strsplit(files[n],split = "_"))[1]
start=c()
for(i in 1:nrow(file)){
  if(file[i,6]=="+"){start[i]=file[i,5]}else
    if(file[i,6]=="-"){start[i]=file[i,5]-20}
}
file$start=start
file$end=file$start+20
order=order(file$Mismatches,decreasing = F)
file=file[order,]


name=paste(guide,"_",seq(1:nrow(file)),"_MM",file$Mismatches,sep="")
bed=data.frame(file$Chromosome,file$start,file$end,name,(5-file$Mismatches),file$Direction)
write.table(bed,paste(guide,"_casoff.bed",sep=""),sep="\t",col.names=F,row.names=F,quote=F)

rep.int = intersectBed(bed,Repmask,opt.string='-wao ',path.to.bedtools=bedtools)


rep.int$V10=as.factor(rep.int$V10)
x=list()
for(i in unique(rep.int$V5)){
  sub=rep.int[rep.int$V5==i,]
  x[[i]]=summary(sub$V10)
}
df=do.call(cbind,x)
colnames(df)=c("MM4","MM3","MM2","MM1","MM0")

write.table(df,paste(guide,"repeat_mismatch_hits.txt",sep="_"),sep="\t",col.names = T,row.names = T,quote=F)
}