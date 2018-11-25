setwd("~/R/Bioinformatic for Thesis/TET_KD/")



##Processed data from Miguel (coordinates are in mm9)
meth=read.delim("./Hon_TE_BSpipeline.txt")
meth$ID=paste("RLTR_",seq(1:nrow(meth)),sep="")


##Writing coord for intersect bed with ROI
#meth.cord=meth[,c(2,3,4,8)]
#meth.cord$Chromosome=paste("chr",meth.cord$Chromosome,sep="")
#write.table(meth.cord,"Hon_TE_BS_coord_mm9.txt",sep="\t",col.names = F,row.names = F,quote = F)



#Intersectbed output file
min=read.delim("Hon_min_in.txt",h=F)


min$RE=unlist(lapply(strsplit(as.character(min$V4),split = "_"),function(x){x[1]}))
min$class=unlist(lapply(strsplit(as.character(min$V4),split = "_"),function(x){x[2]}))
min$hon.id=min$V8

hon.class=c()
for(i in 1:length(unique(min$hon.id))){
  hon.id=as.character(unique(min$hon.id))[i]
  sub=min[min$hon.id==hon.id,]
  if(sum(grepl(sub$class,pattern = "REDE"))>0){hon.class[i]="REDE"}else 
  if(sum(grepl(sub$class,pattern = "NonEnh"))>0){hon.class[i]="NonEnh"}else
  {hon.class[i]="No_Class"}
  }

df=data.frame(ID=unique(min$hon.id),class=hon.class)

mer=merge(meth,df,by="ID")
mer=mer[,c(1,13:19)]


class.df=list()
for (i in 1:length(unique(mer$class))){
  class=as.character(unique(mer$class))[i]
  max=max(summary(as.factor(mer$class)))
  sub=mer[mer$class==class,2:7]
  colnames(sub)=paste(class,colnames(sub),sep="_")
  app.var=rep("",max-nrow(sub))
  df=data.frame(app.var,app.var,app.var,app.var,app.var,app.var)
  colnames(df)=colnames(sub)
  df=rbind(sub,df)
  class.df[[i]]=df
  }

tog.df=do.call(cbind,class.df)
write.table(tog.df,"Hon_BS_graphpad.txt",sep="\t",col.names=T,row.names=F,quote = F)

