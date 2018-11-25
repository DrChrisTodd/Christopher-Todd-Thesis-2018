setwd('~/R/Bioinformatic for Thesis')


directory="~/R/Bioinformatic for Thesis/RLTR_clustering/Cutoffs/"
es.enh=read.delim("ESC_enhancer_TEs.txt",h=F)
ts.enh=read.delim("TSC_enhancer_TEs.txt",h=F)
non.enh=read.delim("NonEnhTE_candidates.txt")
#rltr.list=c("RLTR9","RLTR13D6","RLTR13D5","RLTR13B")
#rltr.list=c("RLTR9E","RLTR9D","RLTR9A3","RLTR13B1","RLTR13B2","RLTR13B3","RLTR13B4","RLTR13D5","RLTR13D6")
rltr.list=c("RLTR9E","RLTR9D","RLTR9A3","RLTR13B2","RLTR13B3","RLTR13B4","RLTR13D5","RLTR13D6")

enh=rbind(es.enh,ts.enh)

#class=rltr.list[2]

library(RColorBrewer)
library(scales)
library(ggplot2)
get.percentage.df=function(class,cutoff){
  
fa = scan(paste(directory,class,"_enh.fa",sep=""),sep='\n', character())
enh.id=gsub(">","",fa[grep(">",fa)])
fa = scan(paste(directory,class,"_non.fa",sep=""),sep='\n', character())
non.id=gsub(">","",fa[grep(">",fa)])

fimo.enh=read.delim(paste("./RLTR_clustering/Cutoffs/",class,"_enh_fimo.txt",sep=""))
fimo.non=read.delim(paste("./RLTR_clustering/Cutoffs/",class,"_non_fimo.txt",sep=""))

tog.fimo=rbind(fimo.enh,fimo.non)
fimo.ids=unique(as.character(append(as.character(fimo.enh$motif_alt_id),as.character(fimo.non$motif_alt_id))))

for(i in 1:length(fimo.ids)){
  id=fimo.ids[i]
  fimo.sub=tog.fimo[tog.fimo$motif_alt_id==id,]  
  id.match=(enh.id %in% fimo.sub$sequence_name)*1
  if(i==1){df=data.frame(id=enh.id,TF=id.match)}else{df=cbind(df,id.match)}
}
colnames(df)=c("ID",fimo.ids)
enh.df=df
for(i in 1:length(fimo.ids)){
  id=fimo.ids[i]
  fimo.sub=tog.fimo[tog.fimo$motif_alt_id==id,]  
  id.match=(non.id %in% fimo.sub$sequence_name)*1
  if(i==1){df2=data.frame(id=non.id,TF=id.match)}else{df2=cbind(df2,id.match)}
}
colnames(df2)=c("ID",fimo.ids)
non.df=df2

sum.enh=colSums(enh.df[,2:ncol(enh.df)])
sum.non=colSums(non.df[,2:ncol(non.df)])

per.enh=(sum.enh/length(enh.id))*100
per.non=(sum.non/length(non.id))*100
per.diff=per.enh-per.non

df=data.frame(per.enh,per.non,per.diff,fimo.ids)


ame=read.delim(paste("./RLTR_clustering/Cutoffs/",class,"_ame.txt",sep=""))
if(nrow(ame)>7){
ame=as.character(ame[8:nrow(ame),1])
ame.id=unlist(lapply(strsplit(ame,split = " |)"),function(x){x[9]}))
ame.pval=as.numeric(unlist(lapply(strsplit(ame,split = " |)"),function(x){x[length(x)]})))
ame.df=data.frame(id=ame.id,pval=ame.pval)}else{
  ame.df=data.frame(id="null",pval="1")
}
pval=c()
for(n in 1:nrow(df)){
  id=as.character(df[n,4])
  ame.sub=ame.df[ame.df$id==id,]  
  if(nrow(ame.sub)==0){pval[n]=1}else{pval[n]=ame.sub[1,2]}
  }
df$pval=pval


return(df)}

for(i in 1:length(rltr.list)){
class=rltr.list[i]

df=get.percentage.df(class)
df$FC=(df$per.enh+0.01)/(df$per.non+0.01)
label=c()
for(n in 1:nrow(df)){
  sub=df[n,]
if(sub$FC>2&sub$per.enh>40&sub$pval<0.05){label[n]=rownames(sub)}else{label[n]=""}
}

df$label=label


jpeg(paste(class,"_motif_percentages.jpeg",sep=""),width = 6,height = 6,units='in',res=600)
gg=ggplot(df,aes(x=per.enh,y=per.non,col=as.factor(df$FC>2&df$per.enh>40&df$pval<0.05)))+geom_point()+ggtitle(class)+
  geom_abline(slope = 1,intercept = c(0,0),linetype=2,colour="red")+ylim(0,100)+xlim(0,100)+
  theme_classic()+theme(legend.position="none",axis.text = element_text(face="bold", size=18),axis.title = element_text(face="bold", size=18), title =element_text(face="bold", size=24) )+
  geom_text(label=df$label, nudge_x = 12, nudge_y=2,color="black", size=5,fontface = "bold")+labs(x="Percentage of Enhancer Copies",y="Percentage of Non-Enhancer Copies")
plot(gg)
dev.off()
}
