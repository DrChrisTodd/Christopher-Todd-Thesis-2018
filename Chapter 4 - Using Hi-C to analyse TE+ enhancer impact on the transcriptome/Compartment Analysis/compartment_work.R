setwd("~/R/Thesis Scripts/Results Chapter 2/Compartment Analysis/")

bait.comp=read.delim("bait_compartment_lookup.txt",h=F)[,c(4,8)]
bait.look=read.delim("~/R/Thesis Scripts/Data/PCHi-C/bait_RNAseq_look.txt")


colnames(bait.comp)=c("Comp","bait.ID")
colnames(bait.look)=c("bait.ID","Gene_ID")
mer=merge(bait.look,bait.comp,by="bait.ID",all.x=T)
mer$Comp=as.character(mer$Comp)
mer[is.na(mer$Comp),3]=paste("N")
mer=mer[,2:3]

gene.look=read.delim("~/R/Thesis Scripts/Results Chapter 2/Gothic Analysis/ES_gothic_gene_lookup.txt")
gene.mer=merge(gene.look,mer,by="Gene_ID")

int.class=c()
for(i in 1:nrow(gene.mer)){
  if(gene.mer[i,17]){int.class[i]="ES_REDE_only"}else
  if(gene.mer[i,18]){int.class[i]="ES_NEDE_only"}else
  if(gene.mer[i,19]){int.class[i]="ES_nonEnhTE_only"}else
  if(gene.mer[i,19]==F&gene.mer[i,20]==T){int.class[i]="No_ES_Enh_only"}else{int.class[i]="Mixed"}
}

gene.mer$class=int.class


##write to graphpad format
#start with gene mer
head(gene.mer)
gene.mer$group=paste(gene.mer$Comp,gene.mer$class,sep="_")
groups=as.character(unique(gene.mer$group))
All.A=gene.mer$Comp=="A"
All.B=gene.mer$Comp=="B"
All.N=gene.mer$Comp=="N"

df.col=list()
max=max(sum(All.A),sum(All.B),sum(All.N))
for(i in 1:length(groups)){
  rows=gene.mer$group==groups[i]
  es=gene.mer$ES_mean[rows]
  rel=gene.mer$Rel_exp[rows]
  min.exp=(gene.mer$ES_mean>(-1)|gene.mer$TS_mean>(-1))[rows]
  rel[!min.exp]=""
  df=data.frame(es,rel)
  colnames(df)=paste(c("ES_mean_","Rel_exp_"),groups[i],sep="")
  app.val=rep("",max-sum(rows))
  df2=data.frame(app.val,app.val)
  colnames(df2)=colnames(df)
  df.col[[i]]=rbind(df,df2)
}

df.col2=list()
all.list=list(All.A,All.B,All.N)
all.names=c("All_A","All_B","All_N")
for(i in 1:length(all.names)){
  rows=all.list[[i]]
  es=gene.mer$ES_mean[rows]
  rel=gene.mer$Rel_exp[rows]
  min.exp=(gene.mer$ES_mean>(-1)|gene.mer$TS_mean>(-1))[rows]
  rel[!min.exp]=""
  df=data.frame(es,rel)
  colnames(df)=paste(c("ES_mean_","Rel_exp_"),all.names[i],sep="")
  app.val=rep("",max-sum(rows))
  df2=data.frame(app.val,app.val)
  colnames(df2)=colnames(df)
  df.col2[[i]]=rbind(df,df2)
  
}
df=do.call(cbind,c(df.col,df.col2))
head(df)
summary(df.col2[[2]])
write.table(df,"Compartments_ES_gene_exp_graphpad.txt",sep="\t",col.names = T,row.names = F,quote = F)
