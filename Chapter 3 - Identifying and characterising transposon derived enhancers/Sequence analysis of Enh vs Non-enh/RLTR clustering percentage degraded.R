
setwd('~/R/Thesis Scripts/Results Chapter 1/Sequence analysis of Enh vs Non-enh/')

#clustalo and meme analysis directory
clust.dir="~/R/Thesis Scripts/Data/Clustal_MEME/"


class.list=c("RLTR9E","RLTR9D","RLTR9A3","RLTR13B1","RLTR13B2","RLTR13B3","RLTR13B4","RLTR13D6","RLTR13D5")

#plot percentage histograms
for(n in 1:length(class.list)){

  ltr=class.list[n]


##cutoff for getting only mostly full copies
fa=scan(paste(clust.dir,ltr,'_enh.fa',sep=''), sep='\n', character()) 
id=seq(from = 1,to = length(fa),by = 2)
seq=seq(from = 2,to = length(fa),by = 2)
df.enh=data.frame(seq.length=unlist(lapply(fa[seq], nchar)),id=rep("Enh"))
fa.non=scan(paste(ltr,'_non.fa',sep=''), sep='\n', character()) 
id=seq(from = 1,to = length(fa.non),by = 2)
seq=seq(from = 2,to = length(fa.non),by = 2)
df.non=data.frame(seq.length=unlist(lapply(fa.non[seq], nchar)),id=rep("Non"))
tog.df=rbind(df.enh,df.non)
n.bin=10
max=max(tog.df$seq.length)
breaks=max/n.bin
break.numbers=seq(from=0,to=max,by=breaks)


break.val.enh=c()
break.val.non=c()
for(i in 1:length(break.numbers)-1){
    break.val.enh[i]=(sum(df.enh$seq.length>break.numbers[i]&df.enh$seq.length<=break.numbers[i+1])/nrow(df.enh))*100
    break.val.non[i]=(sum(df.non$seq.length>break.numbers[i]&df.non$seq.length<=break.numbers[i+1])/nrow(df.non))*100
}
df.enh2=data.frame(bar=break.numbers[1:n.bin+1],percentage=break.val.enh,type=rep("Enh"))
head(df.enh2)
df.non2=data.frame(bar=break.numbers[1:n.bin+1],percentage=break.val.non,type=rep("Non"))
tog.df2=rbind(df.enh2,df.non2)
summary(tog.df2)
#write.table(tog.df2,paste(ltr,"_histogram_vals.txt",sep=""),quote = F,col.names = T,row.names = F,sep="\t")

#gg=ggplot(tog.df,aes(x=seq.length,col=as.factor(tog.df$id)))+
#  geom_freqpoly(aes(y = (..count..)/sum(..count..)),bins=20)+
#  ylim(0,0.5)+labs(title=ltr,x="Sequence Length (bp)",y="Percentage")
mround <- function(x,base){
  base*round(x/base)
} 
gg=ggplot(tog.df2,aes(x=bar,y=percentage,fill=type))+geom_col(position = "dodge")+
  labs(title=ltr,x="Sequence Length (bp)",y="Percentage")+scale_fill_manual(values=c("#00ee40", "#C07D44"))+
  scale_x_continuous(breaks=seq(0,mround(max(tog.df2$bar),100),100))
jpeg(paste(ltr,"_seq_length.jpeg",sep=""),w=6,h=3,units='in',res=600)
plot(gg)
dev.off()


}

#get percentage for pie charts
df.row=list()
for(n in 1:length(class.list)){
  ltr=class.list[n]
  fa=scan(paste(ltr,'_enh.fa',sep=''), sep='\n', character()) 
  id=seq(from = 1,to = length(fa),by = 2)
  seq=seq(from = 2,to = length(fa),by = 2)
  max=max(unlist(lapply(fa[seq], nchar)))
  cutoff=max*0.6
  df.enh=data.frame(seq.length=unlist(lapply(fa[seq], nchar)),id=rep("Enh"))
  fa.non=scan(paste(ltr,'_non.fa',sep=''), sep='\n', character()) 
  id=seq(from = 1,to = length(fa.non),by = 2)
  seq=seq(from = 2,to = length(fa.non),by = 2)
  df.non=data.frame(seq.length=unlist(lapply(fa.non[seq], nchar)),id=rep("Non"))
  
  percentage.enh=(sum(df.enh$seq.length>cutoff)/nrow(df.enh))*100
  percentage.non=(sum(df.non$seq.length>cutoff)/nrow(df.non))*100
  df.row[[n]]=data.frame(class=ltr,percentage.enh,percentage.non)
  
}
df=do.call(rbind,df.row)
head(df)
write.table(df,"Percentage_degraded_seq.txt",sep="\t",quote = F,col.names = T,row.names = F)
