library(circlize)
library(ComplexHeatmap)

rev.comp<-function(dna){
  #reverse compliment function
  seq<-strsplit(dna,split="")[[1]]
  seq.rev<-rev(seq)
  seq.rev<-paste(seq.rev,collapse = "")
  sub<-gsub("C","g",seq.rev)
  sub<-gsub("G","c",sub)
  sub<-gsub("A","t",sub)
  sub<-gsub("T","a",sub)
  revcom<-toupper(sub)
  return(revcom)}
group.guide.lookup<-function(seq,guides){
  for(i in 1:length(guides)){
    guide=guides[i]
    rev.guide=rev.comp(guide)
    guide.options<-c(paste(guide,"AGG",sep=""),paste(guide,"TGG",sep=""),paste(guide,"CGG",sep=""),paste(guide,"GGG",sep=""))
    rev.guide.options<-c(paste("CCA",rev.guide,sep=""),paste("CCT",rev.guide,sep=""),paste("CCC",rev.guide,sep=""),paste("CCG",rev.guide,sep=""))
    guide.hits=grepl(x=seq,pattern = paste(guide.options,rev.guide.options,sep="|",collapse = "|"))
    if(i == 1){df=data.frame(names(seq),guide.hits)}else{df[,i+1]=guide.hits}
    colnames(df)[i+1]=paste("Guide_",i,sep="")
  }
  if(ncol(df)>2){guide.hits=rowSums(df[,2:ncol(df)])}else{guide.hits=df[,2]*1}
  df=data.frame(df[,1],guide.hits)
  colnames(df)=c("names.seq","guide.hits")
  return(df)}


fasta = scan("RLTR13D6.fa", sep='\n', character()) 
newseq = grep('>',fasta)

seq = character(length(newseq))
for (i in 1:length(newseq)) {
  first = newseq[i]+1
  if (i==length(newseq)) {
    last = length(fasta)
  } else {
    last = newseq[i+1]-1
  }
  seq[i] = toupper(paste(fasta[first:last],collapse=''))
}
names(seq) = unlist(lapply(strsplit(fasta[newseq],'[> ]'),function(x) x[2]))



for(n in 1:8){

guide.name=paste("guide",n,sep="")

x=read.delim(paste("./",guide.name,"_Cas_OFFinder.txt",sep=""))

mm0=x[x$Mismatches==0,]
mm0seq=substr(as.character(mm0[1,2]),1,20)

mm1=x[x$Mismatches==1,]
mm1seqs=toupper(unique(substr(as.character(mm1[,3]),1,20)))

mm2=x[x$Mismatches==2,]
mm2seqs=toupper(unique(substr(as.character(mm2[,3]),1,20)))

mm0.hits=group.guide.lookup(seq,guides = mm0seq)
mm1.hits=group.guide.lookup(seq,guides = mm1seqs)
mm2.hits=group.guide.lookup(seq,guides = mm2seqs)

df=data.frame(mm0.hits,mm1.hits[,2],mm2.hits[,2])
colnames(df)=c("names.seq","mm0","mm1","mm2")

score=rep(0,nrow(df))
for(i in 1:nrow(df)){
  if(df[i,4]==1){score[i]=1}else
    if(df[i,3]==1){score[i]=2}else
      if(df[i,2]==1){score[i]=3}
}
df2=data.frame(df[,1],score)
colnames(df2)=c("seq.names",paste(guide.name,"score",sep="_"))
if(n==1){df3=df2}else{df3=cbind(df3,df2[,2])}
}
colnames(df3)=c("seq.names","guide1_score","guide2_score","guide3_score","guide4_score","guide5_score","guide6_score","guide7_score","guide8_score")
head(df3)
Heatmap(as.matrix(df3[,2:5]),cluster_columns = F)
df3[df3==1]=0
order=order(rowSums(df3[,2:5]),decreasing=T)
df3.group1=df3[order,1:5]



group1.rede=df3.group1[grep(df3.group1$seq.names,pattern = "REDE"),]
Heatmap(as.matrix(group1.rede[,2:5]),cluster_columns = F,cluster_rows = T,col = colours)
group1.non=df3.group1[grep(df3.group1$seq.names,pattern = "noclass"),]
Heatmap(as.matrix(group1.non[,2:5]),cluster_columns = F,cluster_rows = T,col = colours)

colours=structure(c("#cccbcb","#cccbcb","#fb8888","#fd0000"),names=c(0,1,2,3))
Heatmap(as.matrix(group1.rede[,2:5]),cluster_columns = F,cluster_rows = F,col = colours)


order=order(rowSums(df3[,6:9]),decreasing=T)
df3.group2=df3[order,c(1,6:9)]
group2.rede=df3.group2[grep(df3.group2$seq.names,pattern = "REDE"),]
Heatmap(as.matrix(group2.rede[,2:5]),cluster_columns = F,cluster_rows = T,col = colours)
group2.non=df3.group2[grep(df3.group2$seq.names,pattern = "noclass"),]
Heatmap(as.matrix(group2.non[,2:5]),cluster_columns = F,cluster_rows = T,col = colours)
