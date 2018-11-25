####Going to adapt my script so it will select ALL possible CRISPR guide options and then process these guides instead of using the clustalo analysis
##This will require getting the fasta sequences of each member of RE family
##For each sequence in the fasta generating a list of CRISPR guides with NGG pam 
##condensing this list by removing all replicates
##proceed as usual after that


setwd('~/chris/Ozgen_guides')
##Bed file of RE coordinates ##make sure its mm10 coordinates
##list format
bed.list = c("~/chris/Ozgen_guides/LTR5Hs.txt","~/chris/Ozgen_guides/LTR5B.txt")
RE.name.list = c("LTR5Hs","LTR5B")
homer.path = '~/Homer'

for(i in 1:length(bed.list)){

bed.input=bed.list[i]
RE.name=RE.name.list[i]
##########
command = paste(homer.path,'/bin/homerTools extract ',bed.input," ",homer.path,'/data/genomes/mm10/ -fa > ',RE.name,'.fa',sep='')
try(system(command))

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

fasta = scan(paste(RE.name,".fa",sep=""), sep='\n', character()) 
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
head(seq)

for.guide<-list()
for(z in 1:length(seq)){
fa.seq<-seq[z]
x<-gregexpr('GG',fa.seq)
len<-length(x[[1]])
guide<-list()
for(i in 1:len){
pam<-x[[1]][i]
log<-pam<23
if(log==T) next
guide[[i]]<-substr(fa.seq,pam-21,pam-2)}
for.guide[[z]]<-unlist(guide)}
rev.guide<-list()
for(z in 1:length(seq)){
  z=1
  fa.seq<-rev.comp(seq[z])
  x<-gregexpr('GG',fa.seq)
  len<-length(x[[1]])
  guide<-list()
  for(i in 1:len){
    pam<-x[[1]][i]
    log<-pam<23
    if(log==T) next
    guide[[i]]<-substr(fa.seq,pam-21,pam-2)}
  rev.guide[[z]]<-unlist(guide)
}
tog.guide<-append(for.guide,rev.guide)

unique.guide<-unique(unname(unlist(tog.guide)))
unique.guide<-unique.guide[nchar(unique.guide)==20]
guide.df<-data.frame(seq(1:length(unique.guide)),unique.guide)
guide.df[,1]<-paste(RE.name,"guide",guide.df[,1],sep="_")
colnames(guide.df)<-c("Guide_ID","Guide_seq")
guide.df$Guide_RevCom_seq<-unlist(lapply(as.character(guide.df$Guide_seq),rev.comp))
head(guide.df)
write.table(paste(RE.name,"_guide_options.txt",sep=""),sep="\t",row.names = F,col.names = T,quote = F)


guide.lookup<-function(seq){
  x<-list()
  for(n in 1:nrow(guide.df)){
    guide.options<-c(paste(substr(guide.df[n,2],1,20),"AGG",sep=""),paste(substr(guide.df[n,2],1,20),"TGG",sep=""),paste(substr(guide.df[n,2],1,20),"CGG",sep=""),paste(substr(guide.df[n,2],1,20),"GGG",sep=""))
    rev.guide.options<-c(paste("CCA",substr(guide.df[n,3],4,23),sep=""),paste("CCT",substr(guide.df[n,3],4,23),sep=""),paste("CCC",substr(guide.df[n,3],4,23),sep=""),paste("CCG",substr(guide.df[n,3],4,23),sep=""))
    x[[n]]<-((grepl(x = seq, pattern = paste(guide.options,rev.guide.options,sep="|",collapse = "|")))*1)
  }
  y<-unlist(x)
  names(y)<-guide.df$guide.name
  return(y)}
mat<-lapply(FUN = guide.lookup,X=seq)
mat2 = matrix(unlist(mat),ncol=length(mat))
row.names(mat2)<-guide.df$Guide_ID
colnames(mat2)<-names(seq)  

head(mat2)
ord<-order(rowSums(mat2),decreasing = T)
mat3<-mat2[ord,]
head(mat3)

library(gplots)
pdf(paste(RE.name,"guides_ordered_heatmap.pdf",sep="_"),width=8,height=8)
heatmap.2(mat3,col=c('black','red'),key=FALSE,trace='none',scale='none',Rowv=F,Colv=F, dendrogram = "none",main=RE.name,srtCol=45,srtRow = 45,cexRow = 0.2,cexCol = 0.2)
dev.off()

write.table(mat3,paste(RE.name,"guides_table.txt",sep="_"),col.names = T,row.names = T,sep="\t",quote = F)


x<-mat3
###Going to attempt a method of picking guides by removing the best guide each time (putting it into a new df) and the reps it targets
###Then will select the next best guide ect to try pick a selection of guides which will cover the highest number of reps

#order by the highest number of REs hit
step<-x[order(rowSums(x),decreasing = T),]

#use top row as best one
guide1<-x[row.names(step)[1],]
#remove this guide and the associated REs
REs.hit<-step[1,]==1
step<-step[-c(1),REs.hit==F]
step<-step[order(rowSums(step),decreasing=T),]
guide2<-x[row.names(step)[1],]

REs.hit<-step[1,]==1
step<-step[-c(1),REs.hit==F]
step<-step[order(rowSums(step),decreasing=T),]
guide3<-x[row.names(step)[1],]

REs.hit<-step[1,]==1
step<-step[-c(1),REs.hit==F]
step<-step[order(rowSums(step),decreasing=T),]
guide4<-x[row.names(step)[1],]


guides<-rbind(guide1,guide2,guide3,guide4)
mat<-as.matrix(guides)
colnames(mat)<-colnames(guides)
rownames(mat)<-rownames(guides)


pdf(paste(RE.name,"_most_hit.pdf",sep=""),width=8,height=8)
heatmap.2(mat,col=c('black','red'),key=FALSE,trace='none',scale='none', dendrogram = "none",main=paste(RE.name,"_most_hit",sep=""),srtCol=45,srtRow = 45,cexRow = 0.4,cexCol = 0.2)
dev.off()

###second approach is to simply take the top 4-6 guides and plot onto heatmap to see how well they do

guides<-x[1:8,]
mat<-as.matrix(guides)
colnames(mat)<-colnames(guides)
rownames(mat)<-rownames(guides)

pdf(paste(RE.name,"_top_guides.pdf",sep=""),width=8,height=8)
heatmap.2(mat,col=c('black','red'),key=FALSE,trace='none',scale='none', dendrogram = "none",main=paste(RE.name,"_top_guides",sep=""),srtCol=45,srtRow = 45,cexRow = 0.4,cexCol = 0.2)
dev.off()
}

