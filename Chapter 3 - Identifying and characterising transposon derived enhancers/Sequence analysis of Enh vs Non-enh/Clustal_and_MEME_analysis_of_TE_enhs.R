setwd('~/MEME_analysis/RLTR_full')

es.enh=read.delim("../ESC_enhancer_TEs.txt",h=F)
ts.enh=read.delim("../TSC_enhancer_TEs.txt",h=F)
non.enh=read.delim("../NonEnhTE_candidates.txt")
rltr.list=c("RLTR9E","RLTR9D","RLTR9A3","RLTR13B1","RLTR13B2","RLTR13B3","RLTR13B4","RLTR13D6","RLTR13D5")

enh=rbind(es.enh,ts.enh)
clustalo = function(ltr) {
  
  ##Get fasta for enh and nonenh groups
  enh.sub=enh[grep(enh[,5],pattern = ltr),c(1:3,8,7,4)]
  non.sub=non.enh[grep(non.enh$repName,pattern = ltr),c(1:3,8,7,4)]
  
  temp.enh.file<-tempfile()
  write.table(enh.sub,file=temp.enh.file,quote=F,sep='\t',col.names=F,row.names=F)
  command = paste('~/Homer/bin/homerTools extract ',temp.enh.file,' ~/Homer/data/genomes/mm10/ -fa > ',ltr,"_enh",".fa",sep='')
  try(system(command))
  unlink(temp.enh.file)
  temp.non.file=tempfile()
  write.table(non.sub,file=temp.non.file,quote=F,sep='\t',col.names=F,row.names=F)
  command = paste('~/Homer/bin/homerTools extract ',temp.non.file,' ~/Homer/data/genomes/mm10/ -fa > ',ltr,"_non",".fa",sep='')
  try(system(command))
  unlink(temp.non.file)
  
  ##cutoff for getting only mostly full copies
  fa=scan(paste(ltr,'_enh.fa',sep=''), sep='\n', character()) 
  id=seq(from = 1,to = length(fa),by = 2)
  seq=seq(from = 2,to = length(fa),by = 2)
  max=max(unlist(lapply(fa[seq], nchar)))
  cutoff=max*0.6
  get.fa.file=function(seq){
    log=nchar(fa[seq])>cutoff
    file.id=list()
    for(i in 1:length(log)){
      if(log[i]==T){
        file.id[[i]]=c(fa[id[i]],fa[seq[i]])
      }else{next}
    }
    new.file=unlist(file.id)
    return(new.file)}
  write.table(get.fa.file(seq), paste(ltr,"_enh.fa",sep=''),sep = "\n",quote = F,col.names = F,row.names = F)
  fa=scan(paste(ltr,'_non.fa',sep=''), sep='\n', character()) 
  id=seq(from = 1,to = length(fa),by = 2)
  seq=seq(from = 2,to = length(fa),by = 2)
  write.table(get.fa.file(seq),paste(ltr,"_non.fa",sep=""),sep="\n",quote = F,col.names = F,row.names = F)
  
  
  
  
  
  
  
  ltr.subs=c(paste(ltr,"_enh",sep=""),paste(ltr,"_non",sep=""))
  for(l.sub in ltr.subs){
  ##clustalo round 1
  
  command = paste('~/clustalo -i ',l.sub,'.fa -t DNA -o ',l.sub,'_omega.fa --threads=2 --output-order=tree-order --force',sep='')
  try(system(command))
  
  
  ##read for sequence filtering
  
  fasta = scan(paste(l.sub,'omega.fa',sep='_'), sep='\n', character())
  
  
  ##Merge sequences
  
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
  
  bp = strsplit(seq,split='')
  bp.mat = matrix(unlist(bp),nrow=length(bp),byrow=T)
  
  
  ##Count gap size in sliding window
  
  n = 5
  gap.size = matrix(nrow=nrow(bp.mat),ncol=ncol(bp.mat)-n+1)
  for (i in 1:ncol(gap.size)) {
    window = bp.mat[,i:(i+n-1)]
    gap.size[,i] = rowSums(window=='-')
  }
  
  
  ##Get frequent gaps
  
  gap.freq = colSums(gap.size==n)
  frequent = gap.freq>=0.97*length(bp)
  
  
  ##Remove sequences making gaps
  
  not.empty = gap.size[,frequent]==0
  spurious = rowSums(not.empty)>0
  filtered = seq[!spurious]
  
  
  #Write
  
  file = paste(l.sub,'filt.fa',sep='_')
  unlink(file)
  for (i in 1:length(filtered)) {
    write(paste('>',names(filtered)[i],sep=''),file,append=T)
    write(filtered[i],file,append=T)
  }
  
  
  ##clustalo round 2
  
  command = paste('~/clustalo -i ',l.sub,'_filt.fa -t DNA -o ',l.sub,'_omega2.fa --threads=2 --output-order=tree-order --force',sep='')
  try(system(command))
  
}

}
fimo=function(ltr){
  ltr.subs=c(paste(ltr,"_enh",sep=""),paste(ltr,"_non",sep=""))
  for(l.sub in ltr.subs){
    command = paste('~/meme/bin/fimo --o ',l.sub,'_fimo ~/motif_databases/JASPAR/JASPAR_CORE_2016_vertebrates.meme ',l.sub,'.fa',sep='')
    try(system(command))
    try(system(paste('mv ',l.sub,'_fimo/fimo.txt ./',l.sub,'_fimo.txt',sep='')))
    try(system(paste('rm -r ',l.sub,'_fimo/',sep='')))
}}
ame=function(ltr){
  ltr.subs=c(paste(ltr,"_enh",sep=""),paste(ltr,"_non",sep=""))
  
  command = paste('~/meme/bin/ame --o ',ltr,'_ame --control ',ltr.subs[2],'.fa  ' ,ltr.subs[1],'.fa',' ~/motif_databases/JASPAR/JASPAR_CORE_2016_vertebrates.meme',sep='')
  try(system(command))
  ##move and rename file
  try(system(paste('mv ',ltr,'_ame/ame.txt ./',ltr,'_ame.txt',sep='')))
  try(system(paste('rm -r ',ltr,'_ame/',sep='')))
  

}

##Clustalo 
for (ltr in rltr.list) clustalo(ltr)

##FIMO and AME of RLTRs
for (ltr in rltr.list) fimo(ltr)
for (ltr in rltr.list) ame(ltr)
