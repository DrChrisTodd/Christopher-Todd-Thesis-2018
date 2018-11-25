setwd('~/MEME_analysis/RLTR')


library(ComplexHeatmap,lib= '~/Rpackages/')
library(circlize,lib= '~/Rpackages/')
heatmap.alignment=function(class){
  
  ltr.enh=paste(class,"enh",sep="_")
  ltr.non=paste(class,"non",sep="_")
  fasta.enh = scan(paste(ltr.enh,'omega2.fa',sep='_'), sep='\n', character()) 
  fasta.non = scan(paste(ltr.non,'omega2.fa',sep='_'), sep='\n', character()) 
  
  
  ##Merge sequences
  get.mer.seq=function(fasta){
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
    return(seq)}
  enh.seq=get.mer.seq(fasta.enh)
  non.seq=get.mer.seq(fasta.non)
  
  get.heatmap.df=function(seq){
    
    bp = strsplit(seq,split='')
    seq.length = length(bp[[1]])
    bp.count = lapply(bp, function(x) cumsum(!grepl('-',x)))
    
    mat.row=list()
    for(i in 1:length(bp)){
      txt = bp[[i]]
      bp.at=(txt!="-")*1
      mat.row[[i]]=bp.at
    }
    
    df=do.call(rbind,mat.row)
    row.names(df)=names(bp)
    return(df)
  }
  enh.heat=get.heatmap.df(enh.seq)
  non.heat=get.heatmap.df(non.seq)
  summary(as.factor(enh.heat))
  #put heatmaps together
  max.col=max(ncol(enh.heat),ncol(non.heat))  
  if(ncol(enh.heat)<ncol(non.heat)){
    fill=matrix(rep(0),ncol=max.col-ncol(enh.heat),nrow = nrow(enh.heat))
    enh.heat=cbind(enh.heat,fill)}
  if(ncol(enh.heat)>ncol(non.heat)){
    fill=matrix(rep(0),ncol=max.col-ncol(non.heat),nrow = nrow(non.heat))
    non.heat=cbind(non.heat,fill)}
  break.mat=matrix(rep(0),nrow=(100),ncol=(max.col))
  tog.mat=rbind(enh.heat,break.mat,non.heat)  
  cols=colorRamp2(c(0,1),c("white","gray17"))
  
  
  hm.plot= Heatmap(tog.mat,cluster_rows=F,cluster_columns=F,
                   show_row_names=F,show_column_names=F,show_heatmap_legend=F,col=cols)
  
  width=5
  height=5
  jpeg(paste(class,"alignment.jpeg",sep="_"),w=width,h=height,units='in',res=600)
  draw(hm.plot)
  dev.off()
}

class.list=c("RLTR9E","RLTR9D","RLTR9A3","RLTR13B1","RLTR13B2","RLTR13B3","RLTR13B4","RLTR13D5","RLTR13D6")
for(d in 1:length(class.list)){
  class=class.list[d]
  heatmap.alignment(class)
  }
  


