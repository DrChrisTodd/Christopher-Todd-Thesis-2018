setwd('~/R/Thesis Scripts/Results Chapter 1/Distance to nearest gene/')

###lookup files
chrm.look=read.delim("~/R/Thesis Scripts/Data/Annotations/chromosome_lengths.txt")
refseq.lookup="~/R/Thesis Scripts/Data/Annotations/RefSeq_mRNA_lookup.txt"
ROI=read.delim("~/R/Thesis Scripts/Data/Output/ROI_coord.txt")
seqmonk.genes=read.delim("~/R/Thesis Scripts/Data/Annotations/seqmonk genes.txt")
ROI_shuffle=read.delim("~/R/Thesis Scripts/Data/Output/ROI_shuffle.txt",h=F)


#prepare chromosome lookup files
len=c()
for(i in 1:nrow(chrm.look)){
  sub=as.character(chrm.look$Total.length..bp.)[i]
  x=unlist(strsplit(sub,split = ","))
  y=as.numeric(paste(x,collapse = ""))
  len[i]=y
}
chr=paste("chr",chrm.look$Chromosome,sep="")             
chrm.look=data.frame(chr,len)
seqmonk.genes$chrom=paste("chr",seqmonk.genes$Chromosome,sep="")
ROI$mid=(ROI[,7]+ROI[,8])/2
ROI_shuffle$mid=(ROI_shuffle[,2]+ROI_shuffle[,3])/2
ROI_shuffle$class=paste(substr(ROI_shuffle[,4],1,3),unlist(lapply(strsplit(as.character(ROI_shuffle[,4]),split = "_"),function(x){x[2]})),sep = "_")
class.list=as.character(unique(ROI$class))


#functions
gen.random.chr.match=function(input.group,chr.col){
  for(i in 1:length(unique(input.group[,chr.col]))){
    chr=as.character(unique(input.group[,chr.col])[i])
    n.chr=nrow(input.group[input.group[,chr.col]==chr,])
    chr.len=chrm.look$len[chrm.look$chr==chr]
    ran.pos=sample(1:chr.len,n.chr)
    df=data.frame(rep(chr),ran.pos)
    if(i==1){tog.df=df}else{tog.df=rbind(tog.df,df)}
  }
  tog.df$id=paste("rand_pos_",seq(1:nrow(tog.df)),sep="")
  return(tog.df)
}
get.nearest.gene=function(input.group,chr.col,mid.col){
  refseq= read.delim(refseq.lookup)
  nearest.gene=c()
  distance=c()
  for(i in 1:nrow(input.group)){
    chr.sub=refseq[refseq$chrom %in% input.group[i,chr.col],]
    nearest=chr.sub[which.min(abs(input.group[i,mid.col]-chr.sub$TSS)),]
    distance[i]=abs(input.group[i,mid.col]-nearest$TSS)
    if(nearest$mRNA_name!="null"){nearest.gene[i]=as.character(nearest$mRNA_name)}else{nearest.gene[i]=as.character(nearest$refseq_name)}
  }
  df=data.frame(input.group,distance,nearest.gene)
  
  return(df)
}



#get distances to the nearest TSS
for(i in 1:length(class.list)){
  len=max(summary(ROI$class))
  sub=ROI[ROI$class==class.list[i],]
  shuffle.sub=ROI_shuffle[ROI_shuffle$class==class.list[i],]
  sub.near=get.nearest.gene(input.group = sub,chr.col = 1,mid.col = 10)
  sub.rand.near=get.nearest.gene(input.group = shuffle.sub,chr.col = 1,mid.col = 5)
  sub.dist=append(sub.near$distance,rep("",len-nrow(sub)))
  sub.rand.dist=append(sub.rand.near$distance,rep("",len-nrow(sub)))
  df.distance=data.frame(class=sub.dist,rand=sub.rand.dist)
  df.nearest.gene=data.frame(class=sub.near$nearest.gene,rand=sub.rand.near$nearest.gene)
  colnames(df.distance)[1]=paste(class.list[i],"distance",sep="_")
  colnames(df.distance)[2]=paste(class.list[i],"randcontrl_distance",sep="_")
  colnames(df.nearest.gene)[1]=paste(class.list[i],"nearestgene",sep="_")
  colnames(df.nearest.gene)[2]=paste(class.list[i],"randcontrlnearestgene",sep="_")
  if(i==1){tog.df=df.distance}else{tog.df=cbind(tog.df,df.distance)}
}
write.table(tog.df,"ROI_distances_to_nearest_TSS.txt",sep="\t",quote = F,col.names = T,row.names = F)


#function for getting number of intronic instances
get.intronic.numbers=function(input.df,chr.col,mid.col,class.col){
  ##figuring out if ROI is intronic or not
  input.df=input.df[,c(chr.col,mid.col,class.col)]
  colnames(input.df)=c("chr","mid","class")
  chr.col=1
  mid.col=2
  class.col=3
  class.list=as.character(unique(input.df[,class.col]))
  
  #Test if intronic
  intronic=c()
  for(i in 1:nrow(input.df)){
    chr.sub=seqmonk.genes[seqmonk.genes$chrom %in% input.df[i,chr.col],]
    any.intron<-sum((input.df[i,mid.col]>chr.sub$Start&input.df[i,mid.col]<chr.sub$End))>1
    intronic[i]=any.intron
  }
  input.df$intronic=intronic
  
  
  for(i in 1:length(unique(input.df[,class.col]))){
    class=as.character(unique(input.df[,class.col]))[i]
    sub=input.df[input.df[,class.col] %in% class, ]
    in.gene=sum(sub$intronic==T)
    out.gene=sum(sub$intronic==F)
    row=data.frame(class=class,in.gene=in.gene,out.gene=out.gene)
    if(i==1){df=data.frame(row)}else{df=rbind(df,row)}
  }
  
  df[,1]=as.character(unique(input.df[,class.col]))
  return(df)}


roi.intronic=get.intronic.numbers(input.df = ROI,chr.col = 1,mid.col = 10,class.col = 9)
ROI_shuffle$temp_class=rep("Rand")
rand.intronic=get.intronic.numbers(input.df = ROI_shuffle,chr.col = 1,mid.col = 5,class.col = 7)
roi.intronic=rbind(roi.intronic,rand.intronic)
write.table(roi.intronic,"Number_intronic.txt",sep = "\t",quote = F,col.names = T,row.names = F)
