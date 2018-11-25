###compiling all the different classes of elements I shall be interested in plotting
#TE coord
bed1=read.delim("Data/Output/TE_coordinates_lookup.txt")

#NonTE coord
bed2=read.delim("Data/Output/NonTE_Enhancer_coord.txt")

#NonEnhTE coord
bed3=read.delim("Data/Output/NonEnhTE_candidates.txt")


#fill in missing columns
bed1$score=rep(0)
colnames(bed2)=c("chr","start","end","ID_class","score","Plot.start","Plot.end")
bed2$strand=rep(".")
bed3$score=rep(0)


##ROI file format will have chr, plot.start, plot.end, ID_class, score, strand, peak.start, peak.end, class
bed1=bed1[,c(1,9:10,8,11,4,2,3)]
bed2=bed2[,c(1,6:7,4:5,8,2:3)]
bed3=bed3[,c(1,9:10,8,11,4,2,3)]

get.class=function(x){
  tissue=substr(as.character(x$ID_class),1,3)
  type=unlist(lapply(strsplit(as.character(x$ID_class),split = "_"), function(x){x[2]}))
  class=paste(tissue,type,sep="_")
  return(class)
  }


tog=rbind(bed1,bed2,bed3)
tog$class=get.class(tog)

write.table(tog,"Data/Output/ROI_coord.txt",sep="\t",quote = F,col.names = T,row.names = F)
