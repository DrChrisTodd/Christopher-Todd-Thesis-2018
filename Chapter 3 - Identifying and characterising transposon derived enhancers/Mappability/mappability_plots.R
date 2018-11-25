setwd('~/R/Thesis Scripts/Results Chapter 1/Mappability/')
library(ggplot2)

te.look=read.delim("~/R/Thesis Scripts/Data/Output/TE_coordinates_lookup.txt")[,c(8,7,5)]

##Heatmap orders 
es.ord=read.delim("~/R/Thesis Scripts/Data/Homer_heatmaps/ESC_noclass_heatmap_order.txt",h=F)
ts.ord=read.delim("~/R/Thesis Scripts/Data/Homer_heatmaps/TSC_noclass_heatmap_order.txt",h=F)



#prepare es heatmap order file
colnames(es.ord)=c("ID_class","order")
mer=merge(es.ord,te.look,by="ID_class")

r9=mer[grep(mer$repName,pattern="RLTR9"),]
r13d6=mer[grep(mer$repName,pattern="RLTR13D6"),]


##plotting es te 
gg=ggplot(r9,aes(y=order,x=mappability))+geom_rect(aes(xmin=-0.05,xmax=0.5,ymin=0,ymax=Inf),fill="grey",alpha=0.5)+
  geom_point()+scale_x_reverse()+scale_y_reverse()+theme_void()+
  labs(y="Order of H3K27ac Signal",x="Mappability",title="RLTR9")
jpeg("RLTR9_mappability.jpg",w=2.5,h=4,units='in',res=600)
gg
dev.off()

jpeg("RLTR13D6_mappability.jpg",w=2.5,h=4,units='in',res=600)
gg=ggplot(r13d6,aes(y=order,x=mappability))+geom_rect(aes(xmin=-0.05,xmax=0.5,ymin=0,ymax=Inf),fill="grey",alpha=0.5)+
  geom_point()+scale_x_reverse()+scale_y_reverse()+theme_void()+
  labs(y="Order of H3K27ac Signal",x="Mappability",title="RLTR13D6")
gg
dev.off()

##prepare ts heatmap order
colnames(ts.ord)=c("ID_class","order")
mer=merge(ts.ord,te.look,by="ID_class")

r13b=mer[grep(mer$repName,pattern="RLTR13B"),]
r13d5=mer[grep(mer$repName,pattern="RLTR13D5"),]

#plot ts te 
gg=ggplot(r13b,aes(y=order,x=mappability))+geom_rect(aes(xmin=-0.05,xmax=0.5,ymin=0,ymax=Inf),fill="grey",alpha=0.5)+
  geom_point()+scale_x_reverse()+scale_y_reverse()+theme_void()+
  labs(y="Order of H3K27ac Signal",x="Mappability",title="RLTR9")

jpeg("RLTR13B_mappability.jpg",w=2.5,h=4,units='in',res=600)
gg
dev.off()

gg=ggplot(r13d5,aes(y=order,x=mappability))+geom_rect(aes(xmin=-0.05,xmax=0.5,ymin=0,ymax=Inf),fill="grey",alpha=0.5)+
  geom_point()+scale_x_reverse()+scale_y_reverse()+theme_void()+
  labs(y="Order of H3K27ac Signal",x="Mappability",title="RLTR9")

jpeg("Output/Mappability/RLTR13D5_mappability.jpg",w=2.5,h=4,units='in',res=600)
gg
dev.off()
