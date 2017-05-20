#!/usr/bin/Rscript
library(cluster)
args <- commandArgs(trailingOnly = TRUE) 
print(args)
file<-args[1]
htxt<-as.numeric(args[2])

inp = read.table(file,header=F,sep="\t",na.strings="-",fill=TRUE)
m = data.matrix(inp)
print(dim(m))
d = as.dist(t(m),upper=FALSE)
hc1 <- agnes(d,diss = TRUE,method="average")
pdf('cluster.pdf')
plot(hc1)
dev.off()

hcac=as.hclust(hc1)

memb=cutree(hcac,h=htxt)
#memb2 = cutree(hcac, h=0.15)
#memb3 = cutree(hcac, h=0.20)
write.table(memb,file=paste(file,"UPGMA","Cl",htxt,sep="_"),sep="\t", quote = FALSE,col.names = FALSE,row.names=FALSE)
#write.table(memb2,file=paste(file,"UPGMA","Cl",0.15,sep="_"),sep="\t", quote = FALSE,col.names = FALSE,row.names=FALSE)
#write.table(memb3,file=paste(file,"UPGMA","Cl",0.20,sep="_"),sep="\t", quote = FALSE,col.names = FALSE,row.names=FALSE)
