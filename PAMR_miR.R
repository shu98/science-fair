# PAMR script for miRNA

library(pamr)
library(hgu133a.db)
library(annotate)


# read training files

source("./Scripts/data_compilation_script_miR.R")
source("./Scripts/data_compilation_script_normal_miR.R")

int.over<-read.table("./miRNA/miRNA_batch 9_tumor/final hits_overexp_batch 9.txt",sep="\t",stringsAsFactors=F,header=T)
int.over<-rownames(int.over)

#
# set up PAMR training samples
# 
x=data.matrix(cbind(b9.data[(int.over),],b9.data.normal[(int.over),]))
y=c(rep("tumor", ncol(b9.data)), rep("normal", ncol(b9.data.normal)))
geneid=rownames(x)

#
# PAMR training
#
mydata<-list(x=x, y=y, geneid=geneid)
mytrain<-pamr.train(mydata)
mycv<-pamr.cv(mytrain, mydata)
pamr.plotcv(mycv)

# assign optimal threshold 
pamSel.41<-pamr.listgenes(mytrain, mydata, threshold= 4.455)


# save table

all.results<-cbind(mycv$threshold,mycv$size,errors)
colnames(all.results)<-c("threshold","number.miRNA","errors")
write.table(all.results,file="./miRNA/pamr results_overexp_batch 9.txt",quote=F,sep="\t")

miR<-read.table("./Downloaded Data/h-miR_V2.txt",sep="\t",header=T,stringsAsFactors=F)
miR<-miR[,c(1,2,5)]
rownames(miR)<-miR[,1]
miR<-miR[,-1]
mir.tar<-miR[pamSel.41[,1],1]
names(mir.tar)<-pamSel.41[,1]
write.table(mir.tar,file="./miRNA/PAMR miRNA_overexp_batch 9.txt",quote=F,sep="\t")

# find target gene IDs 
load("./Downloaded Data/hsa_mir.Rdata")
mir<-hsa.mir[,c(1,2,4,5)]
int.over<-read.table("./miRNA/intersect_miR_overexp.txt",sep="\t",stringsAsFactors=F,header=T)
int<-intersect(pamSel.41[,1],rownames(int.over))
int.over<-int.over[int,]
mir.tar<-int.over
mir.tar$Target<-rownames(mir.tar) %in% int
mir$Target<-mir[,2] %in% mir.tar[,1]
mir.tar.true<-mir.tar[which(mir.tar[,"Target"]==T),]
mir.true<-mir[which(mir[,"Target"]==T),]
mir.true<-unique(mir.true)
write.table(mir.true, file="./temp.txt",quote=F,sep="\t")

#
# set up prediction testing samples
#
all.data<-read.table("./miRNA/miRNA_all batches.txt",sep="\t",stringsAsFactors=F,header=T)
all.data.normal<-b9.data.normal
rownames(all.data)<-all.data[,1]
all.data<-all.data[,-1]
a=data.matrix(cbind(all.data[(int.over),],all.data.normal[(int.over),]))
b=c(rep("tumor", ncol(all.data)), rep("normal", ncol(all.data.normal)))
geneid2=rownames(a)

#
# predictions
#
mydata2<-list(x=a, y=b, geneid=geneid2)
predict<-pamr.predict(mytrain, mydata2$x, threshold= 4.455)

sum(predict=="normal")
which(predict=="normal")