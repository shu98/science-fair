# Script for comparing weighted voting algorithm results

#
# load tumor gene expression predictions 
#

fs<-list.files(pattern="prediction", path="./Gene Expression/Gene underexp Predictions/")
fs<-fs[-9]

temp<-read.table(paste("./Gene Expression/Gene underexp Predictions/", fs[1], sep=""), sep="\t", header=F, stringsAsFactors=F, skip=1)

for(f in fs[-1]) 
 {
 	f1<-read.table(paste("./Gene Expression/Gene Underexp Predictions/", f, sep=""), sep="\t", header=F, stringsAsFactors=F, skip=1)
 	temp<-rbind(temp,f1)
 }
 
rownames(temp)<-temp[,1]
temp<-temp[,-1]
colnames(temp)<-c("final.class","confidence")
 
gene.exp.pred<-temp

#
# load normal gene expression predictions 
#

temp<-read.table("./Gene Expression/Gene Underexp Predictions/prediction_batch 9_normal.txt", sep="\t", header=F, stringsAsFactors=F, skip=1)

rownames(temp)<-temp[,1]
temp<-temp[,-1]
colnames(temp)<-c("final.class","confidence")

gene.exp.normal.pred<-temp

#
# load tumor miRNA expression predictions 
#

fs<-list.files(pattern="prediction", path="./miRNA/miRNA Underexp Predictions/")
fs<-fs[-9]

temp<-read.table(paste("./miRNA/miRNA Underexp Predictions/", fs[1], sep=""), sep="\t", header=F, stringsAsFactors=F, skip=1)

for(f in fs[-1]) 
 {
 	f1<-read.table(paste("./miRNA/miRNA Underexp Predictions/", f, sep=""), sep="\t", header=F, stringsAsFactors=F, skip=1)
 	temp<-rbind(temp,f1)
 }
 
rownames(temp)<-temp[,1]
temp<-temp[,-1]
colnames(temp)<-c("final.class","confidence")

miRNA.pred<-temp

#
# load normal miRNA expression predictions 
#

temp<-read.table("./miRNA/miRNA Underexp Predictions/prediction_batch 9_normal.txt", sep="\t", header=F, stringsAsFactors=F, skip=1)

rownames(temp)<-temp[,1]
temp<-temp[,-1]
colnames(temp)<-c("final.class","confidence")

miRNA.normal.pred<-temp

#
# load tumor DNA methylation predictions 
#

fs<-list.files(pattern="prediction", path="./DNA Methylation/DNA Methylation Underexp Predictions/")
fs<-fs[-9]

temp<-read.table(paste("./DNA Methylation/DNA Methylation Underexp Predictions/", fs[1], sep=""), sep="\t", header=F, stringsAsFactors=F, skip=1)

for(f in fs[-1]) 
 {
 	f1<-read.table(paste("./DNA Methylation/DNA Methylation Underexp Predictions/", f, sep=""), sep="\t", header=F, stringsAsFactors=F, skip=1)
 	temp<-rbind(temp,f1)
 }
 
rownames(temp)<-temp[,1]
temp<-temp[,-1]
colnames(temp)<-c("final.class","confidence")

DNA.methyl.pred<-temp

#
# load normal DNA methylation predictions 
#

temp<-read.table("./DNA Methylation/DNA Methylation Underexp Predictions/prediction_batch 9_normal.txt",sep="\t",header=F,stringsAsFactors=F,skip=1)

rownames(temp)<-temp[,1]
temp<-temp[,-1]
colnames(temp)<-c("final.class","confidence")

DNA.methyl.normal.pred<-temp

#
# calculate sensitivity and specificity
#

sens<-c(sum(gene.exp.pred[,1]==1)/nrow(gene.exp.pred),sum(miRNA.pred[,1]==1)/nrow(miRNA.pred),sum(DNA.methyl.pred[,1]==1)/nrow(DNA.methyl.pred))

names(sens)<-c("Gene.Exp.Under","miRNA.Under","DNA Methyl.Under")

spec<-c(sum(gene.exp.normal.pred[,1]==-1)/sum(gene.exp.normal.pred[,1]==-1),sum(miRNA.normal.pred[,1]==-1)/sum(miRNA.normal.pred[,1]==-1), sum(DNA.methyl.normal.pred[,1]==-1)/sum(DNA.methyl.normal.pred[,1]==-1))

names(spec)<-c("Gene.Exp.Under","miRNA.Under","DNA Methyl.Under")

accuracy.under<-rbind(sens, spec)

rownames(accuracy.under)<-c("sensitivity","specificity")

#
# calculate average confidences for each prediction type
#

conf<-c(apply(gene.exp.pred,2,mean)[2],apply(miRNA.pred,2,mean)[2],apply(DNA.methyl.pred,2,mean)[2])

names(conf)<-c("Gene.Exp.Under","miRNA.Under","DNA Methyl.Under")

conf.norm<-c(apply(gene.exp.normal.pred,2,mean)[2],apply(miRNA.normal.pred,2,mean)[2],apply(DNA.methyl.normal.pred,2,mean)[2])

names(conf.norm)<-c("Gene.Exp.Under","miRNA.Under","DNA Methyl.Under")

conf.under<-rbind(conf,conf.norm)

rownames(conf.under)<-c("confidence.tumor","confidence.normal")
