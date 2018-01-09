# PAMR script for DNA methylation

library(pamr)

#
# read training files
#
source("./Scripts/data_compilation_script_DNA methylation.R")
source("./Scripts/data_compilation_script_normal_DNA methylation.R")

int.under<-read.table("./DNA Methylation/DNA methylation_batch 9_tumor/final hits_overexp_batch 9.txt",sep="\t",stringsAsFactors=F,header=T)
int.under<-rownames(int.under)

#
# set up PAMR training samples
#
x=data.matrix(cbind(b9.data[(int.under),],b9.data.normal[(int.under),]))
y=c(rep("tumor", ncol(b9.data)), rep("normal", ncol(b9.data.normal)))
geneid=rownames(x)

#
# set up prediction testing samples
#
all.data<-read.table("./DNA Methylation/DNA methylation_all batches.txt", stringsAsFactors=F, sep="\t", header=T)
all.data.normal<-b9.data.normal
a=data.matrix(cbind(all.data[(int.under),],all.data.normal[(int.under),]))
b=c(rep("tumor", ncol(all.data)), rep("normal", ncol(all.data.normal)))
geneid2=rownames(a)

#
# remove genes with NA methylation values
#
x<-x[-which(apply(is.na(a),1,sum)>0),]
geneid=rownames(x)

a<-a[-which(apply(is.na(a),1,sum)>0),]
geneid=rownames(a)

#
# PAMR training
#
mydata<-list(x=x, y=y, geneid=geneid)
mytrain<-pamr.train(mydata)
mycv<-pamr.cv(mytrain, mydata)
pamr.plotcv(mycv)

# assign optimal threshold 
pamSel.41<-pamr.listgenes(mytrain, mydata, threshold= 2.737)

#
# save table
#
# manually generate "errors" vector
all.results<-cbind(mycv$threshold,mycv$size,errors)
colnames(all.results)<-c("threshold","number.genes","errors")
write.table(all.results,file="./DNA Methylation/PAMR results_overexp_batch 9.txt",quote=F,sep="\t")

#
# find target gene IDs 
#
gene.ret<-read.table("./DNA Methylation/gene list.txt", sep="\t", header=T, stringsAsFactors=F)
gene.ret<-cbind(pamSel.41[,1],gene.ret[pamSel.41[,1],])
write.table(gene.ret,file="./DNA Methylation/PAMR genes_overexp_batch 9.txt",quote=F,sep="\t")

#
# predictions
#
mydata2<-list(x=a, y=b, geneid=geneid2)
predict<-pamr.predict(mytrain, mydata2$x, threshold= 2.737)

sum(predict=="normal")
which(predict=="normal")