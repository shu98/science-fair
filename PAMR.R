library(pamr)
library(hgu133a2.db)
library(annotate)

#
# read training files
#
source("./Scripts/data_compilation_script.R")
source("./Scripts/data_compilation_script_normal.R")

int.under<-read.table("./Gene Expression/gene exp_batch 9_tumor/final hits_underexp_batch 9.txt",sep="\t",stringsAsFactors=F,header=T)
int.under<-rownames(int.under)

#
# set up PAMR training samples
#
x=data.matrix(cbind(b9.data[(int.under),],b9.data.normal[(int.under),]))
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
pamSel.41<-pamr.listgenes(mytrain, mydata, threshold=8.361)

#
# save table
#
# manually generate "errors" vector
all.results<-cbind(mycv$threshold,mycv$size,errors)
colnames(all.results)<-c("threshold","number.genes","errors")
write.table(all.results,file="./Gene Expression/PAMR results_underexp_batch 9.txt",quote=F,sep="\t")

#
# find target gene IDs 
#
write.table(getSYMBOL(pamSel.41[,1],"hgu133a.db"),file="./Gene Expression/PAMR genes_underexp_batch 9.txt",quote=F,sep="\t")

#
# set up prediction testing samples
#
all.data<-read.table("./Gene Expression/gene exp_all batches.txt", stringsAsFactors=F, sep="\t", header=T)
all.data.normal<-b9.data.normal
a=data.matrix(cbind(all.data[(int.under),],all.data.normal[(int.under),]))
b=c(rep("tumor", ncol(all.data)), rep("normal", ncol(all.data.normal)))
geneid2=rownames(a)

#
# predictions
#
mydata2<-list(x=a, y=b, geneid=geneid2)
predict<-pamr.predict(mytrain, mydata2$x, threshold=8.361)

sum(predict=="normal")
which(predict=="normal")
