source("./Scripts/data_compilation_script_DNA methylation.R")
source("./Scripts/data_compilation_script_normal_DNA methylation.R")

#
# t.test to compare tumor and normal
#
y<-NULL
for (i in 1:nrow(b9.data)){
	x<-as.numeric(b9.data[i,])
	w<-as.numeric(b9.data.normal[i,])
	res<-t.test(x,w)$p.value
	y<-c(y,res)
	}

y<-p.adjust(y,method="BH")

hits<-y<0.05 & regexpr("AFFX", rownames(b9.data))==-1
pvalue.b9.data<-y
names(pvalue.b9.data)<-rownames(b9.data)
final.hits<-pvalue.b9.data[hits==TRUE]

#
# fold change
#
temp<-NULL
for (i in 1:nrow(b9.data)){
	a<-as.numeric(b9.data[i,])
	b<-as.numeric(b9.data.normal[i,])
	fc<-mean(2^a)/mean(2^b)
	temp<-c(temp,fc)
	}

fold.change<-temp
names(fold.change)<-rownames(b9.data)
final.fold.change<-fold.change[fold.change>1.25]

#
# intersect t.test hits and fold change hits
#
fc.t.res<-intersect(names(final.hits), names(final.fold.change))
final.hits.info<-cbind(final.fold.change[fc.t.res], final.hits[fc.t.res])
colnames(final.hits.info)[1:2]<-c("Fold Change", "P Value")

gene.ret<-read.table("./DNA Methylation/gene list.txt", sep="\t", header=T, stringsAsFactors=F)
final.hits.info<-cbind(gene.ret[rownames(final.hits.info),],final.hits.info)
colnames(final.hits.info)[1]<-"Gene Name"

b9.all<-cbind(b9.data, b9.data.normal)
tmp<-b9.all[rownames(final.hits.info),]

final.hits.info<-cbind(final.hits.info,tmp)

#
# add gene annotation
#

b9<-final.hits.info

# compute mean of tumor and normal 
b9.mean<-cbind(apply(b9.data, 1, mean, na.rm=T), apply(b9.data.normal, 1, mean, na.rm=T))
colnames(b9.mean)<-c("Tumor Mean Methylation Level","Normal Mean Methylation Level")
b9<-cbind(b9,b9.mean[rownames(b9),])

write.table(b9, file="./DNA Methylation/DNA methylation_batch 9_tumor/final hits_overexp_batch 9.txt",quote=F, sep="\t", row.names=T, col.names=T)
