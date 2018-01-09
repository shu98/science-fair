source("./Scripts/data_compilation_script.R")
source("./Scripts/data_compilation_script_normal.R")

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
final.fold.change<-fold.change[fold.change>2]

#
# intersect t.test hits and fold change hits
#
fc.t.res<-intersect(names(final.hits), names(final.fold.change))
final.hits.info<-cbind(final.fold.change[fc.t.res], final.hits[fc.t.res])
colnames(final.hits.info)[1:2]<-c("Fold Change", "P Value")

library(hgu133a.db)

#
# add gene IDs
#
final.hits.info<-cbind(unlist(mget(rownames(final.hits.info), hgu133aSYMBOL)), final.hits.info)
colnames(final.hits.info)[1]<-"Gene Name"

b9.all<-cbind(b9.data, b9.data.normal)
tmp<-b9.all[rownames(final.hits.info),]

final.hits.info<-cbind(final.hits.info,tmp)

# #
# # intensity filter for overexpressed genes: remove genes not normally expressed at significant levels in cell
# #
intensity.flt<-NULL
for (i in 1:nrow(final.hits.info)) {
	intensity.flt<-c(intensity.flt, mean(as.numeric(final.hits.info[i, -c(1:3)]))>=log2(100))
	}
names(intensity.flt)<-rownames(final.hits.info)
final.hits.info<-final.hits.info[intensity.flt==TRUE,]

#
# intensity filter for underexpressed genes: remove genes not normally expressed at significant levels in cell
#
# intensity.flt<-NULL
# for (i in 1:nrow(final.hits.info)) {
	# intensity.flt<-c(intensity.flt, mean(as.numeric(final.hits.info[i, -c(1:41)]))>=log2(100))
	# }
# names(intensity.flt)<-rownames(final.hits.info)
# final.hits.info<-final.hits.info[intensity.flt==TRUE,]

b9<-final.hits.info

#
# compute mean of tumor and normal 
#
b9.mean<-cbind(apply(b9.data, 1, mean, na.rm=T), apply(b9.data.normal, 1, mean, na.rm=T))
colnames(b9.mean)<-c("Tumor Mean Exp Level","Normal Mean Exp Level")
b9<-cbind(b9,b9.mean[rownames(b9),])

# for (i in 1:nrow(b9)){
# temp<-as.numeric(b9[i,4:41])
# temp2<-mean(temp)
# b9.mean<-c(b9.mean,temp2)}
# b9<-cbind(b9,b9.mean)
# colnames(b9)[50]<-"Tumor Mean Exp Level"

# b9.mean.normal<-NULL
# for (i in 1:nrow(b9)){
# temp<-as.numeric(b9[i,42:49])
# temp2<-mean(temp)
# b9.mean.normal<-c(b9.mean.normal,temp2)}
# b9<-cbind(b9,b9.mean.normal)
# colnames(b9)[51]<-"Normal Mean Exp Level"

# write.table(b9, file="./Gene Expression/gene exp_batch 9_tumor/final hits_overexp_batch 9.txt",quote=F, sep="\t", row.names=T, col.names=T)

#
# find overrepresented pathways 
#
# library(hgu133a.db)
# all.paths<-sort(table(unlist(mget(int[,1], hgu133aPATH))), decreasing=T)

#
# annotate pathways 
#
# p53<-read.table("./Pathways/kegg_p53 signaling pathway.txt", header=T, stringsAsFactors=F,sep="\t")[,1]
# intersect(b9[,1],p53)
# b9$P53<-b9[,1] %in% p53
# b9[which(b9[,"P53"]==T),]

#
# test for similar genes between gene and miRNA expression and DNA methylation
#
# genes<-read.table("./Gene Expression/gene exp_batch 9_tumor/final hits_underexp_batch # 9.txt",sep="\t",stringsAsFactors=F,header=T)
# miRNA<-read.table("./miRNA/miRNA_batch 9_tumor/mirtarbase_underexp_batch 9.txt", sep="\t", header=T, stringsAsFactors=F)
# DNA<-read.table("./DNA Methylation/DNA methylation_batch 9_tumor/final hits_underexp_batch 9.txt",sep="\t",stringsAsFactors=F,header=T)
# intersect(genes[,1],DNA[,1])
# intersect(miRNA[,3],genes[,1])
# intersect(intersect(miRNA[,3],genes[,1]),intersect(genes[,1],DNA[,1]))