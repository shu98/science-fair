source("./Scripts/data_compilation_script_miR.R")
source("./Scripts/data_compilation_script_normal_miR.R")

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

#
# add miRNA annotations 
#
miR<-read.table("./Downloaded Data/h-miR_V2.txt",sep="\t",header=T,stringsAsFactors=F)
miR<-miR[,c(1,2,5)]
rownames(miR)<-miR[,1]
miR<-miR[,-1]
int<-intersect(rownames(miR),rownames(final.hits.info))
miR<-miR[int,]
final.hits.info<-cbind(miR,final.hits.info)

b9.all<-cbind(b9.data, b9.data.normal)
tmp<-b9.all[rownames(final.hits.info),]

final.hits.info<-cbind(final.hits.info,tmp)

b9<-final.hits.info

# compute mean of tumor and normal 
b9.mean<-cbind(apply(b9.data, 1, mean, na.rm=T), apply(b9.data.normal, 1, mean, na.rm=T))
colnames(b9.mean)<-c("Tumor Mean miRNA Level","Normal Mean miRNA Level")
b9<-cbind(b9,b9.mean[rownames(b9),])

write.table(b9, file="./miRNA/miRNA_batch 9_tumor/final hits_underexp_batch 9.txt",quote=F, sep="\t", row.names=T, col.names=T)

load("./Downloaded Data/hsa_mir.Rdata")
mir<-hsa.mir[,c(1,2,4,5)]
mir.int<-b9[,1:2]
# mir.int<-read.table("./miRNA/intersect_miR_underexp.txt",sep="\t",stringsAsFactors=F,header=T)
int<-intersect(mir[,2],mir.int[,1])
mir.tar<-mir.int
mir.tar$Target<-mir.tar[,1] %in% int
mir$Target<-mir[,2] %in% int
mir.tar.true<-mir.tar[which(mir.tar[,"Target"]==T),]
mir.true<-mir[which(mir[,"Target"]==T),]
mir.true<-unique(mir.true)
# write.table(mir.true, file="./miRNA/mirtarbase_underexp.txt",quote=F,sep="\t")