#
# set up data
#
clin12<-read.table("./Clinical Info/clinical info_batch 12/nationwidechildrens.org_clinical_patient_ov.txt",sep="\t",stringsAsFactors=F,skip=2,header=T)
temp<-clin12[,2]
clin12<-clin12[,38]
names(clin12)<-temp
# manual correction to eliminate N/As
# clin12<-clin12[-19]

source("./Scripts/data_compilation_script.R")

stage3c<-b12.data[,names(which(clin12=="Stage IIIC"))]
stage4<-b12.data[,names(which(clin12=="Stage IV"))]

#
# t.test to compare tumor and normal
#
y<-NULL
for (i in 1:nrow(stage4)){
 	x<-as.numeric(stage4[i,])
 	w<-as.numeric(stage3c[i,])
 	res<-t.test(x,w)$p.value
 	y<-c(y,res)
 	}

hits<-y<0.05 & regexpr("AFFX", rownames(stage4))==-1
pvalue.stage12<-y
names(pvalue.stage12)<-rownames(stage4)
final.hits<-pvalue.stage12[hits==TRUE]

#
# fold change
#
temp<-NULL
for (i in 1:nrow(stage4)){
 	a<-as.numeric(stage4[i,])
 	b<-as.numeric(stage3c[i,])
 	fc<-mean(2^a)/mean(2^b)
 	temp<-c(temp,fc)
 	}
 
fold.change<-temp
names(fold.change)<-rownames(stage4)
final.fold.change<-fold.change[fold.change<0.5]

#
# intersect t.test hits and fold change hits
#
fc.t.res<-intersect(names(final.hits), names(final.fold.change))
final.hits.info<-cbind(final.fold.change[fc.t.res], final.hits[fc.t.res])
colnames(final.hits.info)[1:2]<-c("Fold Change", "P Value")

library(hgu133a.db)

final.hits.info<-cbind(unlist(mget(rownames(final.hits.info), hgu133aSYMBOL)), final.hits.info)
colnames(final.hits.info)[1]<-"Gene Name"

stage12.all<-cbind(stage4, stage3c)
tmp<-stage12.all[rownames(final.hits.info),]

final.hits.info<-cbind(final.hits.info,tmp)

#
# add gene annotation
#

stage12<-final.hits.info

# compute mean of tumor and normal 
stage12.mean<-cbind(apply(stage4, 1, mean, na.rm=T), apply(stage3c, 1, mean, na.rm=T))
colnames(stage12.mean)<-c("Tumor Mean Exp Level","Normal Mean Exp Level")
stage12<-cbind(stage12,stage12.mean[rownames(stage12),])

write.table(stage12, file="./Gene Expression/gene exp_batch 12_tumor/stage comparison_underexp_batch 12.txt",quote=F, sep="\t", row.names=T, col.names=T)