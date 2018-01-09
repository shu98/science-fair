read.table("./Clinical Info/clinical info_batch 11/nationwidechildrens.org_clinical_patient_ov.txt",sep="\t",stringsAsFactors=F,header=T)->clin11info

clin11info<-clin11info[-c(1:2),]
rownames(clin11info)<-clin11info[,2]
clin11info<-clin11info[,c(4,8,14,17:19,38)]

clin11info.dead<-clin11info[clin11info[,4]=="Dead",]
dead.below.clin<-clin11info.dead[as.numeric(clin11info.dead[,6])<=1825,]
dead.above.clin<-rbind(clin11info.dead[as.numeric(clin11info.dead[,6])>=1825,],clin11info[as.numeric(clin11info[,5])>=1825,])
dead.above.clin<-unique(dead.above)

#manually remove NA
is.na(rownames(dead.below))
is.na(rownames(dead.above))

source("./Scripts/data_compilation_script.R")

dead.below<-b11.data[,rownames(dead.below)]
dead.above<-b11.data[,rownames(dead.above)]

y<-NULL
for (i in 1:nrow(dead.below)){
	x<-as.numeric(dead.below[i,])
	w<-as.numeric(dead.above[i,])
	res<-t.test(x,w)$p.value
	y<-c(y,res)
	}
	
hits<-y<0.05 & regexpr("AFFX", rownames(b11.data))==-1
b11surv.pvalue<-y
names(b11surv.pvalue)<-rownames(dead.below)
final.hits.surv<-b11surv.pvalue[hits==TRUE]

temp<-NULL
for (i in 1:nrow(dead.below)){
 	a<-as.numeric(dead.below[i,])
 	b<-as.numeric(dead.above[i,])
 	fc<-mean(2^a)/mean(2^b)
 	temp<-c(temp,fc)
 	}
 	
b11surv.fc<-temp
names(b11surv.fc)<-rownames(dead.below)
final.fc.surv<-b11surv.fc[b11surv.fc>2]

fc.t.res<-intersect(names(final.hits.surv), names(final.fc.surv))
final.hits.info<-cbind(final.fc.surv[fc.t.res], final.hits.surv[fc.t.res])
colnames(final.hits.info)[1:2]<-c("Fold Change", "P Value")

final.hits.info<-cbind(unlist(mget(rownames(final.hits.info), hgu133aSYMBOL)), final.hits.info)
colnames(final.hits.info)[1]<-"Gene Name"

b11all<-cbind(dead.below,dead.above)
tmp<-b11all[rownames(final.hits.info),]

final.hits.info<-cbind(final.hits.info,tmp)
b11surv<-final.hits.info
write.table(b11surv,file="./OV data_batch 11_tumor/new_survival data_batch 11_underexp",row.names=T, col.names=T, quote=F, sep="\t")
