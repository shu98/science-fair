###########################################################
#
# this is the script for the weighted voting algorithm for DNA methylation
#
###########################################################

#read samples
# int<-read.table("./DNA Methylation/PAMR genes_underexp_batch 9.txt", stringsAsFactors=F, header=T, sep="\t")
# rownames(int)<-int[,2]
# int<-read.table("./DNA Methylation/DNA methylation_batch 9_tumor/final hits_overexp_batch 9.txt", stringsAsFactors=F, header=T, sep="\t")
# all.data<-read.table("./DNA Methylation/DNA methylation_all batches.txt", stringsAsFactors=F, header=T, sep="\t")
# all.data<-all.data[rownames((int)),]
# source("./Scripts/data_compilation_script_DNA methylation.R")
# source("./Scripts/data_compilation_script_normal_DNA methylation.R")
# all.data.normal<-all.data.normal[rownames((int)),]
# all<-cbind(all.data,all.data.normal)

# #compute signal-to-noise statistic 
# Sx<-apply(all, 1, function(x) {
	# (mean(x[1:ncol(all.data)], na.rm=T) - mean(x[(ncol(all.data)+1):ncol(all)], na.rm=T))/(sd(x[1:ncol(all.data)], na.rm=T) + sd(x[(ncol(all.data)+1):ncol(all)], na.rm=T))
	# }
	# )

# #compute decision boundaries (half way) between the class means
# all.data.mean<-apply(all.data,1,mean,  na.rm=T)
# all.data.normal.mean<-apply(all.data.normal,1,mean, na.rm=T)
# Bx<-(all.data.mean+all.data.normal.mean)/2

# #predict class of test sample y
# Vx<-Sx*(all-Bx)
# sigmaVx<-apply(Vx,2,sum)
# final.class<-sign(sigmaVx)

# #confidence in the prediction of winning class 
# conf<-apply(Vx,2,function(x){
	# abs((sum(x[sign(x)==1])-abs(sum(x[sign(x)==-1])))/(sum(x[sign(x)==1])+abs(sum(x[sign(x)==-1]))))
	# }
	# )
	
# all.pred<-cbind(final.class,conf)	
# write.table(all.pred, file="./DNA Methylation/wva_PAMR prediction_underexp_batch 9.txt", quote=F, sep="\t")
# write.table(all.pred, file="./DNA Methylation/wva_prediction_underexp_batch 9.txt", quote=F, sep="\t")

#read samples
source("./Scripts/data_compilation_script_DNA methylation.R")
source("./Scripts/data_compilation_script_normal_DNA methylation.R")
int<-read.table("./DNA Methylation/DNA methylation_batch 9_tumor/final hits_overexp_batch 9.txt", stringsAsFactors=F, header=T, sep="\t")

all.data<-all.data[rownames(int),]
all.data.normal<-all.data.normal[rownames(int),]
all<-cbind(all.data, all.data.normal)

#compute signal-to-noise statistic 
Sx<-apply(all, 1, function(x) {
	(mean(x[1:ncol(all.data)], na.rm=T) - mean(x[(ncol(all.data)+1):ncol(all)], na.rm=T))/(sd(x[1:ncol(all.data)], na.rm=T) + sd(x[(ncol(all.data)+1):ncol(all)], na.rm=T))
	}
	)

#compute decision boundaries (half way) between the class means
all.data.mean<-apply(all.data,1,mean)
all.data.normal.mean<-apply(all.data.normal,1,mean)
Bx<-(all.data.mean+all.data.normal.mean)/2

#predict class of test sample y
Vx<-Sx*(all.data-Bx)

#remove NA
Vx<-Vx[which(!is.na(Vx[,1])),]

sigmaVx<-apply(Vx,2,sum)
final.class<-sign(sigmaVx)

#confidence in the prediction of winning class 
conf<-apply(Vx,2,function(x){
	abs((sum(x[sign(x)==1])-abs(sum(x[sign(x)==-1])))/(sum(x[sign(x)==1])+abs(sum(x[sign(x)==-1]))))
	}
	)
	
all.pred<-cbind(final.class,conf)	
#write.table(all.pred, file="./DNA Methylation/DNA Methylation Underexp Predictions/prediction_batch 19.txt", quote=F, sep="\t")