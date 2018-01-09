###########################################################
#
# this is the script for the weighted voting algorithm for miRNA
#
###########################################################

# read samples
# int<-read.table("./miRNA/PAMR miRNA_underexp_batch 9.txt", stringsAsFactors=F, header=T, sep="\t")
int<-read.table("./miRNA/miRNA_batch 9_tumor/final hits_underexp_batch 9.txt", stringsAsFactors=F, header=T, sep="\t")
all.data<-read.table("./miRNA/miRNA_all batches.txt", stringsAsFactors=F, header=T, sep="\t")
rownames(all.data)<-all.data[,1]
all.data<-all.data[,-1]
all.data<-all.data[rownames((int)),]
source("./Scripts/data_compilation_script_normal_miR.R")
all.data.normal<-b9.data.normal[rownames((int)),]
all<-cbind(all.data,all.data.normal)

#compute signal-to-noise statistic 
Sx<-apply(all, 1, function(x) {
	(mean(x[1:ncol(all.data)], na.rm=T) - mean(x[(ncol(all.data)+1):ncol(all)], na.rm=T))/(sd(x[1:ncol(all.data)], na.rm=T) + sd(x[(ncol(all.data)+1):ncol(all)], na.rm=T))
	}
	)

#compute decision boundaries (half way) between the class means
all.data.mean<-apply(all.data,1,mean,  na.rm=T)
all.data.normal.mean<-apply(all.data.normal,1,mean, na.rm=T)
Bx<-(all.data.mean+all.data.normal.mean)/2

#predict class of test sample y
Vx<-Sx*(all-Bx)
sigmaVx<-apply(Vx,2,sum)
final.class<-sign(sigmaVx)

#confidence in the prediction of winning class 
conf<-apply(Vx,2,function(x){
	abs((sum(x[sign(x)==1])-abs(sum(x[sign(x)==-1])))/(sum(x[sign(x)==1])+abs(sum(x[sign(x)==-1]))))
	}
	)
	
all.pred<-cbind(final.class,conf)	
# write.table(all.pred, file="./miRNA/wva_PAMR prediction_underexp_batch 9.txt", quote=F, sep="\t")
# write.table(all.pred, file="./miRNA/wva_prediction_underexp_batch 9.txt", quote=F, sep="\t")


# #read samples
# source("./Scripts/data_compilation_script_miR.R")
# source("./Scripts/data_compilation_script_normal_miR.R")
# int<-read.table("./miRNA/intersect_miR_overexp.txt",stringsAsFactors=F,header=T,sep="\t")

# b19.data<-b19.data[rownames(int),]
# b19.data.normal<-b19.data.normal[rownames(int),]
# b19<-cbind(b19.data, b19.data.normal)

# #compute signal-to-noise statistic 
# Sx<-apply(b19, 1, function(x) {
	# (mean(x[1:ncol(b19.data)], na.rm=T) - mean(x[(ncol(b19.data)+1):ncol(b19)], na.rm=T))/(sd(x[1:ncol(b19.data)], na.rm=T) + sd(x[(ncol(b19.data)+1):ncol(b19)], na.rm=T))
	# }
	# )

# #compute decision boundaries (half way) between the class means
# b19.data.mean<-apply(b19.data,1,mean)
# b19.data.normal.mean<-apply(b19.data.normal,1,mean)
# Bx<-(b19.data.mean+b19.data.normal.mean)/2

# #predict class of test sample y
# Vx<-Sx*(b19.data-Bx)
# sigmaVx<-apply(Vx,2,sum)
# final.class<-sign(sigmaVx)

# #confidence in the prediction of winning class 
# conf<-apply(Vx,2,function(x){
	# abs((sum(x[sign(x)==1])-abs(sum(x[sign(x)==-1])))/(sum(x[sign(x)==1])+abs(sum(x[sign(x)==-1]))))
	# }
	# )

# b19.pred<-cbind(final.class,conf)	
# write.table(b19.pred, file="./miRNA/miRNA Overexp Predictions/prediction_batch 19.txt", quote=F, sep="\t")
