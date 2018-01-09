###########################################################
#
# this is the script for the weighted voting algorithm
#
###########################################################

#read samples
int<-read.table("./Gene Expression/PAMR/PAMR genes_underexp_batch 9.txt", stringsAsFactors=F, header=T, sep="\t")
# int<-read.table("./Gene Expression/gene exp_batch 9_tumor/final hits_underexp_batch 9.txt", stringsAsFactors=F, header=T, sep="\t")
# rownames(int)<-int[,1]
all.data<-read.table("./Gene Expression/gene exp_all batches.txt", stringsAsFactors=F, header=T, sep="\t")
all.data<-all.data[rownames((int)),]
source("./Scripts/data_compilation_script_normal.R")
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
# write.table(all.pred, file="./Gene Expression/PAMR/wva_PAMR prediction_underexp_batch 9.txt", quote=F, sep="\t")
# write.table(all.pred, file="./Gene Expression/wva_prediction_underexp_batch 9.txt", quote=F, sep="\t")


# #read samples
# source("./Scripts/data_compilation_script.R")
# source("./Scripts/data_compilation_script_normal.R")
# int<-read.table("./Gene Expression/intersect_overexp.txt",stringsAsFactors=F,header=T,sep="\t")

# b9.data<-b9.data[rownames(int),]
# b9.data.normal<-b9.data.normal[rownames(int),]
# b9<-cbind(b9.data, b9.data.normal)

# #compute signal-to-noise statistic 
# Sx<-apply(b9, 1, function(x) {
	# (mean(x[1:ncol(b9.data)], na.rm=T) - mean(x[(ncol(b9.data)+1):ncol(b9)], na.rm=T))/(sd(x[1:ncol(b9.data)], na.rm=T) + sd(x[(ncol(b9.data)+1):ncol(b9)], na.rm=T))
	# }
	# )

# #compute decision boundaries (half way) between the class means
# b9.data.mean<-apply(b9.data,1,mean,  na.rm=T)
# b9.data.normal.mean<-apply(b9.data.normal,1,mean, na.rm=T)
# Bx<-(b9.data.mean+b9.data.normal.mean)/2

# #predict class of test sample y
# Vx<-Sx*(b9.data-Bx)
# sigmaVx<-apply(Vx,2,sum)
# final.class<-sign(sigmaVx)

# #confidence in the prediction of winning class 
# conf<-apply(Vx,2,function(x){
	# abs((sum(x[sign(x)==1])-abs(sum(x[sign(x)==-1])))/(sum(x[sign(x)==1])+abs(sum(x[sign(x)==-1]))))
	# }
	# )
	
# b9.pred<-cbind(final.class,conf)	
# write.table(b9.pred, file="./Gene Expression/Gene Overexp Predictions/prediction_batch 9.txt", quote=F, sep="\t")