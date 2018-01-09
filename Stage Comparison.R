###########################################################
#
# this is for comparing gene expression across different stages 
#
###########################################################

test<-read.table("./Gene Expression/gene exp_clinical stages.txt",sep="\t",stringsAsFactors=F,header=T)

all.data<-read.table("./Gene Expression/gene exp_all batches.txt",sep="\t",stringsAsFactors=F,header=T)

temp<-intersect(rownames(test),colnames(all.data))

test<-test[temp,]

stage1<-cbind(all.data[,which(test[,2]=="Stage IA")],all.data[,which(test[,2]=="Stage IB")],all.data[,which(test[,2]=="Stage IC")])

# manual correction
colnames(stage1)[1]<-"C09_452012"

stage2<-cbind(all.data[,which(test[,2]=="Stage IIA")],all.data[,which(test[,2]=="Stage IIB")],all.data[,which(test[,2]=="Stage IIC")])

# manual correction
colnames(stage2)[1]<-"F12_491286"

stage3<-cbind(all.data[,which(test[,2]=="Stage IIIA")],all.data[,which(test[,2]=="Stage IIIB")],all.data[,which(test[,2]=="Stage IIIC")])

stage4<-all.data[,which(test[,2]=="Stage IV")]

stages.test<-cbind(stage1,stage2,stage3,stage4)
test<-test[colnames(stages.test),]
stages.order<-c(rep("stage1",ncol(stage1)),rep("stage2",ncol(stage2)),rep("stage3",ncol(stage3)),rep("stage4",ncol(stage4)))

res<-NULL
for (i in 1:nrow(stages.test)){
	x<-kruskal.test(as.numeric(stages.test[i,])~as.factor(as.character(stages.order)))$p.value
	res<-c(res,x)
	}

names(res)<-rownames(stages.test)
res.sub<-res[which(res<0.05)]

overexp<-read.table("./Gene Expression/gene exp_batch 9_tumor/final hits_overexp_batch 9.txt", stringsAsFactors=F, sep="\t", header=T)
intersect(names(which(stage.kruskal.test<0.05)),rownames(overexp))

underexp<-read.table("./Gene Expression/gene exp_batch 9_tumor/final hits_underexp_batch 9.txt")
intersect(names(which(stage.kruskal.test<0.05)),rownames(underexp))

# alternative method

test.2<-test[colnames(all.data),]
test.2[,"stage"]<-gsub("(A|B|C)", "", test.2[,2])

stage.kruskal.test<-apply(data.matrix(stages.test), 1, function(x) kruskal.test(x~as.factor(test.2[colnames(stages.test),3]))$p.value)

hits<-stage.kruskal.test[which(stage.kruskal.test<0.05)]
hits<-cbind(unlist(mget(names(hits),hgu133aSYMBOL)),hits)
colnames(hits)<-c("Gene","P.Value")
write.table(hits,file="./Gene Expression/gene exp_stages comparison.txt",sep="\t",quote=F)

overexpress.kruskal.test<-apply(data.matrix(stages.test[rownames(overexp), ]), 1, function(x) kruskal.test(x~as.factor(test.2[colnames(stages.test),3]))$p.value)
hits.over<-overexpress.kruskal.test[which(overexpress.kruskal.test<0.05)]
hits.over<-cbind(unlist(mget(names(hits.over), hgu133aSYMBOL)),hits.over)
colnames(hits.over)<-c("Gene","P.Value")
write.table(hits.over,file="./Gene Expression/gene exp_stages comparison_overexp.txt",sep="\t",quote=F)


underexpress.kruskal.test<-apply(data.matrix(stages.test[rownames(underexp), ]), 1, function(x) kruskal.test(x~as.factor(test.2[colnames(stages.test),3]))$p.value)
hits.under<-underexpress.kruskal.test[which(underexpress.kruskal.test<0.05)]
hits.under<-cbind(unlist(mget(names(hits.under), hgu133aSYMBOL)),hits.under)
colnames(hits.under)<-c("Gene","P.Value")
write.table(hits.under,file="./Gene Expression/gene exp_stages comparison_underexp.txt",sep="\t",quote=F)