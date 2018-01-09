#
# correlation between gene underexpression and DNA hypermethylation 
#

# read in files
underexp<-read.table("./Gene Expression/gene exp_batch 9_tumor/final hits_underexp_batch 9.txt",sep="\t",stringsAsFactors=F,header=T)
DNA.hyper<-read.table("./DNA Methylation/DNA methylation_batch 9_tumor/final hits_overexp_batch 9.txt",sep="\t",stringsAsFactors=F,header=T) 

# set up objects
# rownames(underexp)<-underexp[,1]
# underexp<-underexp[,-1]

int<-intersect(DNA.hyper[,1],underexp[,1])
underexp$DNA.hyper<-underexp[,1] %in% int
DNA.hyper$underexp<-DNA.hyper[,1] %in% int
underexp<-underexp[which(underexp[,"DNA.hyper"]==T),]
DNA.hyper<-DNA.hyper[which(DNA.hyper[,"underexp"]==T),]
underexp<-underexp[!duplicated(underexp[,1]),]
DNA.hyper<-DNA.hyper[!duplicated(DNA.hyper[,1]),]

DNA.hyper<-cbind(DNA.hyper,is.na(DNA.hyper[,1]))
DNA.hyper<-DNA.hyper[which(DNA.hyper[,"is.na(DNA.hyper[, 1])"]==F),]
underexp<-cbind(underexp,is.na(underexp[,1]))
underexp<-underexp[which(underexp[,"is.na(underexp[, 1])"]==F),]

rows<-read.table("./names_methylation_gene exp.txt",sep="\t",header=T,stringsAsFactors=F)
rows<-rows[-19,]
rows<-rows[,-1]
test.DNA.hyper<-DNA.hyper[,rows[,1]]
test.underexp<-underexp[,rows[,2]]

names.cor<-cbind(rownames(DNA.hyper),DNA.hyper[,1])
temp<-NULL
for (i in 1:nrow(DNA.hyper)){
	probe<-rownames(underexp)[which(underexp[,1]==DNA.hyper[i,1])]
	temp<-c(temp,probe)}
names.cor<-cbind(names.cor,temp)
colnames(names.cor)<-c("CpG","Gene","Probe")

# perform correlation 
temp<-data.frame(cor.value=character(0),p.value=character(0),stringsAsFactors=F)
for (i in 1:nrow(names.cor)){
	cor.value<-cor(as.numeric(test.DNA.hyper[names.cor[i,1],]),as.numeric(test.underexp[names.cor[i,3],]))
	p.value<-cor.test(as.numeric(test.DNA.hyper[names.cor[i,1],]),as.numeric(test.underexp[names.cor[i,3],]))$p.value
	x<-NULL
	x<-c(cor.value,p.value)
	temp<-rbind(temp,x)}

# which(underexp[,1]=="MEST")
# rownames(underexp)[3]

# DNA.hyper[which(DNA.hyper[,1]=="MEST"),]

# cor(as.numeric(test.DNA.hyper["cg01888566",]),as.numeric(test.underexp["202016_at",]))

# write.table("./Gene Expression/gene underexp_DNA hypermethylation_correlations.txt")

# save data
colnames(temp)<-c("pearson.coeff","p.value")
temp<-cbind(names.cor,temp)
# write.table(temp,file="./Gene Expression/gene underexp_DNA hypermethylation_correlations.txt", quote=F, sep="\t")