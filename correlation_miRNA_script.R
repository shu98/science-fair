#
# correlation between gene underexpression and miRNA overexpression
#
# read in files
underexp<-read.table("./Gene Expression/gene exp_batch 9_tumor/final hits_underexp_batch 9.txt",sep="\t",stringsAsFactors=F,header=T)
miRNA<-read.table("./miRNA/miRNA_batch 9_tumor/final hits_overexp_batch 9.txt",sep="\t",stringsAsFactors=F,header=T)
miRNA.tar<-read.table("./miRNA/miRNA_batch 9_tumor/mirtarbase_overexp_batch 9.txt",sep="\t",stringsAsFactors=F,header=T)

# set up objects
# underexp<-underexp[!duplicated(underexp[,1]),]
# underexp<-underexp[!is.na(underexp[,1]),]
# rownames(underexp)<-underexp[,1]
# underexp<-underexp[,-1]

int2<-intersect(miRNA.tar[,2],miRNA[,1])
miRNA.tar$miRNA.int<-miRNA.tar[,2] %in% int2
miRNA$miRNA.tar<-miRNA[,1] %in% int2
miRNA<-miRNA[which(miRNA[,"miRNA.tar"]==T),]
miRNA.tar<-miRNA.tar[which(miRNA.tar[,"miRNA.int"]==T),]
miRNA<-miRNA[!duplicated(miRNA[,1]),]

int<-intersect(miRNA.tar[,3],underexp[,1])
underexp$miRNA.tar<-underexp[,1] %in% int
miRNA.tar$underexp<-miRNA.tar[,3] %in% int
underexp<-underexp[which(underexp[,"miRNA.tar"]==TRUE),]
miRNA.tar<-miRNA.tar[which(miRNA.tar[,"underexp"]==T),]
underexp<-underexp[!duplicated(underexp[,1]),]

temp<-NULL
for (i in 1:nrow(miRNA.tar)){
	probe<-rownames(miRNA)[which(miRNA[,1]==miRNA.tar[i,2])]
	temp<-c(temp,probe)}
miRNA.tar<-cbind(miRNA.tar,temp)
colnames(miRNA.tar)[8]<-"miRNA.Probe"
names.cor<-miRNA.tar[,c(2,8)]

# int<-intersect(underexp[,1],miRNA.tar[,3])
# miRNA.tar$underexp<-miRNA.tar[,3] %in% int
# miRNA.tar<-miRNA.tar[which(miRNA.tar[,"underexp"]==T),]
# underexp$miRNA.Tar<-underexp[,1] %in% int
# underexp<-underexp[which(underexp[,"miRNA.Tar"]==T),]
# miRNA$miRNA.Tar<-miRNA[,1] %in% miRNA.tar[,2]
# miRNA<-miRNA[which(miRNA[,"miRNA.Tar"]==T),]

rows<-read.table("./names_miRNA_gene exp.txt",sep="\t",header=T,stringsAsFactors=F)
rows<-rows[-19,]
rows<-rows[,-1]
test.miRNA<-miRNA[,rows[,1]]
test.underexp<-underexp[,rows[,2]]

names.cor<-cbind(names.cor,miRNA.tar[,3])
temp<-NULL
for (i in 1:nrow(miRNA.tar)){
	probe<-rownames(underexp)[which(underexp[,1]==miRNA.tar[i,3])]
	temp<-c(temp,probe)}
names.cor<-cbind(names.cor,temp)
colnames(names.cor)[3:4]<-c("Gene","Probe")

# perform correlation 
temp<-data.frame(cor.value=character(0),p.value=character(0),stringsAsFactors=F)
for (i in 1:nrow(names.cor)){
	cor.value<-cor(as.numeric(test.miRNA[names.cor[i,2],]),as.numeric(test.underexp[names.cor[i,4],]))
	p.value<-cor.test(as.numeric(test.miRNA[names.cor[i,2],]),as.numeric(test.underexp[names.cor[i,4],]))$p.value
	x<-NULL
	x<-c(cor.value,p.value)
	temp<-rbind(temp,x)}

# which(underexp[,1]=="UCHL1")
# rownames(underexp)[1]

# miRNA.tar[which(miRNA.tar[,3]=="UCHL1"),]
# rownames(miRNA[which(miRNA[,1]=="hsa-miR-424-5p"),])

# cor(as.numeric(test.miRNA["A_25_P00010986",]),as.numeric(test.underexp["200755_s_at",]))

# save data
colnames(temp)<-c("pearson.coeff","p.value")
temp<-cbind(names.cor,temp)
# write.table(temp,file="./Gene Expression/gene underexp_miRNA overexp_correlations.txt",quote=F,sep="\t")

#
# correlation between gene underexpression and DNA hypermethylation 
#

# read in files
# underexp<-read.table("./Gene Expression/gene exp_batch 9_tumor/final hits_underexp_batch 9.txt",sep="\t",stringsAsFactors=F,header=T)
# DNA.hypo<-read.table("./DNA Methylation/DNA methylation_batch 9_tumor/final hits_underexp_batch 9.txt",sep="\t",stringsAsFactors=F,header=T) 

# # set up objects
# rows<-read.table("./temp2.txt",sep="\t",header=T,stringsAsFactors=F)
# rows<-rows[-19,]
# rows<-rows[,-1]
# test.DNA.hypo<-DNA.hypo[,rows[,1]]
# test.underexp<-underexp[,rows[,2]]

# int<-intersect(DNA.hypo[,1],underexp[,1])
# DNA.hypo$underexp<-DNA.hypo[,1] %in% int
# DNA.hypo<-DNA.hypo[which(DNA.hypo[,"underexp"]==T),]
# underexp$DNA.hypo<-underexp[,1] %in% int
# underexp<-underexp[which(underexp[,"DNA.hypo"]==T),]
# DNA.hypo<-cbind(DNA.hypo,is.na(DNA.hypo[,1]))
# DNA.hypo<-DNA.hypo[which(DNA.hypo[,"is.na(DNA.hypo[, 1])"]==F),]
# underexp<-cbind(underexp,is.na(underexp[,1]))
# underexp<-underexp[which(underexp[,"is.na(underexp[, 1])"]==F),]

# # perform correlation 
# which(underexp[,1]=="MEST")
# rownames(underexp)[3]

# DNA.hypo[which(DNA.hypo[,1]=="MEST"),]

# cor(as.numeric(test.DNA.hypo["cg01888566",]),as.numeric(test.underexp["202016_at",]))