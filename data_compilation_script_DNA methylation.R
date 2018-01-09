###########################################################
#
# this is for merging all files in current working directory into a data matrix
#
###########################################################

fs<-list.files(pattern="vl-3", path="./DNA Methylation/dna methylation_batch 9_tumor/")


temp<-read.table(paste("./DNA Methylation/dna methylation_batch 9_tumor/", fs[1], sep=""), sep="\t", header=T, stringsAsFactors=F, skip=1)

temp<-temp[,c(1,2)]

for(f in fs[-1]) 
{
	f1<-read.table(paste("./DNA Methylation/DNA methylation_batch 9_tumor/", f, sep=""), sep="\t", header=T, stringsAsFactors=F, skip=1)
	stopifnot(identical(f1[,1],temp[,1]))
	temp<-cbind(temp,f1[,-c(1,3,4,5)])
}

# generating column names for "temp"
ft<-substr(fs,nchar("jhu-usc.edu_OV.HumanMethylation27.1.lvl-3.TCGA-09-0"),200)
ft<-gsub("-0359-05.txt","",ft)
# manual correction
ft[39]<-"0791-02A-01D"
# ft[2]<-"1707-02A-01D"
# ft[3]<-"1710-02A-01D"
# ft[4]<-"1770-02A-01D"

rownames(temp)<-temp[,1]
temp<-temp[,-1]
colnames(temp)<-ft

b9.data<-temp

#b9.data<-b9.data[apply(!is.na(data.matrix(b9.data)), 1, sum)>10, ]

#boxplot(b9.data, las=3)