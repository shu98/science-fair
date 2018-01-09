
###########################################################
#
# this is for merging all files in current working directory into a data matrix
#
###########################################################

fs<-list.files(pattern="level2.data.txt", path="./miRNA/miRNA_batch 9_tumor/")


temp<-read.table(paste("./miRNA/miRNA_batch 9_tumor/", fs[1], sep=""), sep="\t", header=T, stringsAsFactors=F, skip=1)

for(f in fs[-1]) 
{
	f1<-read.table(paste("./miRNA/miRNA_batch 9_tumor/", f, sep=""), sep="\t", header=T, stringsAsFactors=F, skip=1)
	stopifnot(identical(f1[,1],temp[,1]))
	temp<-cbind(temp,f1[,-1])
}

# generating column names for "temp"
ft<-substr(fs,nchar("TCGA-13-0"),200)
ft<-gsub("-0364-07.probe.tcga_level2.data.txt","",ft)
# manual correction
ft[27]<-"0791-02A-01T"
# ft[32]<-"1707-02A-01T"
# ft[34]<-"1710-02A-01T"
# ft[35]<-"1711-01A-01T"
# ft[41]<-"1906-01A-01T"
# ft[43]<-"1770-02A-01T"

rownames(temp)<-temp[,1]
temp<-temp[,-1]
colnames(temp)<-ft

b9.data<-temp

#boxplot(b9.data, las=3)