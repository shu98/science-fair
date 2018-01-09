###########################################################
#
# this is for merging all files in current working directory into a data matrix
#
###########################################################

DIR<-"./miRNA/miRNA_batch 9_normal/"
PREFIX<-"TCGA-13-0"

fs<-list.files(pattern="level2.data.txt", path=DIR)

temp<-read.table(paste(DIR, fs[1], sep=""), sep="\t", header=T, stringsAsFactors=F, skip=1)

for(f in fs[-1]) 
{
	f1<-read.table(paste(DIR, f, sep=""), sep="\t", header=T, stringsAsFactors=F, skip=1)
	stopifnot(identical(f1[,1],temp[,1]))
	temp<-cbind(temp,f1[,-1])
}

# generating column names for "temp"
ft<-substr(fs,nchar(PREFIX),200)
ft<-gsub("-0364-07.probe.tcga_level2.data.txt","",ft)

rownames(temp)<-temp[,1]
temp<-temp[,-1]
colnames(temp)<-ft

b9.data.normal<-temp