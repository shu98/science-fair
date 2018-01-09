###########################################################
#
# this is for merging all files in current working directory into a data matrix
#
###########################################################

DIR<-"./DNA Methylation/DNA methylation_batch 9_normal/"
PREFIX<-"jhu-usc.edu_OV.HumanMethylation27.1.lvl-3.TCGA-01-0"

fs<-list.files(pattern="vl-3", path=DIR)

temp<-read.table(paste(DIR, fs[1], sep=""), sep="\t", header=T, stringsAsFactors=F, skip=1)

temp<-temp[,c(1,2)]

for(f in fs[-1]) 
{
	f1<-read.table(paste(DIR, f, sep=""), sep="\t", header=T, stringsAsFactors=F, skip=1)
	stopifnot(identical(f1[,1],temp[,1]))
	temp<-cbind(temp,f1[,-c(1,3,4,5)])
}

# generating column names for "temp"
ft<-substr(fs,nchar(PREFIX),200)
ft<-gsub("-0383-05.txt","",ft)

rownames(temp)<-temp[,1]
temp<-temp[,-1]
colnames(temp)<-ft

b9.data.normal<-temp

b9.data.normal<-b9.data.normal[rownames(b9.data),]