###########################################################
#
# this is for merging all files in current working directory into a data matrix
#
###########################################################

DIR<-"./Gene Expression/gene exp_batch 9_normal/"
PREFIX<-"BONES_p_TCGA_Batch8_9_RNA_HT_HG-U133A_96-HTA_E"

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
ft<-gsub(".level2.data.txt","",ft)

rownames(temp)<-temp[,1]
temp<-temp[,-1]
colnames(temp)<-ft

b9.data.normal<-temp