
###########################################################
#
# this is for merging all files in current working directory into a data matrix
#
###########################################################

fs<-list.files(pattern="level2.data.txt", path="./Gene Expression/gene exp_batch 9_tumor/")


temp<-read.table(paste("./Gene Expression/gene exp_batch 9_tumor/", fs[1], sep=""), sep="\t", header=T, stringsAsFactors=F, skip=1)

for(f in fs[-1]) 
{
	f1<-read.table(paste("./Gene Expression/gene exp_batch 9_tumor/", f, sep=""), sep="\t", header=T, stringsAsFactors=F, skip=1)
	stopifnot(identical(f1[,1],temp[,1]))
	temp<-cbind(temp,f1[,-1])
}

# generating column names for "temp"
ft<-substr(fs,nchar("BONES_p_TCGA_Batch8_9_RNA_HT_HG-U133A_96-HTA_E"),200)
ft<-gsub(".level2.data.txt","",ft)
# manual correction
ft[38]<-"E01_586136"
# ft[49]<-"E03_586140"
# ft[50]<-"E06_586162"
# ft[51]<-"F08_586108"

rownames(temp)<-temp[,1]
temp<-temp[,-1]
colnames(temp)<-ft

b9.data<-temp

#boxplot(b9.data, las=3)