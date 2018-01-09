#
# set up data 
#

# prepare DNA alteration data
dna.alt<-read.table("./Gene Expression/gene underexp DNA alterations.txt", sep="\t", stringsAsFactors=F, header=T)
rownames(dna.alt)<-dna.alt[,1]
dna.alt<-dna.alt[,-1]
dna.del<-dna.alt[,1:20]

# prepare gene underexpression data
b9.data<-read.table("./Gene Expression/gene exp_batch 9_tumor/final hits_underexp_batch 9.txt", sep="\t", stringsAsFactors=F, header=T)
b9.data<-b9.data[!duplicated(b9.data[,1]),]
b9.data<-b9.data[!is.na(b9.data[,1]),]
rownames(b9.data)<-b9.data[,1]
b9.data<-b9.data[,4:41]

# set same colnames for DNA alteration and gene underexpression data 
mapping<-read.table("./tcga samples_gene exp data_names.txt",sep="\t",stringsAsFactors=F,header=T)
underexp.del<-b9.data[,mapping[colnames(dna.del),]]
colnames(underexp.del)<-colnames(dna.del)

# remove all genes that have no sdelles with DNA alterations 
tmp.idx<-apply(dna.del, 1, function(x) length(unique(x))>=2)
underexp.del.sub<-underexp.del[tmp.idx,]
dna.del.sub<-dna.del[tmp.idx,]

# replace NA values in gene underexpression data with "WT"
# dna.del.sub<-apply(dna.del.sub, 2, function(x) ifelse(is.na(x), "WT", x))

#
# execute wilcox test
#
del.res<-NULL

for (i in 1:nrow(underexp.del.sub)){     
    
    if (length(unique(dna.del.sub[i,]))>2){
res<-kruskal.test(as.numeric(underexp.del.sub[i,])~as.factor(as.character(dna.del.sub[i,])))$p.value} else {
res<-wilcox.test(as.numeric(underexp.del.sub[i,])~as.factor(as.character(dna.del.sub[i,])))$p.value}
	
	del.res<-c(del.res,res)
	
	}

# temp<-NULL
# for (i in 1:nrow(dna.alt)){
	# x<-length(which(dna.alt[i,]=="AMP"))
	# temp<-c(temp,x)}

# temp<-NULL
# for (i in 1:nrow(dna.alt)){
	# x<-(alt.res[i]/ncol(dna.alt))
	# temp<-c(temp,x)}
