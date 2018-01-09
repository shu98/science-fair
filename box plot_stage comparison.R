test<-read.table("./Gene Expression/gene exp_clinical stages.txt",sep="\t",stringsAsFactors=F,header=T)

all.data<-read.table("./Gene Expression/gene exp_all batches.txt",sep="\t",stringsAsFactors=F,header=T)

temp<-intersect(rownames(test),colnames(all.data))

test<-test[temp,]

overexp<-read.table("./Gene Expression/gene exp_stages comparison_overexp.txt",stringsAsFactors=F,header=T,sep="\t")
overexp.data<-all.data[rownames(overexp),]
overexp.data<-overexp.data[,rownames(test)]
test.idx<-!is.na(test[,2])
test[, "Status"]<-gsub("(A|B|C)", "", test[,2])

pdf("overexp_plot.pdf")

for ( i in 1:nrow(overexp)) {
	g.name<-overexp[rownames(overexp.data)[i], 1]
	p.value<-signif(overexp[rownames(overexp.data)[i], 2], 3)
	boxplot(as.numeric(overexp.data[i, test.idx])~as.factor(test[test.idx, 3]), ylab="Log2(Expression)", las=3,
	main=paste(g.name, "(p.val=", p.value, ")"))
 }
 
dev.off()

underexp<-read.table("./Gene Expression/gene exp_stages comparison_underexp.txt",stringsAsFactors=F,header=T,sep="\t")
underexp.data<-all.data[rownames(underexp),]
underexp.data<-underexp.data[, rownames(test)]
test.idx<-!is.na(test[,2])

pdf("underexp_plot.pdf")

for ( i in 1:nrow(underexp)) {
	g.name<-underexp[rownames(underexp.data)[i], 1]
	p.value<-signif(underexp[rownames(underexp.data)[i], 2], 3)
	boxplot(as.numeric(underexp.data[i, test.idx])~as.factor(test[test.idx, 3]), ylab="Log2(Expression)", las=3,
	main=paste(g.name, "(p.val=", p.value, ")"))
 }
 
dev.off()
