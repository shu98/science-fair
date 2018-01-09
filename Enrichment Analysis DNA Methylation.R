# Pathway enrichment analysis for DNA methylation 

# set up data
library(annotate)
library(hgu133a2.db)
library(GSEABase)

genes.on.chip<-read.table("./DNA Methylation/gene list.txt", stringsAsFactors=F, header=T, sep="\t")
genes.on.chip<-genes.on.chip[,1]
sets <- getGmt("./Downloaded Data/c2.cp.v4.0.symbols.gmt")
genes.by.pathway <- lapply(sets, geneIds)
names(genes.by.pathway) <- names(sets)

genes.of.interest<-read.table("./DNA Methylation/DNA methylation_batch 9_tumor/final hits_underexp_batch 9.txt", stringsAsFactors=F, header=T, sep="\t")
genes.of.interest<-genes.of.interest[,1]
all.geneIDs<-unique(genes.on.chip)

# enrichment analysis
f1.hyper.enrich<-function(genes.of.interest, all.geneIDs)
 {
     pVals.by.pathway <- t(sapply(names(genes.by.pathway),
                               function(pathway) {
                                  pathway.genes <- genes.by.pathway[[pathway]]
                                  white.balls.drawn <- length(intersect(genes.of.interest, pathway.genes))
                                  white.balls.in.urn <- length(pathway.genes)
                                  total.balls.in.urn <- length(all.geneIDs)
                                  black.balls.in.urn <- total.balls.in.urn - white.balls.in.urn
                                  total.balls.drawn.from.urn <- length(genes.of.interest)
                                  #total.balls.drawn.from.urn <-length(intersect(genes.of.interest, all.geneIDs))
                                  list(
                                     pval=dhyper(white.balls.drawn, white.balls.in.urn,
                                         black.balls.in.urn, total.balls.drawn.from.urn), sig.n=white.balls.drawn, all.n=white.balls.in.urn)}))
                                        
     tmp<-apply(pVals.by.pathway, 2, as.numeric)   
     rownames(tmp)<-rownames(pVals.by.pathway)
     pVals.by.pathway<-tmp
                             
     adj.p<-p.adjust(pVals.by.pathway[, "pval"], "BH", n=sum(!is.na(pVals.by.pathway[, "pval"])))
    
     flt<-adj.p<0.05 & pVals.by.pathway[, "sig.n"] >= 1                                   
     pVals.by.pathway<-cbind(pVals.by.pathway[flt, ], adj.p[flt])[, c(4,3,2,1)]
     pVals.by.pathway[, c(1,4)]<-signif(pVals.by.pathway[,c(1,4)],3)
     colnames(pVals.by.pathway)[1]<-c("Adj.P")
     pVals.by.pathway<-pVals.by.pathway[order(pVals.by.pathway[,"Adj.P"], decreasing=F),]
     pVals.by.pathway
 }
 
enriched.pathway<-f1.hyper.enrich(genes.of.interest, all.geneIDs)

# write.table(enriched.pathway,file="./DNA Methylation/DNA methylation_batch 9_tumor/enriched pathways_overexp.txt",sep="\t",quote=F)

# source("http://www.bioconductor.org/biocLite.R")
# biocLite("hgu133a2.db")

# find overexpressed genes in each pathway 
# temp<-NULL
# for (i in 1:nrow(enriched.pathway)){
	# x<-genes.of.interest %in% genes.by.pathway[[rownames(enriched.pathway)[i]]]
	# temp<-cbind(temp,x)}
	
# colnames(temp)<-rownames(enriched.pathway)
# rownames(temp)<-genes.of.interest

# write.table(temp,file="./Gene Expression/gene exp_batch 9_tumor/enriched pathways_genes_underexp.txt",sep="\t",quote=F)