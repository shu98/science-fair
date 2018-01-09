# heat maps 

f.heatmap <- function(data,
                      cl=rep("1",ncol(data)),
                      h1="yes",
                      h2="yes",
                      index=NULL,
                      n.genes=NULL,
                      probeid=NULL,
                      plot=T,
                      colsidecolors1=matrix(rep("0",2*ncol(data)), ncol=2),
                      col1=greenred.colors(128),
                      margins=c(6,15),
                      affy=F,
                      #scale1=F,
                      platform="hgu133plus2",
                      genenames=NULL,
                      scale="row",
                      ...
                      )
{
    require(geneplotter)
    require(limma)
    require(heatmap.plus)
    if(!affy) col1 <- dChip.colors(128)
    if(!is.null(index)) {
        data <- data[,index]
        cl <- cl[index]
    } else {
        index <- rep(T, ncol(data))
    }
    if(!is.null(n.genes)) {
        filter.tmp <- eBayes(lmFit(data, design=cbind(1, cl=="sens")))$p.value[,2]
        index1 <- rank(filter.tmp)<=n.genes
    } else {
        if(is.null(probeid)) {
            index1 <- rep(T, nrow(data))
        } else {
            index1 <- is.element(row.names(data), probeid)
        }
    }
    if(is.null(genenames)) {
        if(affy) {
            #genenames <- unlist(mget(row.names(data)[index1], hgu133plus2GENENAME))
            eval(parse(text=paste(c("genenames <- unlist(AnnotationDbi::mget(row.names(data)[index1], ", platform, "GENENAME))"), collapse="")))
            genenames <- apply(cbind(names(genenames), genenames), 1, paste, collapse="//")
        } else {
            genenames <- row.names(data[index1,])
        }
    }
    if(is.null(colsidecolors1)) {
        colsidecolors1 <- apply(cbind(subtype1=sampleinfo$subtype1[index], #ER=sampleinfo$ER[index],
                         "res/sens"=cl), 2, function(x) as.character(as.numeric(as.factor(x)) ))
    } else {
        colsidecolors1 <- colsidecolors1[index,]
    }
#    if(scale1) {
#        dat.tmp <- t(apply(data, 1, scale))
#        dimnames(dat.tmp) <- dimnames(data)
#        data <- dat.tmp
#    }
    if(!is.na(h1)) {
        if(class(h1)=="hclust") {
            h1 <- as.dendrogram(h1)
        } else {
            if(affy) {
                h1 <- as.dendrogram(hclust(dist(t(data[index1,])), method="ward"))
            } else {
                h1 <- as.dendrogram(hclust(as.dist(1-cor(data[index1,])), 
                #use="pairwise.complete.obs")), 
                method="ward"))
            }
        }
    }
    if(!is.na(h2)) {
        if(class(h2)=="hclust") {
            h2 <- as.dendrogram(h2)
        } else {
            h2 <- as.dendrogram(hclust(as.dist(1-cor(t(data[index1,]))), 
            #use="pairwise.complete.obs")), 
            method="ward"))
        }
    }
    if(plot) heatmap.plus(data[index1,], margins=margins, ColSideColors=colsidecolors1, scale=scale, col=col1, labRow=genenames,
     Rowv=h2, Colv=h1, ... )
    if(!is.null(n.genes)) return(row.names(data)[index1])
}

# source("http://www.bioconductor.org/biocLite.R")
# biocLite("heatmap.plus")
# library(heatmap.plus)
# source("./Scripts/f.heatmap.R")
# biocLite("limma")
# biocLite("geneplotter")
# library(geneplotter)
# f.heatmap(b11[, -c(1:3)])
# f.heatmap(data.matrix(b11[, 4:49]))
# intersect(rownames(b12), intersect(rownames(b9), rownames(b11)))->com.hits
# f.heatmap(data.matrix(b12[com.hits, 4:58]))
# pdf(file="heatmaps.pdf")
# f.heatmap(data.matrix(b12[com.hits, 4:58]), main="b12, com.hits", genenames="")
# dev.off()