# obtaining TCGA sample identities 
p1.mapping<-sapply(mapping[,1], function(x) strsplit(x, "::")[[1]][1])

# obtaining gene expression colnames 
p2.mapping<-sapply(mapping[,1], function(x) strsplit(x, "HTA_")[[1]][2])
p2.mapping<-gsub(".level2.data.txt", "", sapply(mapping[,1], function(x) strsplit(x, "HTA_")[[1]][2]))
p2.mapping<-sapply(p2.mapping, function(x) strsplit(x, " ")[[1]][1])

# match sample identities with colnames
mapping<-cbind(p1.mapping, p2.mapping)
mapping[,1]<-gsub("-", ".", mapping[,1])