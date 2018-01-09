# script for generating batch bias dendrograph

h1 <- hclust(dist(t(dat)), method="ward")
plotDendroAndColors(h1, colors=as.numeric(as.factor(dat.batch))+3,
                    family="mono", hang=-1)