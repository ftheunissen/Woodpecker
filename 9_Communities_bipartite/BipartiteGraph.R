library(igraph)

setwd("../9_Communities_bipartite/")

rm(list=ls())#empty the workspace
data=read.table(file="AcClust_Communities.csv", header=T, sep=",",dec=".")


twoWay_freqTable <- xtabs(~Community_ID+AcousticClust, data)
rownames(twoWay_freqTable) <- c("Switzerland", "Tykal", "Minnesota", "Malaysia", "Nouragues")
colnames(twoWay_freqTable) <- c("Acceleration", "Double \n Knock", "Steady \n fast", "Steady \n slow", "Steady \n Sequence", "Non-steady \n Sequence")

datFrame <- as.data.frame(twoWay_freqTable)
g <- graph.data.frame(twoWay_freqTable, directed = F)
V(g)
bipartite.mapping(g)
V(g)$label <- V(g)$name
V(g)$type <- bipartite.mapping(g)$type
is.bipartite(g)
plot(g, layout = layout.bipartite, edge.width = E(g)$Freq)


## Plots:

#pdf("communities_acoustic_bipartite.pdf")
plot(g, layout = layout.bipartite, edge.width = E(g)$Freq, main = "Use of drumming strategies \n within communities")
#dev.off()

### if want to remove all labels
#pdf("communities_acoustic_bipartite_nolabels.pdf")
plot(g, layout = layout.bipartite, edge.width = E(g)$Freq, main = "Use of drumming strategies \n within communities", vertex.label="")
#dev.off()

# for some reason, pdf() function in the graphs above messes up the frequency distribution (represented as trait width) of the strategies. 
# Thus make figure using trait width based on 'twoWay_freqTable' raw values
