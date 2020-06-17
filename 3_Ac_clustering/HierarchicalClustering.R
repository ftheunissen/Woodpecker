library(NbClust)
library(car)
library(rgl)
library(dendextend)


# for acoustic groups:
setwd("../3_Ac_clustering/")
rm(list=ls())#empty the workspace
data=read.table(file="Acoustic_species_postDFA.csv", header=T, sep=",",dec=".") # Acoustic data already Z-scored

##### prepare dataset so that only variables are included. Note: file is 'postPCA' hence already centered/scaled
rownames(data) <- as.character(unlist(data$Species))

##### hierarchical clustering using NbClusters --> finds the best number of clusters according to the majority rule
clusterAnalysisWard <- NbClust(data[16:37], method = 'ward.D2', distance = "euclidean") # based on the 22 acoustic variables (Z-scored) for each median specific drum
acousticGroups_clusters <- clusterAnalysisWard$Best.partition

outputdata <- as.data.frame(acousticGroups_clusters)
# write.csv(outputdata, "acousticGroups_clusters_NbClust.csv")


##### hierarchical clustering using hclust and Ward.D2 method (exactly like with NbClusters) --> easier for plotting
hc = hclust(dist(data[16:37]), method="ward.D2")

hc2 <- as.dendrogram(hc, hang=-1) # if wanna place labels at same level in plot, add "hang = -1"
hc2 <- set(hc2, "labels_cex", 0.35)


##### Cutting tree to get 6 clusters (base on NBclust result; achieved by cutting at height = 11 based on visual inspection)
clusterList <- cutree(hc, k=6, h=11)
clusterList <- as.data.frame(clusterList) # check line: similar to 'outputdata2' here above
# write.csv(clusterList, 'acousticGroups_clusters_hclust.csv')

# set dendrogram color labels to match those of the ordered cluster list
# first reorder cluterList according to dendrogram
orderedClusterList <- clusterList[match(labels(hc2),rownames(clusterList)),]
names(orderedClusterList) <- labels(hc2)
orderedClusterList <- as.data.frame(orderedClusterList)
labels_colors(hc2) <- orderedClusterList[,1]

##### plot dendrogram of the clustering
## If want to change legend orientation, use text() function instead of mtext():
pdf("Hierachical_Clustering.pdf")
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 1, 0, 0)) # can be used to defined larger margins for PDF file

plot(hc2, cex=0.3, horiz = T, cex.axis=1.3) # plot tree without colored branches
title("Hierarchical clustering - Ward's method", line = 2)
#abline(v=11, col='red') # if horizontal: 'h' stand for height: value at which want to draw cutting line; if vertical, replace 'h' by 'v'

text(27,90,labels = 'AcCluster 1: acceleration',col = 1,cex = 0.9, xpd = NA, srt = 0)
text(27,86,labels = 'AcCluster 2: double knock',col = 2,cex = 0.9, xpd = NA, srt = 0)
text(27,82,labels = 'AcCluster 3: steady fast',col = 3,cex = 0.9, xpd = NA, srt = 0)
text(27,78,labels = 'AcCluster 4: steady slow',col = 4,cex = 0.9, xpd = NA, srt = 0)
text(27,74,labels = 'AcCluster 5: steady sequences',col = 5,cex = 0.9, xpd = NA, srt = 0)
text(27,70,labels = 'AcCluster 6: non-steady sequences',col = 6,cex = 0.9, xpd = NA, srt = 0)


hc2_colored=color_branches(hc2,k=6, col = c(5,6,2,3,1,4)) #manually adjust color order to match branches and tips
plot(hc2_colored, cex=0.3, horiz=T, cex.axis=1.3)
title("Hierarchical clustering - Ward's method", line = 2)
#abline(v=11, col='red') # if horizontal: 'h' stand for height: value at which want to draw cutting line; if vertical, replace 'h' by 'v'
text(27,90,labels = 'AcCluster 1: acceleration',col = 1,cex = 0.9, xpd = NA, srt = 0)
text(27,86,labels = 'AcCluster 2: double knock',col = 2,cex = 0.9, xpd = NA, srt = 0)
text(27,82,labels = 'AcCluster 3: steady fast',col = 3,cex = 0.9, xpd = NA, srt = 0)
text(27,78,labels = 'AcCluster 4: steady slow',col = 4,cex = 0.9, xpd = NA, srt = 0)
text(27,74,labels = 'AcCluster 5: steady sequences',col = 5,cex = 0.9, xpd = NA, srt = 0)
text(27,70,labels = 'AcCluster 6: non-steady sequences',col = 6,cex = 0.9, xpd = NA, srt = 0)

dev.off()


# now plotting data in a 3D space defined by the first 3 components of the acoustic PCA
acClusters <- as.factor(data$AcousticClust)
PC1 <- data$median_PC1
PC2 <- data$median_PC2
PC3 <- data$median_PC3

scatter3d(PC1, PC2, PC3, groups = acClusters, surface=FALSE, grid = FALSE, ellipsoid = TRUE, level=0.5,
          axis.col = c("black", "black", "black"), sphere.size = 0.1, color=2,
          surface.col=c("black", "red", "green", "blue", "cyan", "purple"))

# Then rotate around the graph created in XQuartz, and once happy with it (angle where groups are most separated):
rgl.snapshot("3DPlot_PCs.png")

# or, if want 2D plot:
#pdf("PC1_PC2.pdf")
ggplot(data, aes(PC1,PC2,color = acClusters)) +
  scale_color_manual(values=c(1,2,3,4,5,6))+
  theme(panel.grid.minor = element_blank(),
         panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x=element_text(angle=0, size=17, hjust=0.5),
        axis.text.y=element_text(angle=90, size=17, hjust=0.5),
        axis.title.y = element_text(colour="grey20",size=22,angle=90,hjust=.5,face="plain"),
        axis.title.x = element_text(colour="grey20",size=22,angle=0,hjust=.5,face="plain"),
        legend.position = "none")+
  geom_point() +
  stat_ellipse(level = 0.90)
#dev.off()


# or using DFs instead of PCs
LD1 <- data$median_LD1
LD2 <- data$median_LD2
LD3 <- data$median_LD3

scatter3d(LD1, LD2, LD3, groups = acClusters, surface=FALSE, grid = FALSE, ellipsoid = TRUE, level=0.5,
          axis.col = c("black", "black", "black"), sphere.size = 0.1, color=2,
          surface.col=c("black", "red", "green", "blue", "cyan", "purple"))

# Then rotate around the graph created in XQuartz, and once happy with it (angle where groups are most separated):
rgl.snapshot("3DPlot_LDs.png")

# Correlation matrix
library(corrplot)
source("http://www.sthda.com/upload/rquery_cormat.r")
mydata <- data[16:37]
colnames(mydata) <- 1:22
rquery.cormat(mydata)



