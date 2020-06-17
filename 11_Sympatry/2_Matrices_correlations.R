#Change working directory to folder with your files (change accordingly)
setwd("../11_Sympatry/")
rm(list=ls())#empty the workspace

library(ggplot2)
library(gplots)
library(car)
library(RVAideMemoire)
library(lsmeans)
library(fitdistrplus)

######################
## This step plots halves (because each half has the same info) of phylo and/or acoustics and/or Classif matrices while color coding for sympatry

AcousticMatrix=read.table(file="Acousticdist_matrix.csv", header=T, sep=",", dec=".", row.names = 1)
PhylMatrix=read.table(file="Phyldist_matrix.csv", header=T, sep=",", dec=".", row.names = 1)
ClassifMatrix = read.table(file="ClassifMatrix_fulldataset.csv", header=T, sep=",", dec=".", row.names = 1)
SympatryMatrix=read.table(file="Sympatry_bypair.csv", header=T, sep=",", dec=".", row.names = 1)


#### Following section aims at comparing the half-matrices
# First, need to format matrices to real matrix
AcousticMatrix <- as.matrix(AcousticMatrix)
PhylMatrix <- as.matrix(PhylMatrix)
SympatryMatrix <- as.matrix(SympatryMatrix)
ClassifMatrix <- as.matrix(ClassifMatrix)


# Then possible to extract all values into a single dimension array (which we can then plot)
Acoustics <- as.vector(AcousticMatrix)
Phylo <- as.vector(PhylMatrix)
Sympatry <- as.vector(SympatryMatrix)
Classif <- as.vector(ClassifMatrix)

# testDF <- data.frame(Acoustics, Phylo, Classif, Sympatry) #control line

# Remove NA values (could have been done in step above but useful to control for vector arrangement before)
TrimmedAcoustics <- as.numeric(na.omit(Acoustics))
TrimmedPhylo <- as.numeric(na.omit(Phylo))
TrimmedSympatry <- as.factor(na.omit(Sympatry))
TrimmedClassif <- as.numeric(na.omit(Classif))

newdata <- data.frame(TrimmedAcoustics, TrimmedPhylo, TrimmedClassif, TrimmedSympatry)

## Correlation tests
# Overall correlation
cor.test(newdata$TrimmedAcoustics, newdata$TrimmedPhylo, method = "spearman")
cor.test(newdata$TrimmedClassif, newdata$TrimmedPhylo, method = "spearman")
cor.test(newdata$TrimmedAcoustics, newdata$TrimmedClassif, method = "spearman")




###################### Key figures ##########
#pdf("Ac_LOGphy_BySymp.pdf")
ggplot(newdata,aes(y = TrimmedAcoustics, x = log10(TrimmedPhylo), colour=TrimmedSympatry,shape=TrimmedSympatry)) +
  scale_color_manual(values=c("coral4", "darkgoldenrod1")) +
  geom_point(size=2) + geom_smooth(method="lm",formula = y~exp(x), fill=NA, lwd =1) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.ticks = element_line(size=0.5), axis.title = element_text(size=25,face = "bold"), axis.text = element_text(size=25), panel.border = element_rect(colour = "black", size=1)) +
  xlab("Phylogenetic Distance") + ylab("Acoustic Distance") + ylim(0, 15)
#dev.off()


#pdf("Ac_Classif_BySymp.pdf")
ggplot(newdata,aes(y = TrimmedClassif, x = TrimmedAcoustics, colour=TrimmedSympatry,shape=TrimmedSympatry)) +
  scale_color_manual(values=c("coral4", "darkgoldenrod1")) +
  geom_point(size=2) + geom_smooth(method="lm", fill=NA, lwd =1) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.ticks = element_line(size=0.5), axis.title = element_text(size=25,face = "bold"), axis.text = element_text(size=25), panel.border = element_rect(colour = "black", size=1)) +
  xlab("Acoustic Distance") + ylab("Classification Index") 
#dev.off()

#pdf("Classif_LogPhyl_BySymp.pdf")
ggplot(newdata,aes(y = TrimmedClassif, x = log10(TrimmedPhylo), colour=TrimmedSympatry,shape=TrimmedSympatry)) +
  scale_color_manual(values=c("coral4", "darkgoldenrod1")) +
  geom_point(size=2) + geom_smooth(method="lm", formula = y~exp(x), fill=NA, lwd =1) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.ticks = element_line(size=0.5), axis.title = element_text(size=25,face = "bold"), axis.text = element_text(size=25), panel.border = element_rect(colour = "black", size=1)) +
  xlab("Phylogenetic Distance") + ylab("Classification Index") 
#dev.off()


#### STATISTICAL TESTS
#pdf('Output_LMs_Sympatry.pdf')

lm_AcPhyl <- lm(TrimmedAcoustics~log10(TrimmedPhylo)+TrimmedSympatry+log10(TrimmedPhylo)*TrimmedSympatry,data = newdata)
summary(lm_AcPhyl)
textplot(capture.output(summary(lm_AcPhyl)))
plotresid(lm_AcPhyl)

lm_AcClassif <- lm(TrimmedClassif~TrimmedAcoustics+TrimmedSympatry+TrimmedAcoustics*TrimmedSympatry,data = newdata)
summary(lm_AcClassif)
textplot(capture.output(summary(lm_AcClassif)))
plotresid(lm_AcClassif)

lm_ClassifPhyl <- lm(TrimmedClassif~log10(TrimmedPhylo)+TrimmedSympatry+log10(TrimmedPhylo)*TrimmedSympatry,data = newdata)
summary(lm_ClassifPhyl)
textplot(capture.output(summary(lm_ClassifPhyl)))
plotresid(lm_ClassifPhyl)

#dev.off()


