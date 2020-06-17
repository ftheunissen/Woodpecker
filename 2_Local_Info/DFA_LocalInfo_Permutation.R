library(MASS)
library(car)
library(FactoMineR)

setwd("../2_Local_Info/")

rm(list=ls())#empty the workspace

# Read the z-scored data
data=read.table(file="AcousticsFINAL_Zscored.csv", header=T, sep=",",dec=".")

(nrow(data))  # Number of measurements
(ncol(data))  # Number of features considered (acoustical + metadata)

# PCA calculation 
resPCA <- PCA (data, quali.sup=1:16, ncp=22, scale.unit = T) # quali.sup= categorical supplementary variables (here species, genus...) not used in PCA calculation; ncp= nb of PC to report (initially, should be equal to the number of variables used in PCA)
resPCA$eig # get percentages of variance explained by each PC and eigenvalues (make sure data are z-scored before this)
resPCA$var$cor # get loading scores of each acoustic variable onto PCs
resPCA$ind$coord # get PC scores for each 'individual' (i.e. here: each drumming file) in the dataset

# Note that 6 componnents camputer 75% of the variance; 10 would give 90%.

PCALoadings <- resPCA$var$cor[,1:6] #based on the PCs for which eigenvalues are > 1
# write.csv(PCALoadings, "PCALoadings_2019.csv")

# adding PCs (with eigenval > 1) to dataframe
data$PC1<-resPCA$ind$coord[,1]
data$PC2<-resPCA$ind$coord[,2]
data$PC3<-resPCA$ind$coord[,3]
data$PC4<-resPCA$ind$coord[,4]
data$PC5<-resPCA$ind$coord[,5]
data$PC6<-resPCA$ind$coord[,6]



################### Run cross-validated DFA based on all files available (ranging from 3 to 30 files per species) 
## DFA with ALL data and CV
ndata=nrow(data)
species.names = levels(data$Species) #'Species' is the variable by which we want to discriminate drums
nbirds=length(species.names)
# There are nbirds = 92 species.

nPermute <- 1000
resLDAPermute <- array(0, dim=nPermute)

for (i in 1:nPermute) {
  data$Species <- sample(data$Species, size = length(data$Species), replace=FALSE)
  resLDA <- lda(Species~PC1+PC2+PC3+PC4+PC5+PC6, prior=rep(1  ,nbirds)/nbirds, data=data, CV=T) # LDA computation. Important for information calculation: equal prior probabilities
  tab <- table(data$Species, resLDA$class)

  resLDAPermute[i] <- sum(tab[row(tab) == col(tab)]) / sum(tab) #correct classification of drums on the classification matrix's diagonal
}

hist(resLDAPermute*100, xlab ='Percent Correct', main='Permutation Test')
abline(v=(1/92)*100, lty='dashed' )

