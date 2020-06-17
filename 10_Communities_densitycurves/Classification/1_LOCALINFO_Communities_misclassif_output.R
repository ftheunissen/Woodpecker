library(MASS)
library(car)
library(ggplot2)
library(FactoMineR)
library(dplyr)

setwd("../10_Communities_densitycurves/Classification/")


rm(list=ls())#empty the workspace

# choose which database to work with depending on community to analyse

#data=read.table(file="Acoustic_alldata_postDFA.csv", header=T, sep=",",dec=".")
#data=read.table(file="Acoustic_alldata_postDFA_community1_Switzerland.csv", header=T, sep=",",dec=".")
#data=read.table(file="Acoustic_alldata_postDFA_community2_Guatemala.csv", header=T, sep=",",dec=".")
#data=read.table(file="Acoustic_alldata_postDFA_community3_Minnesota.csv", header=T, sep=",",dec=".")
#data=read.table(file="Acoustic_alldata_postDFA_community4_Malaysia.csv", header=T, sep=",",dec=".")
#data=read.table(file="Acoustic_alldata_postDFA_community5_FrGuiana.csv", header=T, sep=",",dec=".")
data=read.table(file="Acoustic_alldata_postDFA_AllCommunities.csv", header=T, sep=",",dec=".")


################### Run cross-validated DFA based on all files available (ranging from 3 to 30 files per species) 
## DFA with ALL data and CV
ndata=nrow(data)
species.names = levels(data$Species) #'Species' is the variable by which we want to discriminate drums
nbirds=length(species.names)
# There are nbirds = 92 species.

resLDA=lda(Species~PC1+PC2+PC3+PC4+PC5+PC6, prior=rep(1,nbirds)/nbirds, data=data, CV=T) # LDA computation. Important for information calculation: equal prior probabilities
tab <- table(data$Species, resLDA$class)
tab  # Count table
sum(tab[row(tab) == col(tab)]) / sum(tab) #correct classification of drums on the classification matrix's diagonal

# Probability table
prob.tab <- tab/sum(tab)

#  Conditional Probability table or Confusion Matrix using posterior probabilities
# Rows are actual species and columns the probability of being assigned to another species
# As these are conditional probabilities - the sum of each row sums to 1.
cprob.tab.post=matrix(0,nrow=nbirds,ncol=nbirds)
colnames(cprob.tab.post) <- colnames(prob.tab)
rownames(cprob.tab.post) <- rownames(prob.tab)
# Note colnames(prob.tab) == rownames(prob.tab) == species.names

for (i in 1:nbirds) {
  for (j in 1:nbirds) {
    # Select the rows that correspond to a specific species
    ind.good <- (data$Species == species.names[i])
    cprob.tab.post[i,j] <- mean(na.omit(resLDA$posterior[ind.good,j])) #na.omit to fix one line in community 2 for Col rub where posterior probabilities are not computed. Does not change anything for the other communities
  }
}

# mean classification %age
(pcc <- mean(diag(cprob.tab.post)))

# classification %age at chance level
(chance=1/nbirds)




##################################
# Calculating correct classification by pairs

### example with a pair BD (species B and species D)
### D is wrongly put N% of the time into B, compared to its own classification (D%)
### B is wrongly put NN% of teh time into D, compared to its own classification (B%)

### A way to calculate correct classification is here defined by measuring the distance between own and paired classifications:
### Correct classification within the pair = [(D% - N%) + (B% - NN%)]/2 *100 (multiplying by 100 to put it into 'percentage')
### if D always wrongly placed into B (and then B% = -1) and B always wrongly placed into D (thus D%=-1), correct classif = -100
### if D and B always well classified (thus D% and B% are both = 1), then correct classif = 100
newMat2 <- matrix(0, nrow = nbirds, ncol = nbirds)

for (X in 1:nbirds){ # First go through rows
  for (Y in 1:nbirds){ # then go through columns (first rows then columns has to be in that order because of the layout of the phyl and acoustics matrices)
    if (X >= Y){newMat2[X,Y] <- NA}
    else {newMat2[X,Y] <- mean(c(cprob.tab.post[X,X]-cprob.tab.post[X,Y], cprob.tab.post[Y,Y]-cprob.tab.post[Y,X]))*100}  
  }
}

rownames(newMat2) <- rownames(cprob.tab.post) # changed val_matrix into cprob.tab.post
colnames(newMat2) <- colnames(cprob.tab.post) # changed val_matrix into cprob.tab.post

# choose which classification matrix to write depending on community to analyse

#write.csv(newMat2, 'ClassifMatrix_fulldataset.csv')
#write.csv(newMat2, 'ClassifMatrix_Community1_Switzerland.csv')
#write.csv(newMat2, 'ClassifMatrix_Community2_Guatemala.csv')
#write.csv(newMat2, 'ClassifMatrix_Community3_Minnesota.csv')
#write.csv(newMat2, 'ClassifMatrix_Community4_Malaysia.csv')
#write.csv(newMat2, 'ClassifMatrix_Community5_FrGuiana.csv')
write.csv(newMat2, 'ClassifMatrix_AllCommunities.csv')

