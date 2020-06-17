library(MASS)
library(car)
library(ggplot2)
library(FactoMineR)
library(dplyr)

setwd("../10_Communities_densitycurves/Classification/")
rm(list=ls())#empty the workspace


ClassifMatrix=read.table(file="ClassifMatrix_fulldataset.csv", header=T, sep=",", dec=".", row.names = 1) #overall (92 species)

## new values for communities (as calculated from main analysis of local info only based on species within a community)
ClassifMatrix_Com1=read.table(file="ClassifMatrix_Community1_Switzerland.csv", header=T, sep="," ,dec=".", row.names = 1)
ClassifMatrix_Com2=read.table(file="ClassifMatrix_Community2_Guatemala.csv", header=T, sep="," ,dec=".", row.names = 1)
ClassifMatrix_Com3=read.table(file="ClassifMatrix_Community3_Minnesota.csv", header=T, sep="," ,dec=".", row.names = 1)
ClassifMatrix_Com4=read.table(file="ClassifMatrix_Community4_Malaysia.csv", header=T, sep="," ,dec=".", row.names = 1)
ClassifMatrix_Com5=read.table(file="ClassifMatrix_Community5_FrGuiana.csv", header=T, sep="," ,dec=".", row.names = 1)


######### To plot densities, cannot used a dataframe format
# First, need to format from dataframe to matrix and scale it for visual representation between 0 and 100 (raw data between -100 and 100)
ClassifMatrix <- (as.matrix(ClassifMatrix))
ClassifMatrix_Com1 <- (as.matrix(ClassifMatrix_Com1))
ClassifMatrix_Com2 <- (as.matrix(ClassifMatrix_Com2))
ClassifMatrix_Com3 <- (as.matrix(ClassifMatrix_Com3))
ClassifMatrix_Com4 <- (as.matrix(ClassifMatrix_Com4))
ClassifMatrix_Com5 <- (as.matrix(ClassifMatrix_Com5))

# Then possible to extract all values into a single dimension array (which we can then plot)
ClassifGlobal <- na.omit(as.vector(ClassifMatrix))
min(ClassifGlobal) # checking for plot x-axis
max(ClassifGlobal) # checking for plot x-axis

ClassifCom1 <- na.omit(as.vector(ClassifMatrix_Com1))
min(ClassifCom1) # checking for plot x-axis
max(ClassifCom1) # checking for plot x-axis

ClassifCom2 <- na.omit(as.vector(ClassifMatrix_Com2))
min(ClassifCom2) # checking for plot x-axis
max(ClassifCom2) # checking for plot x-axis

ClassifCom3 <- na.omit(as.vector(ClassifMatrix_Com3))
min(ClassifCom3) # checking for plot x-axis
max(ClassifCom3) # checking for plot x-axis

ClassifCom4 <- na.omit(as.vector(ClassifMatrix_Com4))
min(ClassifCom4) # checking for plot x-axis
max(ClassifCom4) # checking for plot x-axis

ClassifCom5 <- na.omit(as.vector(ClassifMatrix_Com5))
min(ClassifCom5) # checking for plot x-axis
max(ClassifCom5) # checking for plot x-axis


# Putting all communities together --> keep "within community Classification values" instead of running with all 33 community species at the same time oterwise wrong "average within community Classification"
ClassifAllCom <- c(ClassifCom1,ClassifCom2,ClassifCom3,ClassifCom4,ClassifCom5)

# Then compile the kernel density estimation for each dataframe (whole or communities)
densityGlobal <- density(ClassifGlobal)
densityCom1 <- density(ClassifCom1)
densityCom2 <- density(ClassifCom2)
densityCom3 <- density(ClassifCom3)
densityCom4 <- density(ClassifCom4)
densityCom5 <- density(ClassifCom5)

densityAllCom <- density(ClassifAllCom)

# Plot first as a propability density
#pdf("DensityPlot_Classification.pdf")
par(mar=c(5,5,5,5))
plot(densityGlobal$x, densityGlobal$y, xlim=c(-50,100), col="black", type='l', xlab= "Classification Index", ylab= 'Probability Density', lwd=4, cex.axis= 1.8, cex.lab=2.2)
lines(densityCom1$x, densityCom1$y,type= 'l',col="darkorange", lwd=3)
lines(densityCom2$x, densityCom2$y,type= 'l',col="bisque1", lwd=3)     
lines(densityCom3$x, densityCom3$y,type= 'l',col="darkred", lwd=3)     
lines(densityCom4$x, densityCom4$y,type= 'l',col="darkgoldenrod1", lwd=3)     
lines(densityCom5$x, densityCom5$y,type= 'l',col="darkgoldenrod4", lwd=3)  
lines(densityGlobal$x, densityGlobal$y,type='l', col="black", lwd=5)


legend("topright", c("Global","Switzerland","Guatemala", "Minnesota", "Malaysia", "Fr.Guiana"), cex=1, 
       col=c("black","darkorange", "bisque1", "darkred", "darkgoldenrod1","darkgoldenrod4"), lwd = 3, lty=1, bty = "n")

#dev.off()

# Plot next as a cumulative probability distribution
#pdf("CumProb_Classification.pdf")
par(mar=c(5,5,5,5))
cum <- cumsum(densityGlobal$y*(densityGlobal$x[2]-densityGlobal$x[1]))
plot(densityGlobal$x, cum, xlim=c(-25,100), ylim=c(0,1.0), col="black", type='l', xlab= "Classification Index", ylab= 'Cumulative Probability Density', lwd=4, cex.axis= 1.8, cex.lab=2.2)

cum <- cumsum(densityCom1$y*(densityCom1$x[2]-densityCom1$x[1]))
lines(densityCom1$x, cum,type= 'l',col="darkorange", lwd=3)

cum <- cumsum(densityCom2$y*(densityCom2$x[2]-densityCom2$x[1]))
lines(densityCom2$x, cum,type= 'l',col="bisque1", lwd=3)  

cum <- cumsum(densityCom3$y*(densityCom3$x[2]-densityCom3$x[1]))
lines(densityCom3$x, cum,type= 'l',col="darkred", lwd=3) 

cum <- cumsum(densityCom4$y*(densityCom4$x[2]-densityCom4$x[1]))
lines(densityCom4$x, cum,type= 'l',col="darkgoldenrod1", lwd=3)     
cum <- cumsum(densityCom5$y*(densityCom5$x[2]-densityCom5$x[1]))
lines(densityCom5$x, cum,type= 'l',col="darkgoldenrod4", lwd=3)  
legend("topleft", c("Global","Switzerland","Guatemala", "Minnesota", "Malaysia", "Fr.Guiana"), cex=0.8, 
       col=c("black","darkorange", "bisque1", "darkred", "darkgoldenrod1","darkgoldenrod4"), lwd = 3, lty=1, bty = "n")

#dev.off()

wilcoxTest1 <- wilcox.test(x = ClassifGlobal, y = ClassifCom1)
wilcoxTest1
qnorm(wilcoxTest1$p.value/2)
wilcoxTest2 <- wilcox.test(x = ClassifGlobal, y = ClassifCom2)
wilcoxTest2
qnorm(wilcoxTest2$p.value/2)
wilcoxTest3 <- wilcox.test(x = ClassifGlobal, y = ClassifCom3)
wilcoxTest3
qnorm(wilcoxTest3$p.value/2)
wilcoxTest4 <- wilcox.test(x = ClassifGlobal, y = ClassifCom4)
wilcoxTest4
qnorm(wilcoxTest4$p.value/2)
wilcoxTest5 <- wilcox.test(x = ClassifGlobal, y = ClassifCom5)
wilcoxTest5
qnorm(wilcoxTest5$p.value/2)

p.adjust(as.vector(c(wilcoxTest1$p.value,wilcoxTest2$p.value,wilcoxTest3$p.value,wilcoxTest4$p.value,wilcoxTest5$p.value)),method = 'bonferroni', n=5) 



######################################
# Above was the analysis to plot the probability densities from real communities. Now using simulated communities for comparison:

data=read.table(file="Acoustic_alldata_postDFA.csv", header=T, sep=",",dec=".")

Data_AC_Strat <- droplevels(data[data$AcousticClust==1,])
Data_DK_Strat <- droplevels(data[data$AcousticClust==2,])
Data_SF_Strat <- droplevels(data[data$AcousticClust==3,])
Data_SS_Strat <- droplevels(data[data$AcousticClust==4,])
Data_RS_Strat <- droplevels(data[data$AcousticClust==5,])
Data_IS_Strat <- droplevels(data[data$AcousticClust==6,])

# Perform this simulation nBoot times
nBoot <- 100

# Look at range to get an idea
xmin <- min(c(min(densityCom1$x),min(densityCom2$x),min(densityCom3$x),
          min(densityCom4$x),min(densityCom5$x)))
xmax <- max(c(max(densityCom1$x),max(densityCom2$x),max(densityCom3$x),
              max(densityCom4$x),max(densityCom5$x)))

# Need power of 2 for density function
xrange <- -20:107
nx <- length(xrange)

# Make space to store densities
RealDensity <- matrix(nrow=nBoot, ncol=nx)
RandomDensity <- matrix(nrow=nBoot, ncol=nx)
SameDensity <- matrix(nrow=nBoot, ncol=nx)
DiffDensity <- matrix(nrow=nBoot, ncol=nx)

for (iboot in 1:nBoot) {

### Need to create 6 communities (one per drum type) of 5 species where, for each community, the 5 species have the same drum type 
  
# community RS
RS_species_index <- sample(length(levels(Data_RS_Strat$Species)),5, replace = F)
RS_species <- droplevels(Data_RS_Strat[Data_RS_Strat$Species %in% levels(Data_RS_Strat$Species)[RS_species_index],])

# community IS
IS_species_index <- sample(length(levels(Data_IS_Strat$Species)),5, replace = F)
IS_species <- droplevels(Data_IS_Strat[Data_IS_Strat$Species %in% levels(Data_IS_Strat$Species)[IS_species_index],])

# community AC
AC_species_index <- sample(length(levels(Data_AC_Strat$Species)),5, replace = F)
AC_species <- droplevels(Data_AC_Strat[Data_AC_Strat$Species %in% levels(Data_AC_Strat$Species)[AC_species_index],])

# community DK
DK_species_index <- sample(length(levels(Data_DK_Strat$Species)),5, replace = F)
DK_species <- droplevels(Data_DK_Strat[Data_DK_Strat$Species %in% levels(Data_DK_Strat$Species)[DK_species_index],])

# community SF
SF_species_index <- sample(length(levels(Data_SF_Strat$Species)),5, replace = F)
SF_species <- droplevels(Data_SF_Strat[Data_SF_Strat$Species %in% levels(Data_SF_Strat$Species)[SF_species_index],])

# community SS
SS_species_index <- sample(length(levels(Data_SS_Strat$Species)),5, replace = F)
SS_species <- droplevels(Data_SS_Strat[Data_SS_Strat$Species %in% levels(Data_SS_Strat$Species)[SS_species_index],])


### In addition, creating 4 communities with 6 species where, for each community, the 6 species 
### have the same drum type (IS and RS cannot be used as each of these drum types is used only by 5 species)
# community AC_2
AC_species_index2 <- sample(length(levels(Data_AC_Strat$Species)),6, replace = F)
AC_species2 <- droplevels(Data_AC_Strat[Data_AC_Strat$Species %in% levels(Data_AC_Strat$Species)[AC_species_index2],])

# community DK_2
DK_species_index2 <- sample(length(levels(Data_DK_Strat$Species)),6, replace = F)
DK_species2 <- droplevels(Data_DK_Strat[Data_DK_Strat$Species %in% levels(Data_DK_Strat$Species)[DK_species_index2],])

# community SF_2
SF_species_index2 <- sample(length(levels(Data_SF_Strat$Species)),6, replace = F)
SF_species2 <- droplevels(Data_SF_Strat[Data_SF_Strat$Species %in% levels(Data_SF_Strat$Species)[SF_species_index2],])

# community SS_2
SS_species_index2 <- sample(length(levels(Data_SS_Strat$Species)),6, replace = F)
SS_species2 <- droplevels(Data_SS_Strat[Data_SS_Strat$Species %in% levels(Data_SS_Strat$Species)[SS_species_index2],])


communityList <- list(RS_species,IS_species,AC_species,DK_species,SF_species,SS_species,AC_species2,DK_species2,SF_species2,SS_species2)
# communityList <- list(RS_species) #checkup line to run one community at a time


# prepare a vector to storage the correct classification index values
sameStrat_Comm = c()

iteration <- '0' #checkupline
for (community in communityList){
  #print(iteration) #checkup line

  ################### For each community, run cross-validated DFA based on all files available (ranging from 3 to 30 files per species) 
  ## DFA with random community data and CV
  ndata=nrow(community)
  species.names = levels(community$Species) #'Species' is the variable by which we want to discriminate drums
  nbirds=length(species.names)
  # There are nbirds = 92 species.
  
  resLDA=lda(Species~PC1+PC2+PC3+PC4+PC5+PC6, prior=rep(1,nbirds)/nbirds, data=community, CV=T) # LDA computation. Important for information calculation: equal prior probabilities
  tab <- table(community$Species, resLDA$class)
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
      ind.good <- (community$Species == species.names[i])
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
    for (Y in 1:nbirds){ # then go through columns
      if (X >= Y){newMat2[X,Y] <- NA}
      else {newMat2[X,Y] <- mean(c(cprob.tab.post[X,X]-cprob.tab.post[X,Y], cprob.tab.post[Y,Y]-cprob.tab.post[Y,X]))*100}
    }
  }
  
  rownames(newMat2) <- rownames(cprob.tab.post) # changed val_matrix into cprob.tab.post
  colnames(newMat2) <- colnames(cprob.tab.post) # changed val_matrix into cprob.tab.post
  
  # print(newMat2) #checkupline to visualize the correct classification index matrices that are created
  
  newMat2asVector <- na.omit(as.vector(newMat2))
  sameStrat_Comm <- c(sameStrat_Comm, newMat2asVector) # vector of 120 values (6 times 10 values, i.e. 10 values per drum type ; + 4 times 15 values)

}

# Now plotting density curves based on CCI (correct classification index = the values stored in the vector created prior to the for loop)
CCI_sameStrat_Comm <- density(sameStrat_Comm, from=-20, to=107, n=128)
SameDensity[iboot,] <- CCI_sameStrat_Comm$y


# write.csv(sameStrat_Comm,'sameStrat_CommVector.csv')


##############################
### Need to create 6 communities of 5 species and 4 communities of 6 species where each species has a different drum type
data=read.table(file="Acoustic_alldata_postDFA.csv", header=T, sep=",",dec=".")

Data_AC_Strat <- droplevels(data[data$AcousticClust==1,])
Data_DK_Strat <- droplevels(data[data$AcousticClust==2,])
Data_SF_Strat <- droplevels(data[data$AcousticClust==3,])
Data_SS_Strat <- droplevels(data[data$AcousticClust==4,])
Data_RS_Strat <- droplevels(data[data$AcousticClust==5,])
Data_IS_Strat <- droplevels(data[data$AcousticClust==6,])


# prepare a vector to storage the correct classification index values
diffStrat_Comm = c()


# create 2 for loops: one to create 6 communities of 5 species, another to create 4 communities of 6 species 
for (counter in 1:6){
  # community RS
  RS_singlespecies_index <- sample(length(levels(Data_RS_Strat$Species)),1, replace = F)
  RS_singlespecies <- droplevels(Data_RS_Strat[Data_RS_Strat$Species %in% levels(Data_RS_Strat$Species)[RS_singlespecies_index],])

  # community IS
  IS_singlespecies_index <- sample(length(levels(Data_IS_Strat$Species)),1, replace = F)
  IS_singlespecies <- droplevels(Data_IS_Strat[Data_IS_Strat$Species %in% levels(Data_IS_Strat$Species)[IS_singlespecies_index],])

  # community AC
  AC_singlespecies_index <- sample(length(levels(Data_AC_Strat$Species)),1, replace = F)
  AC_singlespecies <- droplevels(Data_AC_Strat[Data_AC_Strat$Species %in% levels(Data_AC_Strat$Species)[AC_singlespecies_index],])

  # community DK
  DK_singlespecies_index <- sample(length(levels(Data_DK_Strat$Species)),1, replace = F)
  DK_singlespecies <- droplevels(Data_DK_Strat[Data_DK_Strat$Species %in% levels(Data_DK_Strat$Species)[DK_singlespecies_index],])

  # community SF
  SF_singlespecies_index <- sample(length(levels(Data_SF_Strat$Species)),1, replace = F)
  SF_singlespecies <- droplevels(Data_SF_Strat[Data_SF_Strat$Species %in% levels(Data_SF_Strat$Species)[SF_singlespecies_index],])

  # community SS
  SS_singlespecies_index <- sample(length(levels(Data_SS_Strat$Species)),1, replace = F)
  SS_singlespecies <- droplevels(Data_SS_Strat[Data_SS_Strat$Species %in% levels(Data_SS_Strat$Species)[SS_singlespecies_index],])


  speciesList <- list(RS_singlespecies, IS_singlespecies, AC_singlespecies, DK_singlespecies, SF_singlespecies, SS_singlespecies)
  FivespeciesListsample <- sample(speciesList,5, replace = F)
  FivespeciesList <- droplevels(do.call(rbind.data.frame, FivespeciesListsample)) # At each iteration of the counter, this is the community with 5 species that we use


  ################### For each community, run cross-validated DFA based on all files available (ranging from 3 to 30 files per species) 
  ## DFA with random community data and CV
  ndata=nrow(FivespeciesList)
  species.names = levels(FivespeciesList$Species) #'Species' is the variable by which we want to discriminate drums
  nbirds=length(species.names)
  # There are nbirds = 92 species.
  
  resLDA=lda(Species~PC1+PC2+PC3+PC4+PC5+PC6, prior=rep(1,nbirds)/nbirds, data=FivespeciesList, CV=T) # LDA computation. Important for information calculation: equal prior probabilities
  tab <- table(FivespeciesList$Species, resLDA$class)
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
      ind.good <- (FivespeciesList$Species == species.names[i])
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
    for (Y in 1:nbirds){ # then go through columns 
      if (X >= Y){newMat2[X,Y] <- NA}
      else {newMat2[X,Y] <- mean(c(cprob.tab.post[X,X]-cprob.tab.post[X,Y], cprob.tab.post[Y,Y]-cprob.tab.post[Y,X]))*100}
    }
  }
  
  rownames(newMat2) <- rownames(cprob.tab.post) # changed val_matrix into cprob.tab.post
  colnames(newMat2) <- colnames(cprob.tab.post) # changed val_matrix into cprob.tab.post
  
  #print(newMat2) #checkupline to visualize the correct classification index matrices that are created
  
  newMat2asVector <- na.omit(as.vector(newMat2))
  diffStrat_Comm <- c(diffStrat_Comm, newMat2asVector) # vector of 60 values (6 times 10 values, i.e. 10 values per community)
  
}

for (counter in 1:4){
  # community RS
  RS_singlespecies_index <- sample(length(levels(Data_RS_Strat$Species)),1, replace = F)
  RS_singlespecies <- droplevels(Data_RS_Strat[Data_RS_Strat$Species %in% levels(Data_RS_Strat$Species)[RS_singlespecies_index],])
  
  # community IS
  IS_singlespecies_index <- sample(length(levels(Data_IS_Strat$Species)),1, replace = F)
  IS_singlespecies <- droplevels(Data_IS_Strat[Data_IS_Strat$Species %in% levels(Data_IS_Strat$Species)[IS_singlespecies_index],])
  
  # community AC
  AC_singlespecies_index <- sample(length(levels(Data_AC_Strat$Species)),1, replace = F)
  AC_singlespecies <- droplevels(Data_AC_Strat[Data_AC_Strat$Species %in% levels(Data_AC_Strat$Species)[AC_singlespecies_index],])
  
  # community DK
  DK_singlespecies_index <- sample(length(levels(Data_DK_Strat$Species)),1, replace = F)
  DK_singlespecies <- droplevels(Data_DK_Strat[Data_DK_Strat$Species %in% levels(Data_DK_Strat$Species)[DK_singlespecies_index],])
  
  # community SF
  SF_singlespecies_index <- sample(length(levels(Data_SF_Strat$Species)),1, replace = F)
  SF_singlespecies <- droplevels(Data_SF_Strat[Data_SF_Strat$Species %in% levels(Data_SF_Strat$Species)[SF_singlespecies_index],])
  
  # community SS
  SS_singlespecies_index <- sample(length(levels(Data_SS_Strat$Species)),1, replace = F)
  SS_singlespecies <- droplevels(Data_SS_Strat[Data_SS_Strat$Species %in% levels(Data_SS_Strat$Species)[SS_singlespecies_index],])
  
  
  SixspeciesList <- droplevels(as.data.frame(rbind(RS_singlespecies, IS_singlespecies, AC_singlespecies, DK_singlespecies, SF_singlespecies, SS_singlespecies))) # At each iteration of the counter, this is the community with 6 species that we use
  
  
  ################### For each community, run cross-validated DFA based on all files available (ranging from 3 to 30 files per species) 
  ## DFA with random community data and CV
  ndata=nrow(SixspeciesList)
  species.names = levels(SixspeciesList$Species) #'Species' is the variable by which we want to discriminate drums
  nbirds=length(species.names)
  # There are nbirds = 92 species.
  
  resLDA=lda(Species~PC1+PC2+PC3+PC4+PC5+PC6, prior=rep(1,nbirds)/nbirds, data=SixspeciesList, CV=T) # LDA computation. Important for information calculation: equal prior probabilities
  tab <- table(SixspeciesList$Species, resLDA$class)
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
      ind.good <- (SixspeciesList$Species == species.names[i])
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
    for (Y in 1:nbirds){ # then go through columns 
      if (X >= Y){newMat2[X,Y] <- NA}
      else {newMat2[X,Y] <- mean(c(cprob.tab.post[X,X]-cprob.tab.post[X,Y], cprob.tab.post[Y,Y]-cprob.tab.post[Y,X]))*100}
    }
  }
  
  rownames(newMat2) <- rownames(cprob.tab.post) # changed val_matrix into cprob.tab.post
  colnames(newMat2) <- colnames(cprob.tab.post) # changed val_matrix into cprob.tab.post
  
  #print(newMat2) #checkupline to visualize the correct classification index matrices that are created
  
  newMat2asVector <- na.omit(as.vector(newMat2)) # vector of 60 values (4 times 15 values, i.e. 15 values per community)
  diffStrat_Comm <- c(diffStrat_Comm, newMat2asVector) # vector of 120 values (6 times 10 values + 4 times 15 values, i.e. 15 values per community)
}

# Now plotting density curves based on CCI (correct classification index = the values stored in the vector created prior to the for loop)
CCI_diffStrat_Comm <- density(diffStrat_Comm, from=-20, to=107, n=128)
DiffDensity[iboot,] <- CCI_diffStrat_Comm$y

#par(mar=c(5,5,5,5))
#plot(CCI_diffStrat_Comm$x, CCI_diffStrat_Comm$y, xlim=c(-50,100), col="blue", type='l', xlab= "Correct classification Index", ylab= 'Probability density', lwd=4, cex.axis= 1.8, cex.lab=2.2)


# write.csv(diffStrat_Comm,'diffStrat_CommVector.csv')


##############################
### Need to create 6 communities of 5 random species and 4 communities of 6 random species (note that random species are chosen from within the pool of 5 communities)
data=read.table(file="Acoustic_alldata_postDFA_AllCommunities.csv", header=T, sep=",", dec=".")
#data=read.table(file="Acoustic_alldata_postDFA.csv", header=T, sep=",", dec=".") #alternative with full dataset

# prepare a vector to storage the correct classification index values
randomStrat_Comm = c()

# create 6 random communities of 5 species
for (counter in 1:6){
  randomCom_index <- sample(length(levels(data$Species)),5, replace = F)
  randomCom <- droplevels(data[data$Species %in% levels(data$Species)[randomCom_index],])

  ################### Run cross-validated DFA based on all files available (ranging from 3 to 30 files per species) 
  ## DFA with random community data and CV
  ndata=nrow(randomCom)
  species.names = levels(randomCom$Species) #'Species' is the variable by which we want to discriminate drums
  nbirds=length(species.names)
  # There are nbirds = 92 species.

  resLDA=lda(Species~PC1+PC2+PC3+PC4+PC5+PC6, prior=rep(1,nbirds)/nbirds, data=randomCom, CV=T) # LDA computation. Important for information calculation: equal prior probabilities
  tab <- table(randomCom$Species, resLDA$class)
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
      ind.good <- (randomCom$Species == species.names[i])
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
    for (Y in 1:nbirds){ # then go through columns 
      if (X >= Y){newMat2[X,Y] <- NA}
      else {newMat2[X,Y] <- mean(c(cprob.tab.post[X,X]-cprob.tab.post[X,Y], cprob.tab.post[Y,Y]-cprob.tab.post[Y,X]))*100}
    }
  }

  rownames(newMat2) <- rownames(cprob.tab.post) # changed val_matrix into cprob.tab.post
  colnames(newMat2) <- colnames(cprob.tab.post) # changed val_matrix into cprob.tab.post

  #print (newMat2)#checkupline to visualize the correct classification index matrices that are created

  newMat2asVector <- na.omit(as.vector(newMat2))
  randomStrat_Comm <- c(randomStrat_Comm, newMat2asVector) # vector of 60 values (6 times 10 values, i.e. 10 values per community)
}

# create 4 random communities of 6 species
for (counter in 1:4){
  randomCom_index <- sample(length(levels(data$Species)),6, replace = F)
  randomCom <- droplevels(data[data$Species %in% levels(data$Species)[randomCom_index],])
  
  ################### Run cross-validated DFA based on all files available (ranging from 3 to 30 files per species) 
  ## DFA with random community data and CV
  ndata=nrow(randomCom)
  species.names = levels(randomCom$Species) #'Species' is the variable by which we want to discriminate drums
  nbirds=length(species.names)
  # There are nbirds = 92 species.
  
  resLDA=lda(Species~PC1+PC2+PC3+PC4+PC5+PC6, prior=rep(1,nbirds)/nbirds, data=randomCom, CV=T) # LDA computation. Important for information calculation: equal prior probabilities
  tab <- table(randomCom$Species, resLDA$class)
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
      ind.good <- (randomCom$Species == species.names[i])
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
    for (Y in 1:nbirds){ # then go through columns 
      if (X >= Y){newMat2[X,Y] <- NA}
      else {newMat2[X,Y] <- mean(c(cprob.tab.post[X,X]-cprob.tab.post[X,Y], cprob.tab.post[Y,Y]-cprob.tab.post[Y,X]))*100}
    }
  }
  
  rownames(newMat2) <- rownames(cprob.tab.post) # changed val_matrix into cprob.tab.post
  colnames(newMat2) <- colnames(cprob.tab.post) # changed val_matrix into cprob.tab.post
  
  #print (newMat2)#checkupline to visualize the correct classification index matrices that are created
  
  newMat2asVector <- na.omit(as.vector(newMat2)) # vector of 60 values (4 times 15 values, i.e. 15 values per community)
  randomStrat_Comm <- c(randomStrat_Comm, newMat2asVector) # vector of 120 values (6*10 + 4*15)
}

# Now plotting density curves based on CCI (correct classification index = the values stored in the vector created prior to the for loop)
CCI_randomStrat_Comm <- density(randomStrat_Comm,from=-20, to=107, n=128)
RandomDensity[iboot,] <- CCI_randomStrat_Comm$y

#par(mar=c(5,5,5,5))
#plot(CCI_randomStrat_Comm$x, CCI_randomStrat_Comm$y, xlim=c(-50,100), col="blue", type='l', xlab= "Correct classification Index", ylab= 'Probability density', lwd=4, cex.axis= 1.8, cex.lab=2.2)


# write.csv(randomStrat_Comm,'randomStrat_CommVector.csv')


##############################
### Finally, need to use real communities, creating 5 communities of 5 species (random sampling of 5 species once from each real community) and 5 communities of 6 species
### Note that, to use a similar number of correct classification values from '5-species communities' and '6-species communities', 
### we will need to downsample the number of values from the 5 communities with 6 species. Yet, those are created to provide a balanced sampling between real communities

Data_com1_CH=read.table(file="Acoustic_alldata_postDFA_community1_Switzerland.csv", header=T, sep=",",dec=".")
Data_com2_GU=read.table(file="Acoustic_alldata_postDFA_community2_Guatemala.csv", header=T, sep=",",dec=".")
Data_com3_MI=read.table(file="Acoustic_alldata_postDFA_community3_Minnesota.csv", header=T, sep=",",dec=".")
Data_com4_MA=read.table(file="Acoustic_alldata_postDFA_community4_Malaysia.csv", header=T, sep=",",dec=".")
Data_com5_FR=read.table(file="Acoustic_alldata_postDFA_community5_FrGuiana.csv", header=T, sep=",",dec=".")

# creating the 5 communities for 5 species:
# realcommunity CH
com1_species_index <- sample(length(levels(Data_com1_CH$Species)),5, replace = F)
com1_species <- droplevels(Data_com1_CH[Data_com1_CH$Species %in% levels(Data_com1_CH$Species)[com1_species_index],])

# realcommunity GU
com2_species_index <- sample(length(levels(Data_com2_GU$Species)),5, replace = F)
com2_species <- droplevels(Data_com2_GU[Data_com2_GU$Species %in% levels(Data_com2_GU$Species)[com2_species_index],])

# realcommunity MI
com3_species_index <- sample(length(levels(Data_com3_MI$Species)),5, replace = F)
com3_species <- droplevels(Data_com3_MI[Data_com3_MI$Species %in% levels(Data_com3_MI$Species)[com3_species_index],])

# realcommunity MA
com4_species_index <- sample(length(levels(Data_com4_MA$Species)),5, replace = F)
com4_species <- droplevels(Data_com4_MA[Data_com4_MA$Species %in% levels(Data_com4_MA$Species)[com4_species_index],])

# realcommunity FR
com5_species_index <- sample(length(levels(Data_com5_FR$Species)),5, replace = F)
com5_species <- droplevels(Data_com5_FR[Data_com5_FR$Species %in% levels(Data_com5_FR$Species)[com5_species_index],])


communityList <- list(com1_species,com2_species,com3_species,com4_species,com5_species)
# communityList <- list(RS_species) #checkup line to run one community at a time


# prepare a vector to storage the correct classification index values
real_Comm = c()

iteration <- '0' #checkupline
for (community in communityList){
  #print(iteration) #checkup line
  
  ################### For each community, run cross-validated DFA based on all files available (ranging from 3 to 30 files per species) 
  ## DFA with random community data and CV
  ndata=nrow(community)
  species.names = levels(community$Species) #'Species' is the variable by which we want to discriminate drums
  nbirds=length(species.names)
  # There are nbirds = 92 species.
  
  resLDA=lda(Species~PC1+PC2+PC3+PC4+PC5+PC6, prior=rep(1,nbirds)/nbirds, data=community, CV=T) # LDA computation. Important for information calculation: equal prior probabilities
  tab <- table(community$Species, resLDA$class)
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
      ind.good <- (community$Species == species.names[i])
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
    for (Y in 1:nbirds){ # then go through columns 
      if (X >= Y){newMat2[X,Y] <- NA}
      else {newMat2[X,Y] <- mean(c(cprob.tab.post[X,X]-cprob.tab.post[X,Y], cprob.tab.post[Y,Y]-cprob.tab.post[Y,X]))*100}
    }
  }
  
  rownames(newMat2) <- rownames(cprob.tab.post) # changed val_matrix into cprob.tab.post
  colnames(newMat2) <- colnames(cprob.tab.post) # changed val_matrix into cprob.tab.post
  
  #print(newMat2) #checkupline to visualize the correct classification index matrices that are created
  
  newMat2asVector <- na.omit(as.vector(newMat2))
  real_Comm <- c(real_Comm, newMat2asVector) # vector of 50 values (5 times 10 values, i.e. 10 values per community)
  
}


### In addition, creating the 5 communities with 6 species (note that those are separated from the loop run on '5-species communities' 
### in order to remove 25 random points from the correct classification values of the 6-species communities only)
# realcommunity CH_2
com1_species_index2 <- sample(length(levels(Data_com1_CH$Species)),6, replace = F)
com1_species2 <- droplevels(Data_com1_CH[Data_com1_CH$Species %in% levels(Data_com1_CH$Species)[com1_species_index2],])

# realcommunity GU_2
com2_species_index2 <- sample(length(levels(Data_com2_GU$Species)),6, replace = F)
com2_species2 <- droplevels(Data_com2_GU[Data_com2_GU$Species %in% levels(Data_com2_GU$Species)[com2_species_index2],])

# realcommunity MI_2
com3_species_index2 <- sample(length(levels(Data_com3_MI$Species)),6, replace = F)
com3_species2 <- droplevels(Data_com3_MI[Data_com3_MI$Species %in% levels(Data_com3_MI$Species)[com3_species_index2],])

# realcommunity MA_2
com4_species_index2 <- sample(length(levels(Data_com4_MA$Species)),6, replace = F)
com4_species2 <- droplevels(Data_com4_MA[Data_com4_MA$Species %in% levels(Data_com4_MA$Species)[com4_species_index2],])

# realcommunity FR_2
com5_species_index2 <- sample(length(levels(Data_com5_FR$Species)),6, replace = F)
com5_species2 <- droplevels(Data_com5_FR[Data_com5_FR$Species %in% levels(Data_com5_FR$Species)[com5_species_index2],])

communityList2 <- list(com1_species2,com2_species2,com3_species2,com4_species2,com5_species2)
# communityList <- list(RS_species) #checkup line to run one community at a time

real_Comm_tmp = c()

iteration <- '0' #checkupline
for (community in communityList2){
  #print(iteration) #checkup line
  
  ################### For each community, run cross-validated DFA based on all files available (ranging from 3 to 30 files per species) 
  ## DFA with random community data and CV
  ndata=nrow(community)
  species.names = levels(community$Species) #'Species' is the variable by which we want to discriminate drums
  nbirds=length(species.names)
  # There are nbirds = 92 species.
  
  resLDA=lda(Species~PC1+PC2+PC3+PC4+PC5+PC6, prior=rep(1,nbirds)/nbirds, data=community, CV=T) # LDA computation. Important for information calculation: equal prior probabilities
  tab <- table(community$Species, resLDA$class)
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
      ind.good <- (community$Species == species.names[i])
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
    for (Y in 1:nbirds){ # then go through columns 
      if (X >= Y){newMat2[X,Y] <- NA}
      else {newMat2[X,Y] <- mean(c(cprob.tab.post[X,X]-cprob.tab.post[X,Y], cprob.tab.post[Y,Y]-cprob.tab.post[Y,X]))*100}
    }
  }
  
  rownames(newMat2) <- rownames(cprob.tab.post) # changed val_matrix into cprob.tab.post
  colnames(newMat2) <- colnames(cprob.tab.post) # changed val_matrix into cprob.tab.post
  
  #print(newMat2) #checkupline to visualize the correct classification index matrices that are created
  
  newMat2asVector <- na.omit(as.vector(newMat2))
  real_Comm_tmp <- c(real_Comm_tmp, newMat2asVector) # vector of 75 values (5 times 15 values, i.e. 15 values per community)

}

newMat2asVector_downsampled <- sample(real_Comm_tmp,50, replace = F)
real_Comm <- c(real_Comm, newMat2asVector_downsampled) # vector of 100 values (5 times 10 values + "5 times 15 values - 25 values")

# Now plotting density curves based on CCI (correct classification index = the values stored in the vector created prior to the for loop)
CCI_real_Comm <- density(real_Comm, from=-20, to=107, n=128)
RealDensity[iboot,] <- CCI_real_Comm$y

#par(mar=c(5,5,5,5))
#plot(CCI_real_Comm$x, CCI_real_Comm$y, xlim=c(-50,100), col="blue", type='l', xlab= "Correct classification Index", ylab= 'Probability density', lwd=4, cex.axis= 1.8, cex.lab=2.2)


# write.csv(real_Comm,'real_CommVector.csv')

}
#   End of boot strap loop

# Make a gg plot for Density Functions
RealDensityMean <- colMeans(RealDensity)
RandomDensityMean <- colMeans(RandomDensity)
SameDensityMean <- colMeans(SameDensity)
DiffDensityMean <- colMeans(DiffDensity)

RealDensitySD <- apply(RealDensity, 2, sd)
RandomDensitySD <- apply(RandomDensity, 2, sd)
SameDensitySD <- apply(SameDensity, 2, sd)
DiffDensitySD <- apply(DiffDensity, 2, sd)

dataPlot <- data.frame(xrange, RealDensityMean, RandomDensityMean, SameDensityMean, DiffDensityMean, RealDensitySD, RandomDensitySD, SameDensitySD, DiffDensitySD)

# Use ggplot to make nice error bar plots
plot.density<- ggplot(dataPlot, aes(x=xrange))

plot.density <- plot.density + geom_line(aes(y=RealDensityMean, colour='Real') ) 
plot.density <- plot.density + geom_ribbon(aes(ymin=RealDensityMean-RealDensitySD, ymax=RealDensityMean-RealDensitySD),alpha=.3, fill = 'black')

plot.density <- plot.density +
  geom_line(aes(y=RandomDensityMean, colour='Random'), lty='dashed') 
plot.density <- plot.density +
  geom_ribbon(aes(ymin=RandomDensityMean-RandomDensitySD, ymax=RandomDensityMean+RandomDensitySD), alpha=.3, fill = 'grey')

plot.density <- plot.density +
  geom_line(aes(x=xrange, y=SameDensityMean, colour='Same'), lty='dashed') 
plot.density <- plot.density +
  geom_ribbon(aes(ymin=SameDensityMean-SameDensitySD, ymax=SameDensityMean+SameDensitySD), alpha=.3, fill = 'blue')

plot.density <- plot.density +
  geom_line(aes(x=xrange, y=DiffDensityMean, colour='Diff'), lty='dashed')
plot.density <- plot.density +
  geom_ribbon(aes(ymin=DiffDensityMean-DiffDensitySD, ymax=DiffDensityMean+DiffDensitySD),alpha=.3, fill = 'red')

plot.density <- plot.density + labs(x="Classification Index", y = "Prob") + ylim(c(0.0, 0.04)) +  xlim(c(-25, 100)) +
  theme_classic() +
  theme(legend.position=c(0.2,0.6)) +
  scale_colour_manual(values = c('red', 'grey', 'black', 'blue'), name=NULL)

plot.density

#ggsave("ClassificationIndexSimulationsProb.pdf", device = "pdf", width = 15, height = 10, units = "cm")

######### Repeat for Cumulative distributions  ###############
dx <- xrange[2]-xrange[1]
# Make space to store cumulative distributions
RealCum <- matrix(nrow=nBoot, ncol=nx)
RandomCum <- matrix(nrow=nBoot, ncol=nx)
SameCum <- matrix(nrow=nBoot, ncol=nx)
DiffCum <- matrix(nrow=nBoot, ncol=nx)

for (i in 1:nBoot) {
RealCum[i,] <- cumsum(RealDensity[i,])*dx
RandomCum[i,] <- cumsum(RandomDensity[i,])*dx
SameCum[i,] <- cumsum(SameDensity[i,])*dx
DiffCum[i,] <- cumsum(DiffDensity[i,])*dx
}

# Make a gg plot for Cum Functions
RealCumMean <- colMeans(RealCum)
RandomCumMean <- colMeans(RandomCum)
SameCumMean <- colMeans(SameCum)
DiffCumMean <- colMeans(DiffCum)

RealCumSD <- apply(RealCum, 2, sd)
RandomCumSD <- apply(RandomCum, 2, sd)
SameCumSD <- apply(SameCum, 2, sd)
DiffCumSD <- apply(DiffCum, 2, sd)

dataPlot <- data.frame(xrange, RealCumMean, RandomCumMean, SameCumMean, DiffCumMean, RealCumSD, RandomCumSD, SameCumSD, DiffCumSD)


# Use ggplot to make nice error bar plots
plot.density<- ggplot(dataPlot, aes(x=xrange))

plot.density <- plot.density + geom_line(aes(y=RealCumMean, colour='real'), lwd=1) 
plot.density <- plot.density + geom_ribbon(aes(ymin=RealCumMean-RealCumSD, ymax=RealCumMean-RealCumSD),alpha=.3, fill = 'black')

plot.density <- plot.density +
  geom_line(aes(y=RandomCumMean, colour='Random'), lty='dashed') 
plot.density <- plot.density +
  geom_ribbon(aes(ymin=RandomCumMean-RandomCumSD, ymax=RandomCumMean+RandomCumSD), alpha=.3, fill = 'grey')

plot.density <- plot.density +
  geom_line(aes(x=xrange, y=SameCumMean, colour='Same'), lty='dashed') 
plot.density <- plot.density +
  geom_ribbon(aes(ymin=SameCumMean-SameCumSD, ymax=SameCumMean+SameCumSD), alpha=.3, fill = 'blue')

plot.density <- plot.density +
  geom_line(aes(x=xrange, y=DiffCumMean, colour='Diff'), lty='dashed')
plot.density <- plot.density +
  geom_ribbon(aes(ymin=DiffCumMean-DiffCumSD, ymax=DiffCumMean+DiffCumSD),alpha=.3, fill = 'red')

plot.density <- plot.density + labs(x="Classification Index", y = "Cumulative Probability Denisty") + ylim(c(0.0, 1.0)) +  xlim(c(-20, 105)) +
  theme_bw() +
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=1)) +
#  theme_classic() +
  theme(axis.title.y = element_text(size = rel(1.8), angle = 90)) +
  theme(axis.title.x = element_text(size = rel(1.8), angle = 00)) +
  theme(axis.text.x = element_text(angle = 00, hjust = 0.5, size=17, color="black")) +
  theme(axis.text.y = element_text(angle = 90, hjust = 0.5, size=17, color="black")) +
  theme(legend.justification=c(0,1), legend.position=c(0,1),legend.box.margin=margin(0,0,0,0), legend.background = element_blank()) +
  scale_colour_manual(values = c('red', 'grey', 'black', 'blue'), labels = c("Different Type", "Random", "Real Communities","Same Type"), name=NULL)

plot.density

#ggsave("ClassificationIndexSimulationsCum.pdf", device = "pdf", width = 12, height = 12, units = "cm")


###############################
# Plotting  graphs 
# First need to check what range to use for y-axis
max_Y_val <- max(CCI_real_Comm$y,CCI_randomStrat_Comm$y, CCI_diffStrat_Comm$y, CCI_sameStrat_Comm$y)

#pdf("DensityPlot_correct_classif.pdf")
par(mar=c(5,5,5,5))
plot(CCI_real_Comm$x, CCI_real_Comm$y, xlim=c(-50,100), ylim=c(0,max_Y_val), col="black", type='l', xlab= "Correct classification Index", ylab= 'Probability density', lwd=5, cex.axis= 1.8, cex.lab=2.2)

lines(CCI_randomStrat_Comm$x, CCI_randomStrat_Comm$y, xlim=c(-50,100), col="darkgray", type='l', lty=2, lwd=5, cex.axis= 1.8, cex.lab=2.2)
lines(CCI_diffStrat_Comm$x, CCI_diffStrat_Comm$y, xlim=c(-50,100), col="red", type='l', lty=2, lwd=5, cex.axis= 1.8, cex.lab=2.2)
lines(CCI_sameStrat_Comm$x, CCI_sameStrat_Comm$y, xlim=c(-50,100), col="blue", type='l', lty=2, lwd=5, cex.axis= 1.8, cex.lab=2.2)

lines(CCI_real_Comm$x, CCI_real_Comm$y, xlim=c(-50,100), type='l', lwd=5, col='black', lty=1)


legend("topleft", c("Real Communities","Random","Same Strategy", "Different Strategy"), cex=1, 
       col=c("black","darkgray", "blue", "red"), lwd = 3, lty=2)
# Repeat legend to overlay previous while leaving Random community in dashed
legend("topleft","Real Communities", cex=1, bty = 'n',
       col=c("black"), lwd = 3, lty=1)
#dev.off()

  
#pdf("CumProb_correct_classif.pdf")
par(mar=c(5,5,5,5))
cum <- cumsum(CCI_real_Comm$y*(CCI_real_Comm$x[2]-CCI_real_Comm$x[1]))
plot(CCI_real_Comm$x, cum, xlim=c(-25,100), ylim=c(0,1.0), col="black", type='l', xlab= "Correct classification Index", ylab= 'Cum Prob', lwd=5, cex.axis= 1.8, cex.lab=2.2)

cum <- cumsum(CCI_randomStrat_Comm$y*(CCI_randomStrat_Comm$x[2]-CCI_randomStrat_Comm$x[1]))
lines(CCI_randomStrat_Comm$x, cum, xlim=c(-25,100), col="darkgray", type='l', lty=2, lwd=5, cex.axis= 1.8, cex.lab=2.2)

cum <- cumsum(CCI_diffStrat_Comm$y*(CCI_diffStrat_Comm$x[2]-CCI_diffStrat_Comm$x[1]))
lines(CCI_diffStrat_Comm$x, cum, xlim=c(-25,100), col="red", type='l', lty=2, lwd=5, cex.axis= 1.8, cex.lab=2.2)

cum <- cumsum(CCI_sameStrat_Comm$y*(CCI_sameStrat_Comm$x[2]-CCI_sameStrat_Comm$x[1]))
lines(CCI_sameStrat_Comm$x, cum, xlim=c(-25,100), col="blue", type='l', lty=2, lwd=5, cex.axis= 1.8, cex.lab=2.2)


legend("topleft", c("Real Communities","Random","Same Strategy", "Different Strategy"), cex=1, 
       col=c("black","darkgray", "blue", "red"), lwd = 3, lty=2)
# Repeat legend to overlay previous while leaving Random community in dashed
legend("topleft","Real Communities", cex=1, bty = 'n',
       col=c("black"), lwd = 3, lty=1)
#dev.off()

  