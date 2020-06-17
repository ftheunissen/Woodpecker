# Set working directory to folder with your files (change accordingly)
setwd("../10_Communities_densitycurves/AcousticDistance/")
rm(list=ls())#empty the workspace

library(ggplot2)
# Import all datafiles
AcousticMatrix=read.table(file="Acousticdist_matrix.csv", header=T, sep=",", dec=".", row.names = 1)

AcousticMatrix_Com1=read.table(file="Acousticdist_matrix_community1_Switzerland.csv", header=T, sep=",", dec=".", row.names = 1)
AcousticMatrix_Com2=read.table(file="Acousticdist_matrix_community2_Guatemala.csv", header=T, sep=",", dec=".", row.names = 1)
AcousticMatrix_Com3=read.table(file="Acousticdist_matrix_community3_Minnesota.csv", header=T, sep=",", dec=".", row.names = 1)
AcousticMatrix_Com4=read.table(file="Acousticdist_matrix_community4_Malaysia.csv", header=T, sep=",", dec=".", row.names = 1)
AcousticMatrix_Com5=read.table(file="Acousticdist_matrix_community5_FrGuiana.csv", header=T, sep=",", dec=".", row.names = 1)
AcousticMatrixAllCom=read.table(file="Acousticdist_matrix_ALLcommunities.csv", header=T, sep=",", dec=".", row.names = 1)

######### To plot densities, cannot used a dataframe format
# First, need to format from dataframe to matrix
AcousticMatrix <- as.matrix(AcousticMatrix)
AcousticMatrix_Com1 <- as.matrix(AcousticMatrix_Com1)
AcousticMatrix_Com2 <- as.matrix(AcousticMatrix_Com2)
AcousticMatrix_Com3 <- as.matrix(AcousticMatrix_Com3)
AcousticMatrix_Com4 <- as.matrix(AcousticMatrix_Com4)
AcousticMatrix_Com5 <- as.matrix(AcousticMatrix_Com5)
AcousticMatrixAllCom <- as.matrix(AcousticMatrixAllCom)

# Then possible to extract all values into a single dimension array (which we can then plot)
AcousticsGlobal <- na.omit(as.vector(AcousticMatrix))
AcousticsCom1 <- na.omit(as.vector(AcousticMatrix_Com1))
AcousticsCom2 <- na.omit(as.vector(AcousticMatrix_Com2))
AcousticsCom3 <- na.omit(as.vector(AcousticMatrix_Com3))
AcousticsCom4 <- na.omit(as.vector(AcousticMatrix_Com4))
AcousticsCom5 <- na.omit(as.vector(AcousticMatrix_Com5))
AcousticsAllCom <- na.omit(as.vector(AcousticMatrixAllCom))

# Then compile the kernel density estimation for each dataframe (whole or communities)
densityGlobal <- density(AcousticsGlobal)
densityCom1 <- density(AcousticsCom1)
densityCom2 <- density(AcousticsCom2)
densityCom3 <- density(AcousticsCom3)
densityCom4 <- density(AcousticsCom4)
densityCom5 <- density(AcousticsCom5)
densityAllCom <- density(AcousticsAllCom)

#pdf("DensityPlot_AcousticDistance.pdf")
par(mar=c(5,5,5,5))
plot(densityGlobal$x, densityGlobal$y, col="black", type='l', xlab= "Acoustic Distance", ylim = c(0,0.5), ylab= 'Probability Density', lwd=4, cex.axis= 1.8, cex.lab=2.2)
lines(densityCom1$x, densityCom1$y,type= 'l',col="darkorange", lwd=3)
lines(densityCom2$x, densityCom2$y,type= 'l',col="bisque1", lwd=3)     
lines(densityCom3$x, densityCom3$y,type= 'l',col="darkred", lwd=3)     
lines(densityCom4$x, densityCom4$y,type= 'l',col="darkgoldenrod1", lwd=3)     
lines(densityCom5$x, densityCom5$y,type= 'l',col="darkgoldenrod4", lwd=3)
#lines(densityAllCom$x, densityAllCom$y/max(densityAllCom$y),type = 'l', col = 'darkgray', lwd=4)
lines(densityGlobal$x, densityGlobal$y,type='l', col="black", lwd=5)

legend("topright", c("Global","Switzerland","Guatemala", "Minnesota", "Malaysia", "Fr.Guiana"), cex=1, bty = 'n',
       col=c("black","darkorange", "bisque1", "darkred", "darkgoldenrod1","darkgoldenrod4"), lwd = 3, lty=1)
#dev.off()


# Plot next as a cumulative probability distribution
#pdf("CumProb_AcousticDistance.pdf")
par(mar=c(5,5,5,5))
cum <- cumsum(densityGlobal$y*(densityGlobal$x[2]-densityGlobal$x[1]))
plot(densityGlobal$x, cum, xlim=c(0,15), ylim=c(0,1.0), col="black", type='l', xlab= "Acoustic Distance", ylab= 'Cumulative Probability Density', lwd=4, cex.axis= 1.8, cex.lab=2.2)

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
legend("bottomright", c("Global","Switzerland","Guatemala", "Minnesota", "Malaysia", "Fr.Guiana"), cex=0.8, 
       col=c("black","darkorange", "bisque1", "darkred", "darkgoldenrod1","darkgoldenrod4"), lwd = 3, lty=1, bty = "n")
#dev.off()




wilcoxTest1 <- wilcox.test(x = AcousticsGlobal, y = AcousticsCom1)
qnorm(wilcoxTest1$p.value/2)
wilcoxTest2 <- wilcox.test(x = AcousticsGlobal, y = AcousticsCom2)
qnorm(wilcoxTest2$p.value/2)
wilcoxTest3 <- wilcox.test(x = AcousticsGlobal, y = AcousticsCom3)
qnorm(wilcoxTest3$p.value/2)
wilcoxTest4 <- wilcox.test(x = AcousticsGlobal, y = AcousticsCom4)
qnorm(wilcoxTest4$p.value/2)
wilcoxTest5 <- wilcox.test(x = AcousticsGlobal, y = AcousticsCom5)
qnorm(wilcoxTest5$p.value/2)

p.adjust(as.vector(c(wilcoxTest1$p.value,wilcoxTest2$p.value,wilcoxTest3$p.value,wilcoxTest4$p.value,wilcoxTest5$p.value)),method = 'bonferroni', n=5) 



######################################
# Above was the analysis to plot the probability densities from real communities. Now using simulated communities for comparison:

data=read.table(file="Acoustic_species_postDFA.csv", header=T, sep=",",dec=".")
AcousticMatrix=read.table(file="Acousticdist_matrix.csv", header=T, sep=",", dec=".", row.names = 1)

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
xrange <- seq(-2.85, 16.3, by=0.15) ######## here get 128 values to match what the density function outputs
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


# prepare a vector to store the correct classification index values
sameStrat_Comm = c()

for (community in communityList){
  ACMatrix <- as.matrix(AcousticMatrix[rownames(AcousticMatrix)%in% community$Species, colnames(AcousticMatrix)%in% community$Species]) # Rownames(AcousticMatrix) should be = to colnames(AcousticMatrix)
  acousticDistanceVector <- na.omit(as.vector(ACMatrix))
  sameStrat_Comm <- c(sameStrat_Comm, acousticDistanceVector) # vector of 120 values (6 times 10 values, i.e. 10 values per drum type ; + 4 times 15 values)
  
  print(ACMatrix)
  
}

# Now plotting density curves based on CCI (correct classification index = the values stored in the vector created prior to the for loop)
#AcDist_sameStrat_Comm <- density(sameStrat_Comm)
AcDist_sameStrat_Comm <- density(sameStrat_Comm, from=-2.85, to=16.3, n=128)
SameDensity[iboot,] <- AcDist_sameStrat_Comm$y


#par(mar=c(5,5,5,5))
#plot(AcDist_sameStrat_Comm$x, AcDist_sameStrat_Comm$y, xlim=c(0,15), col="red", type='l', xlab= "Correct classification Index", ylab= 'Probability density', lwd=4, cex.axis= 1.8, cex.lab=2.2)

#write.csv(sameStrat_Comm,'sameStrat_CommVector_AcDist.csv')






##############################
### Then need to create 6 communities of 5 species and 4 communities of 6 species where each species has a different drum type
data=read.table(file="Acoustic_species_postDFA.csv", header=T, sep=",",dec=".")
AcousticMatrix=read.table(file="Acousticdist_matrix.csv", header=T, sep=",", dec=".", row.names = 1)

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
  
  ACMatrix <- as.matrix(AcousticMatrix[rownames(AcousticMatrix)%in% FivespeciesList$Species, colnames(AcousticMatrix)%in% FivespeciesList$Species]) # Rownames(AcousticMatrix) should be = to colnames(AcousticMatrix)
  acousticDistanceVector <- na.omit(as.vector(ACMatrix))
  diffStrat_Comm <- c(diffStrat_Comm, acousticDistanceVector) # vector of 60 values (6 times 10 values, i.e. 10 values per drum type)
  
  print(ACMatrix)
  
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
  
  ACMatrix <- as.matrix(AcousticMatrix[rownames(AcousticMatrix)%in% SixspeciesList$Species, colnames(AcousticMatrix)%in% SixspeciesList$Species]) # Rownames(AcousticMatrix) should be = to colnames(AcousticMatrix)
  acousticDistanceVector <- na.omit(as.vector(ACMatrix))
  diffStrat_Comm <- c(diffStrat_Comm, acousticDistanceVector) # vector of 120 values (6 times 10 values, i.e. 10 values per drum type ; + 4 times 15 values)
  
  print(ACMatrix)
  
}

# Now plotting density curves based on CCI (correct classification index = the values stored in the vector created prior to the for loop)
#AcDist_diffStrat_Comm <- density(diffStrat_Comm)
AcDist_diffStrat_Comm <- density(diffStrat_Comm, from=-2.85, to=16.3, n=128)
DiffDensity[iboot,] <- AcDist_diffStrat_Comm$y

#par(mar=c(5,5,5,5))
#plot(AcDist_diffStrat_Comm$x, AcDist_diffStrat_Comm$y, xlim=c(0,15), col="red", type='l', xlab= "Correct classification Index", ylab= 'Probability density', lwd=4, cex.axis= 1.8, cex.lab=2.2)

#write.csv(diffStrat_Comm,'diffStrat_CommVector_AcDist.csv')


  



##############################
### Need to create 6 communities of 5 random species and 4 communities of 6 random species (note that random species are chosen from within the pool of 5 communities)
data=read.table(file="Acoustic_species_postDFA_AllComm.csv", header=T, sep=",", dec=".")
#data=read.table(file="Acoustic_species_postDFA.csv", header=T, sep=",", dec=".") # alternative with full dataset

AcousticMatrix=read.table(file="Acousticdist_matrix.csv", header=T, sep=",", dec=".", row.names = 1)


# prepare a vector to storage the correct classification index values
randomStrat_Comm = c()

# create 6 random communities of 5 species
for (counter in 1:6){
  randomCom_index <- sample(length(levels(data$Species)),5, replace = F)
  randomCom <- droplevels(data[data$Species %in% levels(data$Species)[randomCom_index],])
  
  ACMatrix <- as.matrix(AcousticMatrix[rownames(AcousticMatrix)%in% randomCom$Species, colnames(AcousticMatrix)%in% randomCom$Species]) # Rownames(AcousticMatrix) should be = to colnames(AcousticMatrix)
  acousticDistanceVector <- na.omit(as.vector(ACMatrix))
  randomStrat_Comm <- c(randomStrat_Comm, acousticDistanceVector) # vector of 60 values (6 times 10 values, i.e. 10 values per community)
  
  print (ACMatrix)#checkupline to visualize the correct classification index matrices that are created
  
}

# create 4 random communities of 6 species
for (counter in 1:4){
  randomCom_index <- sample(length(levels(data$Species)),6, replace = F)
  randomCom <- droplevels(data[data$Species %in% levels(data$Species)[randomCom_index],])
  
  ACMatrix <- as.matrix(AcousticMatrix[rownames(AcousticMatrix)%in% randomCom$Species, colnames(AcousticMatrix)%in% randomCom$Species]) # Rownames(AcousticMatrix) should be = to colnames(AcousticMatrix)
  acousticDistanceVector <- na.omit(as.vector(ACMatrix))
  randomStrat_Comm <- c(randomStrat_Comm, acousticDistanceVector) # vector of 60 values (4 times 15 values, i.e. 15 values per community)
  
  print (ACMatrix)#checkupline to visualize the correct classification index matrices that are created
}

# Now plotting density curves based on CCI (correct classification index = the values stored in the vector created prior to the for loop)
AcDist_randomStrat_Comm <- density(randomStrat_Comm)
AcDist_randomStrat_Comm <- density(randomStrat_Comm,from=-2.85, to=16.3, n=128)
RandomDensity[iboot,] <- AcDist_randomStrat_Comm$y

#par(mar=c(5,5,5,5))
#plot(AcDist_randomStrat_Comm$x, AcDist_randomStrat_Comm$y, xlim=c(0,15), col="red", type='l', xlab= "Correct classification Index", ylab= 'Probability density', lwd=4, cex.axis= 1.8, cex.lab=2.2)


#write.csv(randomStrat_Comm,'randomStrat_CommVector_AcDist.csv')


  
  
  




##############################
### Finally, need to use real communities, creating 5 communities of 5 species (random sampling of 5 species once from each real community) and 5 communities of 6 species
### Note that, to use a similar number of correct classification values from '5-species communities' and '6-species communities', 
### we will need to downsample the number of values from the 5 communities with 6 species. Yet, those are created to provide a balanced sampling between real communities

AcousticMatrix_Com1=read.table(file="Acousticdist_matrix_community1_Switzerland.csv", header=T, sep=",", dec=".", row.names = 1)
AcousticMatrix_Com2=read.table(file="Acousticdist_matrix_community2_Guatemala.csv", header=T, sep=",", dec=".", row.names = 1)
AcousticMatrix_Com3=read.table(file="Acousticdist_matrix_community3_Minnesota.csv", header=T, sep=",", dec=".", row.names = 1)
AcousticMatrix_Com4=read.table(file="Acousticdist_matrix_community4_Malaysia.csv", header=T, sep=",", dec=".", row.names = 1)
AcousticMatrix_Com5=read.table(file="Acousticdist_matrix_community5_FrGuiana.csv", header=T, sep=",", dec=".", row.names = 1)

real_Comm <- c()
# first choosing 5 species among each of the real communities. Thus creating 50 values (5 times 10 values, i.e. 10 values per community)
com1_5sp_index <- sample(length(colnames(AcousticMatrix_Com1)),5, replace = F)
com1_5sp <- na.omit(as.vector(as.matrix(AcousticMatrix_Com1[c(sort(com1_5sp_index)),c(sort(com1_5sp_index))])))
print(AcousticMatrix_Com1[c(sort(com1_5sp_index)),c(sort(com1_5sp_index))])
  
com2_5sp_index <- sample(length(colnames(AcousticMatrix_Com2)),5, replace = F)
com2_5sp <- na.omit(as.vector(as.matrix(AcousticMatrix_Com2[c(sort(com2_5sp_index)),c(sort(com2_5sp_index))])))
print(AcousticMatrix_Com2[c(sort(com2_5sp_index)),c(sort(com2_5sp_index))])

com3_5sp_index <- sample(length(colnames(AcousticMatrix_Com3)),5, replace = F)
com3_5sp <- na.omit(as.vector(as.matrix(AcousticMatrix_Com3[c(sort(com3_5sp_index)),c(sort(com3_5sp_index))])))
print(AcousticMatrix_Com3[c(sort(com3_5sp_index)),c(sort(com3_5sp_index))])

com4_5sp_index <- sample(length(colnames(AcousticMatrix_Com4)),5, replace = F)
com4_5sp <- na.omit(as.vector(as.matrix(AcousticMatrix_Com4[c(sort(com4_5sp_index)),c(sort(com4_5sp_index))])))
print(AcousticMatrix_Com4[c(sort(com4_5sp_index)),c(sort(com4_5sp_index))])  

com5_5sp_index <- sample(length(colnames(AcousticMatrix_Com5)),5, replace = F)
com5_5sp <- na.omit(as.vector(as.matrix(AcousticMatrix_Com5[c(sort(com5_5sp_index)),c(sort(com5_5sp_index))])))
print(AcousticMatrix_Com5[c(sort(com5_5sp_index)),c(sort(com5_5sp_index))])

real_Comm <- c(com1_5sp, com2_5sp, com3_5sp, com4_5sp, com5_5sp)


# then choosing 6 species among each of the real communities. Thus creating 50 values (5 times 15 values, minus random selection of 25 points to balance wight of 5specie and 6species in the set of acoustic distance values)
com1_6sp_index <- sample(length(colnames(AcousticMatrix_Com1)),6, replace = F)
com1_6sp <- na.omit(as.vector(as.matrix(AcousticMatrix_Com1[c(sort(com1_6sp_index)),c(sort(com1_6sp_index))])))
print(AcousticMatrix_Com1[c(sort(com1_6sp_index)),c(sort(com1_6sp_index))])

com2_6sp_index <- sample(length(colnames(AcousticMatrix_Com2)),6, replace = F)
com2_6sp <- na.omit(as.vector(as.matrix(AcousticMatrix_Com2[c(sort(com2_6sp_index)),c(sort(com2_6sp_index))])))
print(AcousticMatrix_Com2[c(sort(com2_6sp_index)),c(sort(com2_6sp_index))])

com3_6sp_index <- sample(length(colnames(AcousticMatrix_Com3)),6, replace = F)
com3_6sp <- na.omit(as.vector(as.matrix(AcousticMatrix_Com3[c(sort(com3_6sp_index)),c(sort(com3_6sp_index))])))
print(AcousticMatrix_Com3[c(sort(com3_6sp_index)),c(sort(com3_6sp_index))])

com4_6sp_index <- sample(length(colnames(AcousticMatrix_Com4)),6, replace = F)
com4_6sp <- na.omit(as.vector(as.matrix(AcousticMatrix_Com4[c(sort(com4_6sp_index)),c(sort(com4_6sp_index))])))
print(AcousticMatrix_Com4[c(sort(com4_6sp_index)),c(sort(com4_6sp_index))])

com5_6sp_index <- sample(length(colnames(AcousticMatrix_Com5)),6, replace = F)
com5_6sp <- na.omit(as.vector(as.matrix(AcousticMatrix_Com5[c(sort(com5_6sp_index)),c(sort(com5_6sp_index))])))
print(AcousticMatrix_Com5[c(sort(com5_6sp_index)),c(sort(com5_6sp_index))])

real_Comm_tmp <- c(com1_6sp, com2_6sp, com3_6sp, com4_6sp, com5_6sp)

real_Comm_tmp_downsampled <- sample(real_Comm_tmp,50, replace = F)
real_Comm <- c(real_Comm, real_Comm_tmp_downsampled) # vector of 100 values (5 times 10 values + "5 times 15 values - 25 values")

# Now plotting density curves based on CCI (correct classification index = the values stored in the vector created prior to the for loop)
#Ac_Dist_real_Comm <- density(real_Comm)
Ac_Dist_real_Comm <- density(real_Comm, from=-2.85, to=16.3, n=128)
RealDensity[iboot,] <- Ac_Dist_real_Comm$y

#par(mar=c(5,5,5,5))
#plot(Ac_Dist_real_Comm$x, Ac_Dist_real_Comm$y, xlim=c(0,15), col="red", type='l', xlab= "Correct classification Index", ylab= 'Probability density', lwd=4, cex.axis= 1.8, cex.lab=2.2)


#write.csv(real_Comm,'real_CommVector_AcDist.csv')

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

plot.density <- plot.density + geom_line(aes(y=RealDensityMean, colour='Real'),lwd=1.5 ) 
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

plot.density <- plot.density + labs(x="Acoustic distance", y = "Probability density") + ylim(c(-0.01, 0.35)) +  xlim(c(0, 15)) +
  theme_classic() +
  theme(legend.position=c(0.8,0.8)) +
  scale_colour_manual(values = c('red', 'grey', 'black', 'blue'), name=NULL)

plot.density

#ggsave("AcDistanceSimulationsProb.pdf", device = "pdf", width = 15, height = 10, units = "cm")

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

plot.density <- plot.density + geom_line(aes(y=RealCumMean, colour='Real'), lwd=1) 
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


plot.density <- plot.density + labs(x="Acoustic Distance", y = "Cumulative Probability Denisty") + ylim(c(0.0, 1.0)) +  xlim(c(0, 15)) +
  theme_bw() +
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=1)) +
  #  theme_classic() +
  theme(axis.title.y = element_text(size = rel(1.8), angle = 90)) +
  theme(axis.title.x = element_text(size = rel(1.8), angle = 00)) +
  theme(axis.text.x = element_text(angle = 00, hjust = 0.5, size=17, color="black")) +
  theme(axis.text.y = element_text(angle = 90, hjust = 0.5, size=17, color="black")) +
  theme(legend.justification=c(1,0), legend.position=c(1,0),legend.box.margin=margin(0,0,0,0), legend.background = element_blank()) +
  scale_colour_manual(values = c('red', 'grey', 'black', 'blue'), labels = c("Different Type", "Random", "Real Communities","Same Type"), name=NULL)

plot.density

#ggsave("AcDistanceSimulationsCum.pdf", device = "pdf", width = 12, height = 12, units = "cm")





###############################
# Plotting  graphs 
# First need to check what range to use for y-axis
max_Y_val <- max(AcDist_sameStrat_Comm$y,AcDist_diffStrat_Comm$y, AcDist_randomStrat_Comm$y, Ac_Dist_real_Comm$y)

#pdf("DensityPlot_AcousticDistance_simCom.pdf")
par(mar=c(5,5,5,5))
plot(Ac_Dist_real_Comm$x, Ac_Dist_real_Comm$y, xlim=c(0,15), ylim=c(0,max_Y_val), col="black", type='l', xlab= "Acoustic distance", ylab= 'Probability density', lwd=5, cex.axis= 1.8, cex.lab=2.2)

lines(AcDist_randomStrat_Comm$x, AcDist_randomStrat_Comm$y, xlim=c(0,15), col="darkgray", type='l', lty=2, lwd=5, cex.axis= 1.8, cex.lab=2.2)
lines(AcDist_diffStrat_Comm$x, AcDist_diffStrat_Comm$y, xlim=c(0,15), col="red", type='l', lty=2, lwd=5, cex.axis= 1.8, cex.lab=2.2)
lines(AcDist_sameStrat_Comm$x, AcDist_sameStrat_Comm$y, xlim=c(0,15), col="blue", type='l', lty=2, lwd=5, cex.axis= 1.8, cex.lab=2.2)

lines(Ac_Dist_real_Comm$x, Ac_Dist_real_Comm$y, xlim=c(0,15), type='l', lwd=5, col='black', lty=1)


legend("topright", c("Real Communities","Random","Same Strategy", "Different Strategy"), cex=1, 
       col=c("black","darkgray", "blue", "red"), lwd = 3, lty=2)
# Repeat legend to overlay previous while leaving Random community in dashed
legend("topright","Real Communities", cex=1, bty = 'n',
       col=c("black"), lwd = 3, lty=1)
#dev.off()

  
