
setwd("../10_Communities_densitycurves/AcousticDistance/")
rm(list=ls())#empty the workspace
data=read.table(file="Acoustic_species_postDFA.csv", header=T, sep=",",dec=".", row.names = 2)

# Creating acoustic distance matrix based on euclidean distances
AcousticDistances <- dist(data[,15:36])
AcousticDistances <- as.matrix(AcousticDistances)

for (X in 1:dim(AcousticDistances)[1]){ # First go through rows
  for (Y in 1:dim(AcousticDistances)[2]){ # then go through columns 
    if (X >= Y){AcousticDistances[X,Y] <- 'NA'} # set one side (in this case the lower) of the diagonal to NAs (two mirroed half-matrices here thus repetition of the same data)
  }
}
write.csv(AcousticDistances, "Acousticdist_matrix.csv")


##############################
# Control lines to see how eucliden distances are calculated from vectors/lists and make sure the computed distances are the right ones

TESTmean <- rowMeans(data[, 1:22])
TESTmean <- as.matrix(TESTmean)
TESTdistmean <- dist(TESTmean)
TESTdistmean <- as.matrix(TESTdistmean)

TESTsum <- rowSums(data[, 1:22])
TESTsum <- as.matrix(TESTsum)
TESTdistsum <- dist(TESTsum)
TESTdistsum <- as.matrix(TESTdistsum)

# computing euclidean distances between vectors:
testVector <- as.numeric(data[1,1:22])
testVector2 <- as.numeric(data[2,1:22])
dist(rbind(testVector, testVector2))

# checkup line about euclidean distance calculation
sqrt(sum((testVector-testVector2)^2)) # formula for euclidean distance between 2 points (here 2 vectors) in n-dimensional space (n being the number of indices in the vectors)



