library(MASS)
library(car)
library(ggplot2)
library(FactoMineR)
library(dplyr)
library(lme4)

setwd("../4_LI_Subset_acStrat/")

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

# Note that 6 componnents explain 75% of the variance; 10 would give 90%.

PCALoadings <- resPCA$var$cor[,1:6] #based on the PCs for which eigenvalues are > 1
# write.csv(PCALoadings, "PCALoadings.csv")

# adding PCs (with eigenval > 1) to dataframe
data$PC1<-resPCA$ind$coord[,1]
data$PC2<-resPCA$ind$coord[,2]
data$PC3<-resPCA$ind$coord[,3]
data$PC4<-resPCA$ind$coord[,4]
data$PC5<-resPCA$ind$coord[,5]
data$PC6<-resPCA$ind$coord[,6]


output_LIsubsetAcStrat <- matrix(0, nrow = 100, ncol = 8) # 8 columns,1 for overall classif %age, 1 for each acoustic strategy + 1 for full dataset; 1 row by iteration
colnames(output_LIsubsetAcStrat) <- c("Classification %age", "MeanLI_Overall","MeanLI_AcStrat1", "MeanLI_AcStrat2", "MeanLI_AcStrat3", "MeanLI_AcStrat4", "MeanLI_AcStrat5", "MeanLI_AcStrat6")
for (iteration in 1:100){

rm(list=setdiff(ls(), c("data","iteration","output_LIsubsetAcStrat","resPCA","PCALoadings")))


################### Run cross-validated DFA based on all files available (ranging from 3 to 30 files per species) 
# Edit for this part: need subset of same nb of species per acoustic strategies

AcStrat1 <- droplevels(data[data$AcousticClust==1,])
AcStrat2 <- droplevels(data[data$AcousticClust==2,])
AcStrat3 <- droplevels(data[data$AcousticClust==3,])
AcStrat4 <- droplevels(data[data$AcousticClust==4,])
AcStrat5 <- droplevels(data[data$AcousticClust==5,])
AcStrat6 <- droplevels(data[data$AcousticClust==6,])


randAcStrat1 <- sample(length(levels(AcStrat1$Species)),5, replace = F)
subAcStrat1 <- droplevels(AcStrat1[AcStrat1$Species %in% levels(AcStrat1$Species)[randAcStrat1],])

randAcStrat2 <- sample(length(levels(AcStrat2$Species)),5, replace = F)
subAcStrat2 <- droplevels(AcStrat2[AcStrat2$Species %in% levels(AcStrat2$Species)[randAcStrat2],])

randAcStrat3 <- sample(length(levels(AcStrat3$Species)),5, replace = F)
subAcStrat3 <- droplevels(AcStrat3[AcStrat3$Species %in% levels(AcStrat3$Species)[randAcStrat3],])

randAcStrat4 <- sample(length(levels(AcStrat4$Species)),5, replace = F)
subAcStrat4 <- droplevels(AcStrat4[AcStrat4$Species %in% levels(AcStrat4$Species)[randAcStrat4],])

randAcStrat5 <- sample(length(levels(AcStrat5$Species)),5, replace = F)
subAcStrat5 <- droplevels(AcStrat5[AcStrat5$Species %in% levels(AcStrat5$Species)[randAcStrat5],])

randAcStrat6 <- sample(length(levels(AcStrat6$Species)),5, replace = F)
subAcStrat6 <- droplevels(AcStrat6[AcStrat6$Species %in% levels(AcStrat6$Species)[randAcStrat6],])


subdata <- rbind(subAcStrat1,subAcStrat2,subAcStrat3,subAcStrat4,subAcStrat5,subAcStrat6)
#subdata2 <- subdata[order(as.character(subdata$Species)),]
#subdata2 <- subdata[order(row.names(subdata$Species)),]



## DFA with ALL data and CV
ndata=nrow(subdata)
species.names = levels(subdata$Species) #'Species' is the variable by which we want to discriminate drums
nbirds=length(species.names)
# There are nbirds = 92 species.

resLDA=lda(Species~PC1+PC2+PC3+PC4+PC5+PC6, prior=rep(1,nbirds)/nbirds, data=subdata, CV=T) # LDA computation. Important for information calculation: equal prior probabilities
tab <- table(subdata$Species, resLDA$class)
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
    ind.good <- (subdata$Species == species.names[i])
    cprob.tab.post[i,j] <- mean(resLDA$posterior[ind.good,j])
  }
}

# mean classification %age
(pcc <- mean(diag(cprob.tab.post)))
output_LIsubsetAcStrat[iteration,1] <- pcc


# classification %age at chance level
(chance=1/nbirds)


# ----- Define a function for plotting a matrix -----found here: http://www.phaget4.org/R/image_matrix.html #
myImagePlot <- function(x, ...){
  min <- min(x)
  max <- max(x)
  yLabels <- rownames(x)
  xLabels <- colnames(x)
  title <-c()
  
  # check for additional function arguments
  if( length(list(...)) ){
    Lst <- list(...)
    if( !is.null(Lst$zlim) ){
      min <- Lst$zlim[1]
      max <- Lst$zlim[2]
    }
    if( !is.null(Lst$yLabels) ){
      yLabels <- c(Lst$yLabels)
    }
    if( !is.null(Lst$xLabels) ){
      xLabels <- c(Lst$xLabels)
    }
    if( !is.null(Lst$title) ){
      title <- Lst$title
    }
  }
  # check for null values
  if( is.null(xLabels) ){
    xLabels <- c(1:ncol(x))
  }
  if( is.null(yLabels) ){
    yLabels <- c(1:nrow(x))
  }
  
  
  # Red, green and blue range from 0 to 1 (to get white background; if all 3 ranging from 1 to 0 then black background)
  ColorRamp <- rgb( seq(1,0,length=7),  # Red
                    seq(1,0,length=7),  # Green
                    seq(1,1,length=7))  # Blue
  ColorLevels <- seq(min, max, length=length(ColorRamp))
  
  #### example if only want greyscale:
  # ColorRamp <- rgb( seq(1,0,length=20),  # Red
  #                   seq(1,0,length=20),  # Green
  #                   seq(1,0,length=20))  # Blue
  # ColorLevels <- seq(min, max, length=length(ColorRamp))
  #### example if only want redscale:
  # ColorRamp <- rgb( seq(1,1,length=20),  # Red
  #                   seq(1,0,length=20),  # Green
  #                   seq(1,0,length=20))  # Blue
  # ColorLevels <- seq(min, max, length=length(ColorRamp))
  
  # Reverse Y axis
  # reverse <- nrow(x) : 1
  # yLabels <- yLabels[reverse]
  # x <- x[reverse,]
  
  # Data Map
  par(mar = c(5,5,5,2), cex=0.7, las = 2)
  image(1:length(xLabels), 1:length(yLabels), t(x), col=ColorRamp, xlab="",
        ylab="", axes=FALSE, zlim=c(min,max))
  if( !is.null(title) ){
    title(main=title,line = 3, cex.main=2)
  }
  axis(BELOW<-1, at=1:length(xLabels), labels=xLabels, cex.axis=0.7)
  axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, las= HORIZONTAL<-1,
       cex.axis=0.7)
  
}
#pdf("Confusion_matrix_DFA_CondProb.pdf")
#image(1:nbirds, 1:nbirds, -cprob.tab.post, xlim = c(1,nbirds), ylim=c(1,nbirds), col = terrain.colors(18), xlab='Predicted Species', ylab = 'Acutal Species')
myImagePlot(cprob.tab.post,title=sprintf('Posterior Confusion Matrix %.1f %%', pcc*100.0)) # Calling function just defined on validating dataset
#dev.off()


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

##################################### The following section is used to retrieve lda scores to later reconstruct ancestral states
# first this requires running the lda without CV to get LDs
resLDA_noCV=lda(Species~PC1+PC2+PC3+PC4+PC5+PC6, prior=rep(1,nbirds)/nbirds, data=subdata, CV=F) # LDA computation

# Get loadings of acoustic variables onto linear discriminants by multiplying two matrixes: scalings of PC1-PC6 onto the linear
# discriminants, and scalings of acoustic variables onto Principal Components
LD_loadings <- t(resLDA_noCV$scaling) %*% t(resPCA$var$cor[,1:6]) #express LDs as a function of acoustic variables

# write.csv(resLDA_noCV$scaling, "LD_Loadings_PCs.csv")
# write.csv(t(LD_loadings), "LD_Loadings_AcousticVars.csv")

# Create variables to stores LD scores
subdata$LD1 <- 0
subdata$LD2 <- 0
subdata$LD3 <- 0
subdata$LD4 <- 0
subdata$LD5 <- 0
subdata$LD6 <- 0

# Could be done with matrix multiply but more explicit like this
for (drum in 1:length(subdata[,1])) {
  subdata$LD1[drum] <- sum(subdata[drum, 17:38]*LD_loadings[1,])
  subdata$LD2[drum] <- sum(subdata[drum, 17:38]*LD_loadings[2,])
  subdata$LD3[drum] <- sum(subdata[drum, 17:38]*LD_loadings[3,])
  subdata$LD4[drum] <- sum(subdata[drum, 17:38]*LD_loadings[4,])
  subdata$LD5[drum] <- sum(subdata[drum, 17:38]*LD_loadings[5,])
  subdata$LD6[drum] <- sum(subdata[drum, 17:38]*LD_loadings[6,])
  }

##################################### End of linear discriminants computation

# Calculation of local information
# Go from conditional probability to joint probability distribution
normval_matrix <- cprob.tab.post/nbirds  

# Then retrieve the size of the matrix
matrixdimensions <- dim(normval_matrix)

#preparing matrix to store Localinfo
localInfoMatrix<-matrix(0,nrow = nbirds, ncol=1)

for(k in 1:matrixdimensions[1]){ # for each line
  infoLoc = 0
  for (kk in 1:matrixdimensions[2]){ #for each column 
    if (normval_matrix[k,kk] !=0){
      infoLoc <- infoLoc + nbirds*normval_matrix[k,kk]*(log2((nbirds*normval_matrix[k,kk])/(colSums(normval_matrix)[kk])))
    }
  }
  localInfoMatrix[k,1] = infoLoc
}

rownames(localInfoMatrix) <- rownames(normval_matrix)
localInfoMatrix



#####################################
(GlobalInfo <- mean(localInfoMatrix))
(mean(localInfoMatrix))
cat(c('global information from the classification matrix =',GlobalInfo,'bits'))
cat("Careful if the result is close to the ceiling value (here",log2(nbirds),"bits) ; result to take with caution ; similarly if the nb of analyzed signals is low.")
#####################################


# insert local info values into full dataset (careful, one value per species so need to apply same value to all drums within a species)
subdata$LI <- 0
subdata$LI[1] = localInfoMatrix[1]

l=1 #initiate counter
for (i in 2:length(subdata$Species)) 
{ if (subdata$Species[i]==subdata$Species[i-1]) {
  subdata$LI[i]=localInfoMatrix[l]
  }
  else {
    l=l+1
    subdata$LI[i]=localInfoMatrix[l]}
}

# similarly, retrieve specific classification %ages and insert them into full dataset
classif_species = 0
for (i in 1:nbirds){
  for (j in 1:nbirds){
    if (i == j){classif_species [i] <- cprob.tab.post[i,i]}
  }
}

subdata$SP_classif <- 0
subdata$SP_classif[1] = classif_species[1]
l=1 #initiate counter
for (i in 2:length(subdata$Species)) 
{ if (subdata$Species[i]==subdata$Species[i-1]) {subdata$SP_classif[i]=classif_species[l] }
  else {
    l=l+1
    subdata$SP_classif[i]=classif_species[l] }
}

# Also need to produce a species-specific medians table for further analyses (reconstructions, HierarchicalClustering, etc)
# initiate matrix (nrows = nb of species + 1 (header); ncol = Nb of variables + 1 (speciesName))
SpMatrix <- data.frame(matrix(0, nrow = nbirds, ncol = length(subdata[1,])))

firstlineofeachspecies <- subdata[ !duplicated(subdata$Species), ] # Retrieves 1st line of each species --> Do not use for numeric values (only factors)
SpMatrix[,1:16] <- firstlineofeachspecies[1:16]
colnames(SpMatrix)[1:16] <- colnames(subdata[,1:16])

# start a loop going through the variables of interest in the dataset (i.e. not the explanatory variables)
for (i in 17:length(subdata[1,])){
  
  # grouping the dataset by species levels, extract descriptives values for each variable
  # (with 'i' the column number of the variable of interest)  
  
  Test <-
    subdata %>%
    group_by(Species) %>%
    summarise_at(.vars = names(.)[i],funs(max, min, mean, median, sd))
  
  colnames(SpMatrix)[i] <- paste("median_", names(subdata[i]), sep='')
  SpMatrix[, i] <- Test$median
}


# Plot the species-specific normalized local information to show their distribution relative to ceiling
# pdf('SpeciesInfo_normed.pdf',width = 17, height = 7)
ggplot(SpMatrix, aes(x = Species_FullName, y= median_LI/log2(nbirds))) +
  geom_point(cex=2) +
  theme( panel.grid.minor = element_blank(),
         panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab('') +
  ylab('\nNormalized Mutual Information (bits)\n') +
  theme(axis.text.x=element_text(angle=50, size=11, hjust=1, face ="italic"),
        axis.text.y=element_text(angle=90, size=15, hjust=0.5),
        axis.title.y = element_text(colour="grey20",size=20,angle=90,hjust=.5,face="plain"))+
  geom_hline(yintercept=mean(localInfoMatrix/log2(nbirds)),color=2, lty=2)+
  geom_hline(yintercept=1,color=2, lty=1)
# dev.off()

mean(localInfoMatrix/log2(nbirds)) # to put in the MS
sd(localInfoMatrix/log2(nbirds)) # to put in the MS




##### Running stats comparing the mean Local information between acoustic strategies
### according stats to the graph (Anova, comparing between group means)
SpMatrix$AcousticClust <- as.factor(SpMatrix$AcousticClust)
SpMatrix$AcousticClust2 <- car::recode(SpMatrix$AcousticClust,"'1'='AC';'2'='DK';'3'='SF'; '4'='SS'; '5'='RS'; '6'='IS'") # car::recode written as is because conflicting with dplyr::recode
res.aov <- aov(median_LI/log2(nbirds) ~ AcousticClust2, data=SpMatrix)
summary(res.aov)
TukeyHSD(res.aov)




#### control lines before plotting
subdata1 <- SpMatrix[SpMatrix$AcousticClust==1,]
subdata2 <- SpMatrix[SpMatrix$AcousticClust==2,]
subdata3 <- SpMatrix[SpMatrix$AcousticClust==3,]
subdata4 <- SpMatrix[SpMatrix$AcousticClust==4,]
subdata5 <- SpMatrix[SpMatrix$AcousticClust==5,]
subdata6 <- SpMatrix[SpMatrix$AcousticClust==6,]

mean(SpMatrix$median_LI)/log2(nbirds)
mean(subdata1$median_LI)/log2(nbirds)
mean(subdata2$median_LI)/log2(nbirds)
mean(subdata3$median_LI)/log2(nbirds)
mean(subdata4$median_LI)/log2(nbirds)
mean(subdata5$median_LI)/log2(nbirds)
mean(subdata6$median_LI)/log2(nbirds)

output_LIsubsetAcStrat[iteration,2] <- mean(SpMatrix$median_LI)/log2(nbirds)
output_LIsubsetAcStrat[iteration,3] <- mean(subdata1$median_LI)/log2(nbirds)
output_LIsubsetAcStrat[iteration,4] <- mean(subdata2$median_LI)/log2(nbirds)
output_LIsubsetAcStrat[iteration,5] <- mean(subdata3$median_LI)/log2(nbirds)
output_LIsubsetAcStrat[iteration,6] <- mean(subdata4$median_LI)/log2(nbirds)
output_LIsubsetAcStrat[iteration,7] <- mean(subdata5$median_LI)/log2(nbirds)
output_LIsubsetAcStrat[iteration,8] <- mean(subdata6$median_LI)/log2(nbirds)

} # end of main loop selecting 5 species per acoustic strategy 100 times

write.csv(output_LIsubsetAcStrat, 'LI_byAcStrat_equalsizes.csv', row.names=F)

# testing whether the mean LI value are different across acoustic strategy when equal species sample considered (5 species maximum)
# thus 5 species randomly selected per acoustic strategy, after which calculation of information has been computed 100 times (script hereabove)
library(lsmeans)
data4stats=read.table(file="LI_byAcStrat_equalsizes_rearranged.csv", header=T, sep=",",dec=".")  # template from 100 simulations. Adjust according to new simulations
output_LIsubsetAcStrat=read.table(file="LI_byAcStrat_equalsizes.csv", header=T, sep=",",dec=".")  


res.anova <- aov(Mean_LI ~ Strategy, data=data4stats)
summary(res.anova)
TukeyHSD(res.anova)

# checkup lines
mean(output_LIsubsetAcStrat[,3]) # Strategy AC
mean(output_LIsubsetAcStrat[,4]) # Strategy DK
mean(output_LIsubsetAcStrat[,5]) # Strategy SF
mean(output_LIsubsetAcStrat[,6]) # Strategy SS
mean(output_LIsubsetAcStrat[,7]) # Strategy RS
mean(output_LIsubsetAcStrat[,8]) # Strategy IS


### Plot the acoustic strategies' mean local information and order graph based on above control lines
# Create function to get mean ± SE (otherwise boxplots only use quartiles and we want the graph to 
# match the stats, hence why plotting mean±SE with max/min, rather than median±Q25/Q75 with min/max)

MinMeanSEMMax <- function(x) {
  v <- c(min(x) , mean(x) - sd(x)/sqrt(length(x)), mean(x), mean(x) + sd(x)/sqrt(length(x)), max(x) )
  names(v) <- c("ymin", "lower", "middle", "upper", "ymax")
  v
}

# pdf("LI_byAcStrategies_equalSizes.pdf")
ggplot(data = data4stats, aes(reorder(factor(data4stats$Strategy), data4stats$Mean_LI*100, FUN=mean), data4stats$Mean_LI*100, fill = factor(data4stats$Strategy))) +
  stat_summary(fun.data=MinMeanSEMMax, geom="errorbar", position=position_dodge(1), width = .3, lwd=0.8) +
  stat_summary(fun.data=MinMeanSEMMax, geom="boxplot", position=position_dodge(1)) +
  scale_fill_manual(values=c(1:6)) +
  theme_bw() +
  ylab('Normalized Mutual Information (%) \n with equal number of species \n') +
  ylim(10, 100) +
  xlab('\nDrumming type') +
  scale_x_discrete(labels=c("SF\n(n = 5)", "SS\n(n = 5)", "DK\n(n = 5)","AC\n(n = 5)", "RS\n(n = 5)", "IS\n(n = 5)")) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_text(angle=45, size=15, vjust=0.5, face ="bold.italic"),
        axis.text.y=element_text(angle=90, size=15, hjust =0.5, vjust=1),
        axis.title.x = element_text(colour="grey20",size=20,angle=0,hjust=.5,face="plain"),
        axis.title.y = element_text(colour="grey20",size=20,angle=90,hjust=.5,face="plain"),
        legend.position="none")
# dev.off()    
