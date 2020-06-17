################################################################################
# LOAD PACKAGES
################################################################################

library(ape)
library(geiger)
library(nlme)
library(AICcmodavg)

################################################################################
# READ DATA & TREE
################################################################################

# READ DATA
setwd("../6_PGLS_Models/")
rm(list=ls())#empty the workspace
tree <- read.nexus("TreePicidae.txt") # the tree
tree$edge.length <- tree$edge.length/max(branching.times(tree)) # standardize the tree to unit height. By doing so, standardize the estimate of rates => easier for optimization and representation

data <- read.csv("Full_Species_data.csv", header=T, row.names = 1)

# reorder
data <- data[tree$tip.label,]

# #Check that names match
name.check(tree, data)
# subset only life history variables 
lifehistorydata <- as.data.frame(c(data[1], data[5], data[8], data[14:15], data[18]))
rownames(lifehistorydata) <- rownames(data)

######################## computing z-scores of life history variables (raw values until now)
# adjust depending on variables of interest
for (nbvariables in 2:6) # adjust depending on variables of interest
{
  meanoriginal=mean(lifehistorydata[[nbvariables]],na.rm=TRUE)
  sdoriginal=sd(lifehistorydata[[nbvariables]],na.rm=TRUE)
  for(i in 1:length(lifehistorydata[[nbvariables]])) # for each column (i.e. for each variable)
  {
    lifehistorydata[[i,nbvariables]]=(lifehistorydata[[i,nbvariables]]-meanoriginal)/sdoriginal
  }
}
# Replacing (if any) missing values by the mean z-score of the variable
for (nbvariables in 2:6) #adjust depending on variables of interest
{
  meanoriginal=mean(lifehistorydata[[nbvariables]],na.rm=TRUE)
  for(i in 1:length(lifehistorydata[[nbvariables]])) # for each column (i.e. for each variable)
  {
    if (is.na(lifehistorydata[[i,nbvariables]])==TRUE) lifehistorydata[[i,nbvariables]]=meanoriginal
  }
}


# Preparing traits: 
# structure
PC1 <- data$median_PC1
PC2 <- data$median_PC2
PC3 <- data$median_PC3
PC4 <- data$median_PC4
PC5 <- data$median_PC5
PC6 <- data$median_PC6

LD1 <- data$median_LD1
LD2 <- data$median_LD2
LD3 <- data$median_LD3
LD4 <- data$median_LD4
LD5 <- data$median_LD5
LD6 <- data$median_LD6

# function
localInfo <- data$median_LI/log2(92) #Using normalized LI

# life history:
bodyMass <- lifehistorydata$BodyMass #body Mass preferred, as used in Miles et al. (2018) and bodySize redundant with winglength
wingLength <- lifehistorydata$median_Wing_length 
beakWingRatio <- lifehistorydata$median_BeakLength_WingLength 

sympatry <- lifehistorydata$Overall_Sympatry
distribArea <- lifehistorydata$DistribArea_num


# attribute names as rownames for all variables created
names(PC1) = rownames(data)
names(PC2) = rownames(data)
names(PC3) = rownames(data)
names(PC4) = rownames(data)
names(PC5) = rownames(data)
names(PC6) = rownames(data)

names(LD1) = rownames(data)
names(LD2) = rownames(data)
names(LD3) = rownames(data)
names(LD4) = rownames(data)
names(LD5) = rownames(data)
names(LD6) = rownames(data)

names(localInfo) = rownames(data)

names(bodyMass) <- rownames(lifehistorydata)
names(wingLength) <- rownames(lifehistorydata)
names(beakWingRatio) <- rownames(lifehistorydata)

names(sympatry) <- rownames(lifehistorydata)
names(distribArea) <- rownames(lifehistorydata)

######################### Looking at models of PC1 (PGLS Pagel): ######################### 
#NULL model
PC1_pglsPagel0 <- gls(PC1 ~ 1, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(PC1_pglsPagel0)
AICc(PC1_pglsPagel0,second.ord = TRUE)

# Models with single variables (for stepwise forward selection)
PC1_pglsPagel1 <- gls(PC1 ~ bodyMass, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(PC1_pglsPagel1)
AICc(PC1_pglsPagel1,second.ord = TRUE)

PC1_pglsPagel2 <- gls(PC1 ~ wingLength, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(PC1_pglsPagel2)
AICc(PC1_pglsPagel2,second.ord = TRUE)

PC1_pglsPagel3 <- gls(PC1 ~ beakWingRatio, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(PC1_pglsPagel3)
AICc(PC1_pglsPagel3,second.ord = TRUE)

PC1_pglsPagel4 <- gls(PC1 ~ sympatry, correlation = corPagel(1, phy=tree, fixed=F), data=data, method="REML")
summary(PC1_pglsPagel4)
AICc(PC1_pglsPagel4,second.ord = TRUE)

PC1_pglsPagel5 <- gls(PC1 ~ distribArea, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(PC1_pglsPagel5)
AICc(PC1_pglsPagel5,second.ord = TRUE)
 

######################### Looking at models of PC2 (PGLS Pagel): ######################### 
#NULL model
PC2_pglsPagel0 <- gls(PC2 ~ 1, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(PC2_pglsPagel0)
AICc(PC2_pglsPagel0,second.ord = TRUE)

# Models with single variables (for stepwise forward selection)
PC2_pglsPagel1 <- gls(PC2 ~ bodyMass, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(PC2_pglsPagel1)
AICc(PC2_pglsPagel1,second.ord = TRUE)

PC2_pglsPagel2 <- gls(PC2 ~ wingLength, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(PC2_pglsPagel2)
AICc(PC2_pglsPagel2,second.ord = TRUE)

PC2_pglsPagel3 <- gls(PC2 ~ beakWingRatio, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(PC2_pglsPagel3)
AICc(PC2_pglsPagel3,second.ord = TRUE)

PC2_pglsPagel4 <- gls(PC2 ~ sympatry, correlation = corPagel(1, phy=tree, fixed=F), data=data, method="REML")
summary(PC2_pglsPagel4)
AICc(PC2_pglsPagel4,second.ord = TRUE)

PC2_pglsPagel5 <- gls(PC2 ~ distribArea, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(PC2_pglsPagel5)
AICc(PC2_pglsPagel5,second.ord = TRUE)


######################### Looking at models of PC3 (PGLS Pagel): ######################### 
#NULL model
PC3_pglsPagel0 <- gls(PC3 ~ 1, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(PC3_pglsPagel0)
AICc(PC3_pglsPagel0,second.ord = TRUE)

# Models with single variables (for stepwise forward selection)
PC3_pglsPagel1 <- gls(PC3 ~ bodyMass, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(PC3_pglsPagel1)
AICc(PC3_pglsPagel1,second.ord = TRUE)

PC3_pglsPagel2 <- gls(PC3 ~ wingLength, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(PC3_pglsPagel2)
AICc(PC3_pglsPagel2,second.ord = TRUE)

PC3_pglsPagel3 <- gls(PC3 ~ beakWingRatio, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(PC3_pglsPagel3)
AICc(PC3_pglsPagel3,second.ord = TRUE)

PC3_pglsPagel4 <- gls(PC3 ~ sympatry, correlation = corPagel(1, phy=tree, fixed=F), data=data, method="REML")
summary(PC3_pglsPagel4)
AICc(PC3_pglsPagel4,second.ord = TRUE)

PC3_pglsPagel5 <- gls(PC3 ~ distribArea, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(PC3_pglsPagel5)
AICc(PC3_pglsPagel5,second.ord = TRUE)


######################### Looking at models of PC4 (PGLS Pagel): ######################### 
#NULL model
PC4_pglsPagel0 <- gls(PC4 ~ 1, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(PC4_pglsPagel0)
AICc(PC4_pglsPagel0,second.ord = TRUE)

# Models with single variables (for stepwise forward selection)
PC4_pglsPagel1 <- gls(PC4 ~ bodyMass, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(PC4_pglsPagel1)
AICc(PC4_pglsPagel1,second.ord = TRUE)

PC4_pglsPagel2 <- gls(PC4 ~ wingLength, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(PC4_pglsPagel2)
AICc(PC4_pglsPagel2,second.ord = TRUE)

PC4_pglsPagel3 <- gls(PC4 ~ beakWingRatio, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(PC4_pglsPagel3)
AICc(PC4_pglsPagel3,second.ord = TRUE)

PC4_pglsPagel4 <- gls(PC4 ~ sympatry, correlation = corPagel(1, phy=tree, fixed=F), data=data, method="REML")
summary(PC4_pglsPagel4)
AICc(PC4_pglsPagel4,second.ord = TRUE)

PC4_pglsPagel5 <- gls(PC4 ~ distribArea, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(PC4_pglsPagel5)
AICc(PC4_pglsPagel5,second.ord = TRUE)


######################### Looking at models of PC5 (PGLS Pagel): ######################### 
#NULL model
PC5_pglsPagel0 <- gls(PC5 ~ 1, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(PC5_pglsPagel0)
AICc(PC5_pglsPagel0,second.ord = TRUE)

# Models with single variables (for stepwise forward selection)
PC5_pglsPagel1 <- gls(PC5 ~ bodyMass, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(PC5_pglsPagel1)
AICc(PC5_pglsPagel1,second.ord = TRUE)

PC5_pglsPagel2 <- gls(PC5 ~ wingLength, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(PC5_pglsPagel2)
AICc(PC5_pglsPagel2,second.ord = TRUE)

PC5_pglsPagel3 <- gls(PC5 ~ beakWingRatio, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(PC5_pglsPagel3)
AICc(PC5_pglsPagel3,second.ord = TRUE)

PC5_pglsPagel4 <- gls(PC5 ~ sympatry, correlation = corPagel(1, phy=tree, fixed=F), data=data, method="REML")
summary(PC5_pglsPagel4)
AICc(PC5_pglsPagel4,second.ord = TRUE)

PC5_pglsPagel5 <- gls(PC5 ~ distribArea, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(PC5_pglsPagel5)
AICc(PC5_pglsPagel5,second.ord = TRUE)


######################### Looking at models of PC6 (PGLS Pagel): ######################### 
#NULL model
PC6_pglsPagel0 <- gls(PC6 ~ 1, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(PC6_pglsPagel0)
AICc(PC6_pglsPagel0,second.ord = TRUE)

# Models with single variables (for stepwise forward selection)
PC6_pglsPagel1 <- gls(PC6 ~ bodyMass, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(PC6_pglsPagel1)
AICc(PC6_pglsPagel1,second.ord = TRUE)

PC6_pglsPagel2 <- gls(PC6 ~ wingLength, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(PC6_pglsPagel2)
AICc(PC6_pglsPagel2,second.ord = TRUE)

PC6_pglsPagel3 <- gls(PC6 ~ beakWingRatio, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(PC6_pglsPagel3)
AICc(PC6_pglsPagel3,second.ord = TRUE)

PC6_pglsPagel4 <- gls(PC6 ~ sympatry, correlation = corPagel(1, phy=tree, fixed=F), data=data, method="REML")
summary(PC6_pglsPagel4)
AICc(PC6_pglsPagel4,second.ord = TRUE)

PC6_pglsPagel5 <- gls(PC6 ~ distribArea, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(PC6_pglsPagel5)
AICc(PC6_pglsPagel5,second.ord = TRUE)

# Models with two variables (for stepwise forward selection), since PC6_pglsPagel1 and PC6_pglsPagel2 lead to
# a significant decrease of AICc compared with the NULL model
PC6_pglsPagel6 <- gls(PC6 ~ bodyMass+wingLength, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(PC6_pglsPagel6)
AICc(PC6_pglsPagel6,second.ord = TRUE)

PC6_pglsPagel7 <- gls(PC6 ~ bodyMass+beakWingRatio, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(PC6_pglsPagel7)
AICc(PC6_pglsPagel7,second.ord = TRUE)

PC6_pglsPagel8 <- gls(PC6 ~ bodyMass+sympatry, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(PC6_pglsPagel8)
AICc(PC6_pglsPagel8,second.ord = TRUE)

PC6_pglsPagel9 <- gls(PC6 ~ bodyMass+distribArea, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(PC6_pglsPagel9)
AICc(PC6_pglsPagel9,second.ord = TRUE)

PC6_pglsPagel10 <- gls(PC6 ~ wingLength+beakWingRatio, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(PC6_pglsPagel10)
AICc(PC6_pglsPagel10,second.ord = TRUE)

PC6_pglsPagel11 <- gls(PC6 ~ wingLength+sympatry, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(PC6_pglsPagel11)
AICc(PC6_pglsPagel11,second.ord = TRUE)

PC6_pglsPagel12 <- gls(PC6 ~ wingLength+distribArea, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(PC6_pglsPagel12)
AICc(PC6_pglsPagel12,second.ord = TRUE)

## At the end of this section, only two models improve the AICc compared with the NULL model for PC6 only (No improvement for
## the models using PC1-PC5): PC6_pglsPagel1 and PC6_pglsPagel2, although PC6_pglsPagel2 only very slightly

## We can now test whether the effect of body Mass (in model PC6_pglsPagel1) or Wing Length (in model PC6_pglsPagel2) is marginal 
## or really important in the model by doing an anova between the full and the reduced models:
# To do this, we fix the phylogenetic signal obtained in REML models

Fixed_PC6_pglsPagel1<- gls(PC6 ~ bodyMass, correlation = corPagel(0.8457416, phy=tree, fixed=TRUE), data=data, method="ML")
Fixed_PC6_pglsPagel0 <- gls(PC6 ~ 1, correlation = corPagel(0.9060799, phy=tree, fixed=TRUE), data=data, method="ML")
anova(Fixed_PC6_pglsPagel1,Fixed_PC6_pglsPagel0)

Fixed_PC6_pglsPagel2<- gls(PC6 ~ wingLength, correlation = corPagel(0.8442068, phy=tree, fixed=TRUE), data=data, method="ML")
Fixed_PC6_pglsPagel0 <- gls(PC6 ~ 1, correlation = corPagel(0.9060799, phy=tree, fixed=TRUE), data=data, method="ML")
anova(Fixed_PC6_pglsPagel2,Fixed_PC6_pglsPagel0)

#######################################################################################
#######################################################################################

################ REPEATING THE APPROACH ABOVE USING LDs INSTEAD OF PCs ################

#######################################################################################
#######################################################################################


######################### Looking at models of LD1 (PGLS Pagel): ######################### 
#NULL model
LD1_pglsPagel0 <- gls(LD1 ~ 1, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(LD1_pglsPagel0)
AICc(LD1_pglsPagel0,second.ord = TRUE)

# Models with single variables (for stepwise forward selection)
LD1_pglsPagel1 <- gls(LD1 ~ bodyMass, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(LD1_pglsPagel1)
AICc(LD1_pglsPagel1,second.ord = TRUE)

LD1_pglsPagel2 <- gls(LD1 ~ wingLength, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(LD1_pglsPagel2)
AICc(LD1_pglsPagel2,second.ord = TRUE)

LD1_pglsPagel3 <- gls(LD1 ~ beakWingRatio, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(LD1_pglsPagel3)
AICc(LD1_pglsPagel3,second.ord = TRUE)

LD1_pglsPagel4 <- gls(LD1 ~ sympatry, correlation = corPagel(1, phy=tree, fixed=F), data=data, method="REML")
summary(LD1_pglsPagel4)
AICc(LD1_pglsPagel4,second.ord = TRUE)

LD1_pglsPagel5 <- gls(LD1 ~ distribArea, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(LD1_pglsPagel5)
AICc(LD1_pglsPagel5,second.ord = TRUE)


######################### Looking at models of LD2 (PGLS Pagel): ######################### 
#NULL model
LD2_pglsPagel0 <- gls(LD2 ~ 1, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(LD2_pglsPagel0)
AICc(LD2_pglsPagel0,second.ord = TRUE)

# Models with single variables (for stepwise forward selection)
LD2_pglsPagel1 <- gls(LD2 ~ bodyMass, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(LD2_pglsPagel1)
AICc(LD2_pglsPagel1,second.ord = TRUE)

LD2_pglsPagel2 <- gls(LD2 ~ wingLength, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(LD2_pglsPagel2)
AICc(LD2_pglsPagel2,second.ord = TRUE)

LD2_pglsPagel3 <- gls(LD2 ~ beakWingRatio, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(LD2_pglsPagel3)
AICc(LD2_pglsPagel3,second.ord = TRUE)

LD2_pglsPagel4 <- gls(LD2 ~ sympatry, correlation = corPagel(1, phy=tree, fixed=F), data=data, method="REML")
summary(LD2_pglsPagel4)
AICc(LD2_pglsPagel4,second.ord = TRUE)

LD2_pglsPagel5 <- gls(LD2 ~ distribArea, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(LD2_pglsPagel5)
AICc(LD2_pglsPagel5,second.ord = TRUE)


######################### Looking at models of LD3 (PGLS Pagel): ######################### 
#NULL model
LD3_pglsPagel0 <- gls(LD3 ~ 1, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(LD3_pglsPagel0)
AICc(LD3_pglsPagel0,second.ord = TRUE)

# Models with single variables (for stepwise forward selection)
LD3_pglsPagel1 <- gls(LD3 ~ bodyMass, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(LD3_pglsPagel1)
AICc(LD3_pglsPagel1,second.ord = TRUE)

LD3_pglsPagel2 <- gls(LD3 ~ wingLength, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(LD3_pglsPagel2)
AICc(LD3_pglsPagel2,second.ord = TRUE)

LD3_pglsPagel3 <- gls(LD3 ~ beakWingRatio, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(LD3_pglsPagel3)
AICc(LD3_pglsPagel3,second.ord = TRUE)

LD3_pglsPagel4 <- gls(LD3 ~ sympatry, correlation = corPagel(1, phy=tree, fixed=F), data=data, method="REML")
summary(LD3_pglsPagel4)
AICc(LD3_pglsPagel4,second.ord = TRUE)

LD3_pglsPagel5 <- gls(LD3 ~ distribArea, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(LD3_pglsPagel5)
AICc(LD3_pglsPagel5,second.ord = TRUE)


######################### Looking at models of LD4 (PGLS Pagel): ######################### 
#NULL model
LD4_pglsPagel0 <- gls(LD4 ~ 1, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(LD4_pglsPagel0)
AICc(LD4_pglsPagel0,second.ord = TRUE)

# Models with single variables (for stepwise forward selection)
LD4_pglsPagel1 <- gls(LD4 ~ bodyMass, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(LD4_pglsPagel1)
AICc(LD4_pglsPagel1,second.ord = TRUE)

LD4_pglsPagel2 <- gls(LD4 ~ wingLength, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(LD4_pglsPagel2)
AICc(LD4_pglsPagel2,second.ord = TRUE)

LD4_pglsPagel3 <- gls(LD4 ~ beakWingRatio, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(LD4_pglsPagel3)
AICc(LD4_pglsPagel3,second.ord = TRUE)

LD4_pglsPagel4 <- gls(LD4 ~ sympatry, correlation = corPagel(1, phy=tree, fixed=F), data=data, method="REML")
summary(LD4_pglsPagel4)
AICc(LD4_pglsPagel4,second.ord = TRUE)

LD4_pglsPagel5 <- gls(LD4 ~ distribArea, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(LD4_pglsPagel5)
AICc(LD4_pglsPagel5,second.ord = TRUE)


######################### Looking at models of LD5 (PGLS Pagel): ######################### 
#NULL model
LD5_pglsPagel0 <- gls(LD5 ~ 1, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(LD5_pglsPagel0)
AICc(LD5_pglsPagel0,second.ord = TRUE)

# Models with single variables (for stepwise forward selection)
LD5_pglsPagel1 <- gls(LD5 ~ bodyMass, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(LD5_pglsPagel1)
AICc(LD5_pglsPagel1,second.ord = TRUE)

LD5_pglsPagel2 <- gls(LD5 ~ wingLength, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(LD5_pglsPagel2)
AICc(LD5_pglsPagel2,second.ord = TRUE)

LD5_pglsPagel3 <- gls(LD5 ~ beakWingRatio, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(LD5_pglsPagel3)
AICc(LD5_pglsPagel3,second.ord = TRUE)

LD5_pglsPagel4 <- gls(LD5 ~ sympatry, correlation = corPagel(1, phy=tree, fixed=F), data=data, method="REML")
summary(LD5_pglsPagel4)
AICc(LD5_pglsPagel4,second.ord = TRUE)

LD5_pglsPagel5 <- gls(LD5 ~ distribArea, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(LD5_pglsPagel5)
AICc(LD5_pglsPagel5,second.ord = TRUE)


######################### Looking at models of LD6 (PGLS Pagel): ######################### 
#NULL model
LD6_pglsPagel0 <- gls(LD6 ~ 1, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(LD6_pglsPagel0)
AICc(LD6_pglsPagel0,second.ord = TRUE)

# Models with single variables (for stepwise forward selection)
LD6_pglsPagel1 <- gls(LD6 ~ bodyMass, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(LD6_pglsPagel1)
AICc(LD6_pglsPagel1,second.ord = TRUE)

LD6_pglsPagel2 <- gls(LD6 ~ wingLength, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(LD6_pglsPagel2)
AICc(LD6_pglsPagel2,second.ord = TRUE)

LD6_pglsPagel3 <- gls(LD6 ~ beakWingRatio, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(LD6_pglsPagel3)
AICc(LD6_pglsPagel3,second.ord = TRUE)

LD6_pglsPagel4 <- gls(LD6 ~ sympatry, correlation = corPagel(1, phy=tree, fixed=F), data=data, method="REML")
summary(LD6_pglsPagel4)
AICc(LD6_pglsPagel4,second.ord = TRUE)

LD6_pglsPagel5 <- gls(LD6 ~ distribArea, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(LD6_pglsPagel5)
AICc(LD6_pglsPagel5,second.ord = TRUE)

## At the end of this section, NO model improves the AICc compared with the NULL model for any of the LD1-LD6 considered


#######################################################################################
#######################################################################################

############# REPEATING THE APPROACH ABOVE for LI INSTEAD OF PCs OR LDs ###############

#######################################################################################
#######################################################################################

######################### Looking at models of Local Info (PGLS Pagel): ######################### 
#NULL model
localInfo_pglsPagel0 <- gls(localInfo ~ 1, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(localInfo_pglsPagel0)
AICc(localInfo_pglsPagel0,second.ord = TRUE)

# Models with single variables (for stepwise forward selection)
localInfo_pglsPagel1 <- gls(localInfo ~ bodyMass, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(localInfo_pglsPagel1)
AICc(localInfo_pglsPagel1,second.ord = TRUE)

localInfo_pglsPagel2 <- gls(localInfo ~ wingLength, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(localInfo_pglsPagel2)
AICc(localInfo_pglsPagel2,second.ord = TRUE)

localInfo_pglsPagel3 <- gls(localInfo ~ beakWingRatio, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(localInfo_pglsPagel3)
AICc(localInfo_pglsPagel3,second.ord = TRUE)

localInfo_pglsPagel4 <- gls(localInfo ~ sympatry, correlation = corPagel(1, phy=tree, fixed=F), data=data, method="REML")
summary(localInfo_pglsPagel4)
AICc(localInfo_pglsPagel4,second.ord = TRUE)

localInfo_pglsPagel5 <- gls(localInfo ~ distribArea, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(localInfo_pglsPagel5)
AICc(localInfo_pglsPagel5,second.ord = TRUE)

### In addition to testing life history variable, testing whether struture is correlated with LI (PCs and LDs)

localInfo_pglsPagel6 <- gls(localInfo ~ PC1, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(localInfo_pglsPagel6)
AICc(localInfo_pglsPagel6,second.ord = TRUE)

localInfo_pglsPagel7 <- gls(localInfo ~ PC2, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(localInfo_pglsPagel7)
AICc(localInfo_pglsPagel7,second.ord = TRUE)

localInfo_pglsPagel8 <- gls(localInfo ~ PC3, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(localInfo_pglsPagel8)
AICc(localInfo_pglsPagel8,second.ord = TRUE)

localInfo_pglsPagel9 <- gls(localInfo ~ PC4, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(localInfo_pglsPagel9)
AICc(localInfo_pglsPagel9,second.ord = TRUE)

localInfo_pglsPagel10 <- gls(localInfo ~ PC5, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(localInfo_pglsPagel10)
AICc(localInfo_pglsPagel10,second.ord = TRUE)

localInfo_pglsPagel11 <- gls(localInfo ~ PC6, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(localInfo_pglsPagel11)
AICc(localInfo_pglsPagel11,second.ord = TRUE)



localInfo_pglsPagel12 <- gls(localInfo ~ LD1, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(localInfo_pglsPagel12)
AICc(localInfo_pglsPagel12,second.ord = TRUE)

localInfo_pglsPagel13 <- gls(localInfo ~ LD2, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(localInfo_pglsPagel13)
AICc(localInfo_pglsPagel13,second.ord = TRUE)

localInfo_pglsPagel14 <- gls(localInfo ~ LD3, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(localInfo_pglsPagel14)
AICc(localInfo_pglsPagel14,second.ord = TRUE)

localInfo_pglsPagel15 <- gls(localInfo ~ LD4, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(localInfo_pglsPagel15)
AICc(localInfo_pglsPagel15,second.ord = TRUE)

localInfo_pglsPagel16 <- gls(localInfo ~ LD5, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(localInfo_pglsPagel16)
AICc(localInfo_pglsPagel16,second.ord = TRUE)

localInfo_pglsPagel17 <- gls(localInfo ~ LD6, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(localInfo_pglsPagel17)
AICc(localInfo_pglsPagel17,second.ord = TRUE)

# Models with two variables (for stepwise forward selection), since localInfo_pglsPagel6 and localInfo_pglsPagel12 lead to
# a significant decrease of AICc compared with the NULL model

localInfo_pglsPagel18 <- gls(localInfo ~ PC1+PC2, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(localInfo_pglsPagel18)
AICc(localInfo_pglsPagel18,second.ord = TRUE)

localInfo_pglsPagel19 <- gls(localInfo ~ PC1+PC3, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(localInfo_pglsPagel19)
AICc(localInfo_pglsPagel19,second.ord = TRUE)

localInfo_pglsPagel20 <- gls(localInfo ~ PC1+PC4, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(localInfo_pglsPagel20)
AICc(localInfo_pglsPagel20,second.ord = TRUE)

localInfo_pglsPagel21 <- gls(localInfo ~ PC1+PC5, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(localInfo_pglsPagel21)
AICc(localInfo_pglsPagel21,second.ord = TRUE)

localInfo_pglsPagel22 <- gls(localInfo ~ PC1+PC6, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(localInfo_pglsPagel22)
AICc(localInfo_pglsPagel22,second.ord = TRUE)



localInfo_pglsPagel23 <- gls(localInfo ~ LD1+LD2, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(localInfo_pglsPagel23)
AICc(localInfo_pglsPagel23,second.ord = TRUE)

localInfo_pglsPagel24 <- gls(localInfo ~ LD1+LD3, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(localInfo_pglsPagel24)
AICc(localInfo_pglsPagel24,second.ord = TRUE)

localInfo_pglsPagel25 <- gls(localInfo ~ LD1+LD4, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(localInfo_pglsPagel25)
AICc(localInfo_pglsPagel25,second.ord = TRUE)

localInfo_pglsPagel26 <- gls(localInfo ~ LD1+LD5, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(localInfo_pglsPagel26)
AICc(localInfo_pglsPagel26,second.ord = TRUE)

localInfo_pglsPagel27 <- gls(localInfo ~ LD1+LD6, correlation = corPagel(1, phy=tree, fixed=FALSE), data=data, method="REML")
summary(localInfo_pglsPagel27)
AICc(localInfo_pglsPagel27,second.ord = TRUE)

## At the end of this section, only two models improve the AICc compared with the NULL model: localInfo_pglsPagel6 and localInfo_pglsPagel12

## We can now test whether the effect of PC1s (in model localInfo_pglsPagel6) or LD1 (in model localInfo_pglsPagel12) is marginal 
## or really important in the model by doing an anova between the full and the reduced models:
# To do this, we fix the phylogenetic signal obtained in REML models


Fixed_localInfo_pglsPagel6 <- gls(localInfo ~ PC1, correlation = corPagel(0.9390454, phy=tree, fixed=TRUE), data=data, method="ML")
Fixed_localInfo_pglsPagel0 <- gls(localInfo ~ 1, correlation = corPagel(0.9456714, phy=tree, fixed=TRUE), data=data, method="ML")
anova(Fixed_localInfo_pglsPagel6,Fixed_localInfo_pglsPagel0)

Fixed_localInfo_pglsPagel12 <- gls(localInfo ~ LD1, correlation = corPagel(0.8888651, phy=tree, fixed=TRUE), data=data, method="ML")
Fixed_localInfo_pglsPagel0 <- gls(localInfo ~ 1, correlation = corPagel(0.9456714, phy=tree, fixed=TRUE), data=data, method="ML")
anova(Fixed_localInfo_pglsPagel12,Fixed_localInfo_pglsPagel0)



