setwd("../14_Info_timestep_Reconstruction/")
rm(list=ls())#empty the workspace

# Packages
source("fun.r")
library(ape)
library(phytools)
library(diagram)

# Read in tree and dataset
tree <- read.nexus("TreePicidae.txt") # the tree
data <- read.csv("Full_Species_data.csv", header=T, row.names = 1)

# reorder
data <- data[tree$tip.label,]
n <- Ntip(tree)
# Extract variable os interest to reconstruct: here drumming strategy is this variables (discrete variables with 6 states)
Strategies <- data$AcousticClust
names(Strategies) = tree$tip.label

##
Strategies[which(Strategies==1)] = "AC"
Strategies[which(Strategies==2)] = "DK"
Strategies[which(Strategies==3)] = "SF"
Strategies[which(Strategies==4)] = "SS"
Strategies[which(Strategies==5)] = "RS"
Strategies[which(Strategies==6)] = "IS"

# List of transition models (see the "fun.r" file)
modelList <- list(ER=ER_mod, SYM=SYM_mod, ARD=ARD_mod, Mod1=transMod1, Mod2=transMod2, Mod3=transMod3)

# Model fit using rerooting method
fit <- lapply(modelList,function(x){rerootingMethod(tree, Strategies, model=x)})


# Extract the AIC using aic() (see fun.r)
fit_aic <- sapply(1:length(fit),function(x){aicc(fit[[x]], n, modelList[[x]])})
colnames(fit_aic)=names(modelList)
fit_aic

# Extract the AIC weights (see fun)
fit_aicw <- aicw(unlist(fit_aic[1,]))

# The averaged transition matrix
Q <- Model_averaging(fit, fit_aicw)



## ------------------------------------------------------------ ##
##    Compute the transitions probabilities at the nodes        ##
## -------------------------------------------------------------##

# Here: reconstructing the transition probabilities at nodes based on (previously inferred) averaged transition matrix 
# Thus simulating discrete traits history based on this inferred matrix, to estimate transition frequencies at nodes 
# Note:it is also possible to fix this matrix using fitMk and estimate probabilities at the nodes too.

# Sample with the Q matrix
mtrees <- make.simmap(tree, Strategies, Q=Q, nsim=1000, pi="equal")

# Retrieve inferred nodes probabilities
YY <- describe.simmap(mtrees, all.equal.phylo=TRUE)
# cols <- setNames(c("green","lightblue","cornflowerblue", "orange", "darkgoldenrod", "lightcoral"),
#                 sort(unique(getStates(mtrees[[1]],"tips"))))

#edit colors to match drumming types
cols <- setNames(c("black", "red", "magenta", "cyan", "green", "blue"),
                 sort(unique(getStates(mtrees[[1]],"tips"))))

# Ploting of stochastic map
#pdf('stochastic_map.pdf')
plot(YY, fsize=0.4, cex=c(0.6,0.3), ftype="i", colors=cols, mar=c(2,2,2,2))
add.simmap.legend(rownames(Q), colors=cols, fsize=0.8, prompt=F, x=1, y=15)
axisPhylo()
#dev.off()


## ------------------------------------------------------------ ##
##                     Plot the transitions                     ##
## -------------------------------------------------------------##

M <- t(Q) # take the transpose with rerooting, but check for other methods

M <- M[order(factor(rownames(M), levels = c('IS', 'RS', 'DK', 'AC', 'SS', 'SF'))), 
       +     order(factor(colnames(M), levels = c('IS', 'RS', 'DK', 'AC', 'SS', 'SF')))]

diag(M) <- 0

# remove ~unsignificantly small?
M[M<5e-3]=0

# weight the lines
arrlwd <- 100*(M) # arbitrary scaling to display the links
arrlwd[arrlwd==-Inf]=0

#pdf('TransitionPlot_StochasticReconstruction.pdf')
pp <- plotmat(M, pos = NULL, curve=0.1, arr.lwd = arrlwd, name = rownames(M), 
              lwd = 1, box.lwd = 2,cex.txt = 0, box.size = 0.1, box.type = "ellipse", 
              box.prop = 0.5, box.cex=2,box.col=c('magenta','cyan','red','black','blue','green'),txt.col = 'white', txt.font = 2)
#dev.off()

## ------------------------------------------------------------ ##
##                     Generate Intermediate values             ##
## -------------------------------------------------------------##
# Get the duration in terms of branch length of the tree
source('treeLengths.R')
tLengths <- treeLengths(tree)

maxLen <- max(tLengths)
minLen <- min(tLengths)

# These should be equal (or very close)
print(c('Max Lenght=', maxLen))
print(c('Min Lenght=', minLen))

# Values to get number of branches and probs
tVals <- seq(from = maxLen*0.1, to = maxLen*.9, length.out = 20)

# Returns function with number of branches and composition at each time point
source('treeCutVals2.R')
xVals <- treeCutVals2(tree, tVals, YY, Strategies)


# Plot
cols <- setNames(c("black", "red", "magenta", "cyan", "green", "blue"),
                 sort(unique(getStates(mtrees[[1]],"tips"))))

plot(YY, fsize=0.4, cex=c(0.6,0.3), ftype="i", colors=cols, mar=c(2,2,2,2))
add.simmap.legend(rownames(Q), colors=cols, fsize=0.8, prompt=F, x=1, y=25)
for (i in 1:length(tVals)) {
  abline(v=tVals[i], lty = 'dashed')
}

save(data, tree, YY, cols, Q, tVals, xVals, file='stratMapp2.RData')

