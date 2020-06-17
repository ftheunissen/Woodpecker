#Load relevant packages
library(ape)
library(geiger)


################################################################################
# READ DATA & TREES
################################################################################
#Change working directory to folder with your files (change accordingly)
rm(list=ls())
setwd("../7_EvolRate_Models/")
source('../7_EvolRate_Models/ratesDependent_gen.r')


tree <- read.nexus("TreePicidae.txt") # the tree
tree$edge.length <- tree$edge.length/max(branching.times(tree)) # I standardize the tree t

data <- read.csv("Acoustic_species_postDFA.csv", header=T, row.names = 1)

# reorder
data <- data[tree$tip.label,]

# #Check that names match
name.check(tree, data)

# Preparing traits: 
# structure
PC1 <- data$median_PC1/sd(data$median_PC1)
PC2 <- data$median_PC2/sd(data$median_PC2)
PC3 <- data$median_PC3/sd(data$median_PC3)
PC4 <- data$median_PC4/sd(data$median_PC4)
PC5 <- data$median_PC5/sd(data$median_PC5)
PC6 <- data$median_PC6/sd(data$median_PC6)

LD1 <- data$median_LD1/sd(data$median_LD1)
LD2 <- data$median_LD2/sd(data$median_LD2)
LD3 <- data$median_LD3/sd(data$median_LD3)
LD4 <- data$median_LD4/sd(data$median_LD4)
LD5 <- data$median_LD5/sd(data$median_LD5)
LD6 <- data$median_LD6/sd(data$median_LD6)

# function
normLI <- data$median_LI/log2(92) #Using normalized LI
localInfo <- normLI/sd(normLI)

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

# Bayesian analysis using Monte-Carlo's method by Markov chains

dependentVariables <- list(PC1, PC2, PC3, PC4, PC5, PC6, LD1, LD2, LD3, LD4, LD5, LD6, localInfo)
names(dependentVariables) <- c('PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'LD1', 'LD2', 'LD3', 'LD4', 'LD5', 'LD6', 'localInfo')
# dependentVariables[1]  # check line
# names(dependentVariables[1]) # check line


for (var in 1:length(dependentVariables)){
  # MCMC analysis: use reversible jump mcmc sampling.
  # number of iterations for the mcmc chain
  niter <- 1e6
  # sampling the chain every:
  sampling <- 1000
  # burnin (proportion of the chain to reject as burnin) 
  burnin <- 0.25
  # name for the analysis folder
  folderName <- names(dependentVariables[var])
  
  # Number of branches in the tree
  nn <- 2*Ntip(tree)-2
  ## PRIORS
  # Half-Cauchy distribution for the variance/standard error (see also Gelman 2006)
  dhalfCauchy <- function(x, scale=1, log=F){ 
    if(any(x<0)) return(-Inf)
    density <- 2 /(pi * (1+(x/scale)^2))
    if(log==TRUE) density<-log(density/scale)
    return(density/scale)
  }
  # Priors
  ratePrior <- function(x) dhalfCauchy(sqrt(x), 25, log=TRUE) 
  sePrior <- function(x) dhalfCauchy(sqrt(x), 25, log=TRUE) 
  rootPrior <- function(x) dnorm(x, sd=10, mean=0, log=TRUE) 
  shiftPrior <- dcount(0:(nn-1), FUN=dpois, lambda=log(2))

# running models with relaxed brownian models.

  rjmcmc.bm(tree, dependentVariables[[var]], prop.width=4, ngen=niter, samp=sampling, filebase=folderName, 
            simple.start=TRUE, type="rbm",  ## adjust model type here
            dlnRATE=ratePrior,
            #dlnROOT=rootPrior,
            dlnSE=sePrior, 
            dlnSHIFT=shiftPrior)  
  
  # Load the results ### Adjust according to model type
  outdir <- paste("relaxedBM", folderName, sep=".")
  ps <- load.rjmcmc(outdir)
  
  # Check the traces of the mcmc run:
  require(coda)
  
  ## Loading required package: coda
  range_burnin <- (niter/sampling*burnin):(niter/sampling) 
  
  #Plots
  # Plot the shifts and their probabilities on the tree, as well as traces for the simulation/model
  pdf(paste0('Bayesian_',names(dependentVariables[var]),'_relaxedBM.pdf'))
  plotval <- plot(x=ps, par="shifts", burnin=burnin, 
                  legend=TRUE, show.tip=TRUE, edge.width=2)

  plot(mcmc(ps$log[range_burnin,c("min","median","max")]))
  plot(mcmc(ps$log[range_burnin,c("shifts","lnL")]))
  plot(mcmc(ps$log[range_burnin,c("SE","root","ppos")]))
  dev.off()
}


##################################
##################################

# If want to read in a folder already obtained from a run, while plotting clearer shifts (with a black circle), run the part below (for some
# reason, modifying the source code implies running everything from scratch every time)

### First run lines 1-70

var <- 13 # replace variable number as needed
niter <- 1e6
sampling <- 1000
burnin <- 0.25
folderName <- names(dependentVariables[var])
outdir <- paste("relaxedBM", folderName, sep=".")
outdir #control line for what folder will be loaded
ps <- load.rjmcmc(outdir) #folderName ## make sure folder is in correct working dir

# Check the traces of the mcmc run:
require(coda)
## Loading required package: coda
range_burnin <- (niter/sampling*burnin):(niter/sampling) 
plot(mcmc(ps$log[range_burnin,c("min","median","max")]))
plot(mcmc(ps$log[range_burnin,c("shifts","lnL")]))
plot(mcmc(ps$log[range_burnin,c("SE","root","ppos")]))

source('../7_EvolRate_Models/Geiger_util.r')
source('../7_EvolRate_Models/rjmcmcPlot.r')

#Plots
# Plot the shifts and their probabilities on the tree, as well as traces for the simulation/model
pdf(paste0('Bayesian_',names(dependentVariables[var]),'relaxedBM.pdf'))
plotval <- plot(x=ps, par="shifts", burnin=burnin, label.offset = 0.02,
                legend=TRUE, show.tip=TRUE, edge.width=2, x.lim=5,adj=0.5)
dev.off()

###### edit for correlations between shifts/rates 
write.csv(plotval$median.rates,paste0('rates_',names(dependentVariables[var]),'.csv'))
###### end of edit for correlations


## note: the two source files above are needed as the source code of the plot.rjmcmc has been changed.
## in particular, the function edgelabels.auteur is modified so that shift probabilities are not displayed on the tree anymore:
# original: edgelabels.auteur(text=NULL, pch=ifelse(ll==1, 21, NA), cex=4*cc, col=.transparency(rr,0.95), bg=.transparency(rr,0.5), lwd=0.5)
# modified: edgelabels.auteur(text=NULL, pch=21, cex=4*cc, col=.transparency("black",1), bg=.transparency(rr,0.5), lwd=0.5)
# in addition, an edit regarding color scaling ahs been made, line 564 of rjmcmc.plot code:
#original:   cce <- colorspace::diverge_hcl(2 * colors$branches + 1, power = 0.5);
#modified:   cce <- colorspace::sequential_hcl(2 * colors$branches + 1, rev = T, h = 0, power = 2);







