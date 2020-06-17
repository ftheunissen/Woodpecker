#Change working directory to folder with your files (change accordingly)
setwd("../12_StrategiesReconstruction/")
rm(list=ls())#empty the workspace

#Load relevant packages
library(ape)
library(geiger)
library(phytools)

# Read in tree and dataset
tree <- read.nexus("TreePicidae.txt") # the tree
data <- read.csv("Full_Species_data.csv", header=T, row.names = 1)

# reorder
data <- data[tree$tip.label,]

# #Check that names match
name.check(tree, data)

# Extract variable os interest to reconstruct: here drumming strategy is this variables (discrete variables with 6 states)
Strategies <- data$AcousticClust

# Test models with different transistion probability matrices
fitER<-ace(Strategies, tree,model='ER', type='discrete')
fitSYM<-ace(Strategies, tree,model='SYM', type='discrete') # Bsed on Log-likelihood, this model (SYM) is better than ER

############## plot tree
# with tip pie charts inside
#pdf("Reconstructed_Strategies_dotsInside.pdf")
par(lwd=0.5)
plot(tree, type='fan', cex=0.5, label.offset=1.5, x.lim = c(-30,30), y.lim = c(-40,40), edge.width = 0.5)
cols<-setNames(palette()[1:length(unique(Strategies))],sort(unique(Strategies)))
nodelabels(node=1:tree$Nnode+Ntip(tree),pie=fitSYM$lik.anc,piecol=cols,cex=0.3) #pie charts at nodes
tiplabels(pie=to.matrix(Strategies,sort(unique(Strategies))),piecol=cols,cex=0.25, offset= 0.75) # Pie charts at tips
#dev.off()

# or with tip pie charts outside
# plot tree
#pdf("Reconstructed_Strategies_dotsOutside.pdf")
par(lwd=0.5)
plot(tree, type='fan', cex=0.5, label.offset=0.5, x.lim = c(-30,30), y.lim = c(-40,40), edge.width = 0.5)
cols<-setNames(palette()[1:length(unique(Strategies))],sort(unique(Strategies)))
nodelabels(node=1:tree$Nnode+Ntip(tree),pie=fitSYM$lik.anc,piecol=cols,cex=0.3) #pie charts at nodes
tiplabels(pie=to.matrix(Strategies,sort(unique(Strategies))),piecol=cols,cex=0.25, offset= 19) # Pie charts at tips
#dev.off()


fitSYM$lik.anc # outputs scaled likelihoods of each ancetstral states
nodelabels() # CAREFUL: nodelabels start at 93 --> node labelled '93' is the first line in 'fitSYM$lik.anc'
branching.times(tree) #checkup line.


# Another step is to get a quantification of the phylogenetic signal reconstructed with from the acoustic strategies
# Using geiger's fitDiscrete function for this (makes the same reconstruction as method above, but it gives Lambda in addition)
names(Strategies) <- rownames(data)
fitDiscrete(tree,Strategies,model = 'SYM', transform="lambda")


# Creating Sup Fig 4a with tip pie charts outside ## note: issues creating a pdf from console --> run what's below and export as pdf from plot window
# need to also run "Mapping_Accoustic_Phylo.R" for isolated variable (nb #38) before, to get the obj$tree object to plot

#pdf("Reconstructed_Strategies_dotsOutside.pdf")
par(lwd=0.5)
plot(tree, type='fan', cex=0.5, label.offset=0.5, x.lim = c(-30,30), y.lim = c(-40,40), edge.width = 0.5, edge.color = "transparent", tip.color = "transparent", rotate.tree = 360/92) #only way to integrate tip 'round dots' on tree. Careful, this alignement is not good --> need to rotate by one species, hence 360/92
cols<-setNames(palette()[1:length(unique(Strategies))],sort(unique(Strategies)))
#nodelabels(node=1:tree$Nnode+Ntip(tree),pie=fitSYM$lik.anc,piecol=cols,cex=0.3) #pie charts at nodes
tiplabels(pie=to.matrix(Strategies,sort(unique(Strategies))),piecol=cols,cex=0.25, offset= 16) # Pie charts at tips
plot(obj$tree,type="fan",fsize=0.55,ftype="i",add=TRUE,
     xlim=obj2$x.lim,ylim=obj2$y.lim,lwd=3,colors=obj$cols)
#dev.off()
