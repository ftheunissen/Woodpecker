library(phytools)

rm(list=ls())#empty the workspace
setwd("../1_Drum_nonDrum/")
load("../1_Drum_nonDrum/MG10012018.Rdata")
View(PP) # probability distribution for each species of being drummer vs non-drummer vs occasionnal

#extract probabilities on second node from the simulation run used for data in the manuscript
obj$ace[2,]
# result:
# a     b     c 
# 0.075 0.340 0.585 



# As this is simulation-based, can do it endlessly using the following code. Will give slightly different values but altogether very similar
# to what reported in the MS, thus leaving the conclusions and interpretations unchanged (ancestral woodpecker drummed).

Picidae <- read.nexus("TreePicidae.txt")

# # model ER --> not chosen as not as well supported as the 'SYM' model: see below
# btrees<-make.simmap(Picidae,PP,nsim=200,model="ER")
# BB<-describe.simmap(btrees,plot=TRUE)
# nodelabels(cex=0.5)
# BB$ace # common ancestor of woodpeckers (exclusing Jynx) has label #211
# BB$ace[2,] # probabilities on second node (label #211)


#OR: model SYM
symtrees2<-make.simmap(Picidae,PP,nsim=200,model="SYM", Q="mcmc") #200simulations is sufficient
obj<-describe.simmap(symtrees,plot=T)
nodelabels(cex=0.5)
obj$ace # common ancestor of woodpeckers (exclusing Jynx) has label #211
obj$ace[2,] # probabilities on second node (label #211)


