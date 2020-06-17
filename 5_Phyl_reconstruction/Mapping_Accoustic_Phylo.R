#Change working directory to folder with your files (change accordingly)
setwd("../5_Phyl_reconstruction/")
rm(list=ls())#empty the workspace

#Load relevant packages
library(ape)
library(geiger)
library(phytools)

#L Read in the tree (requires ape)
Picidae <-read.nexus(file = 'TreePicidae.txt')

states <-read.csv("Acoustic_species_postDFA.csv", header=T)
states$median_LI <- states$median_LI/log2(92)
states[is.na(states)] <- 0


# Ancestral state reconstruction (requires phytools)
# with extraction of phylogenetic signal (pagel's lambda values) and most ancestral state
phylSigMatrix <- matrix(0,nrow=35,ncol=6) #22 acoustic variables + 6 PCs + 6LDs + LI
colnames(phylSigMatrix) <- c("Variable", "lambda", "p val lambda signal", "Ancestral_state","low_CI95", "high_CI95")

pdf(paste("PhylogeneticMapping.pdf", sep=''))

for (i in 16:50){ #variables of interest: acoustic variables + PCs + LDs + LI
  states2 <- data.frame(states[,i])
  rownames(states2)<- states$Species_FullName
  states3 <- states2[match(Picidae$tip.label,rownames(states2)),]
  names(states3) <- Picidae$tip.label
  fit<-fastAnc(Picidae, states3, vars=TRUE,CI=TRUE)

  # K-value/lambda calculation
  lambdaVal <- phylosig(Picidae, states3, method="lambda", test=TRUE, nsim=1000)
  roundedLambda <- round(lambdaVal$lambda, digits=3)
  phylSigMatrix[i-15,1] <- names(states[i])
  phylSigMatrix[i-15,2] <- lambdaVal$lambda
  phylSigMatrix[i-15,3] <- lambdaVal$P
  
  # Ancestral state calculation
  ancestralState <- round(fit$ace[names(fit$ace)==93], digits=3) # value 93 gives the label of the ancestral node, which can be plotted with 'nodelabels()'
  
  phylSigMatrix[i-15,4] <- ancestralState
  phylSigMatrix[i-15, 5:6] <- round(fit$CI95[1,], digits=3)
  
  obj<-contMap(Picidae,states3,plot=FALSE)
  n<-length(obj$cols)

#### if want to use defined colorRamp instead of color ramp:
  obj$cols[1:n]<-colorRampPalette(c("gray100", "gray50", "yellow","goldenrod1", "orange","chocolate", "chocolate4"), space="Lab", interpolate = c("linear"))(n)

  h<-max(nodeHeights(obj$tree))
  offset.factor<-1.01 ## increase this for greater offset
  plotTree(rescale(obj$tree,model="depth",depth=offset.factor*h),
           color="transparent",ftype="i",type="fan",fsize=0.7,lwd=3)

  par(fg="transparent")
  obj2<-get("last_plot.phylo",envir=.PlotPhyloEnv)
#  plot.new()
  plot(obj$tree,type="fan",fsize=0.7,ftype="i",add=TRUE,
       xlim=obj2$x.lim,ylim=obj2$y.lim,lwd=3,colors=obj$cols)
  par(fg="black")
  add.color.bar(12, obj$cols, title='', lims=round(obj$lims,digits = 2), prompt=FALSE,
                x=0.9*par()$usr[1],y=0.95*par()$usr[3], lwd = 6, outline = F, fsize= 1)    
  
  mtext(paste("Var =",names(states[i])), at=30, line=-1, cex=1, las=1)
  mtext(paste("AncState =",ancestralState), at=-35, line=-1, cex=1, las=1)
  mtext(paste("     lambda =",roundedLambda), at=30, line=-35, cex=1, las=1)

}
dev.off()
write.csv(phylSigMatrix, "LambdaVal_AncState_CI95.csv")
#write.csv(fit$ace, paste0('AncStates_',names(states[i]),'.csv'))

# adding an extra plot for PC1 to prepare Sup Fig 4a (showing strategies around tree)
states2 <- data.frame(states[,38]) # PC1
rownames(states2)<- states$Species_FullName
states3 <- states2[match(Picidae$tip.label,rownames(states2)),]
names(states3) <- Picidae$tip.label
fit<-fastAnc(Picidae, states3, vars=TRUE,CI=TRUE)
obj<-contMap(Picidae,states3,plot=F)

#### if want to use greyscale instead of color ramp:
n<-length(obj$cols)
# obj$cols[1:n]<-grey(0:(n-1)/(n-1))
#### if want to use defined colorRamp instead of color ramp:
obj$cols[1:n]<-colorRampPalette(c("gray100", "gray50", "yellow","goldenrod1", "orange","chocolate", "chocolate4"), space="Lab", interpolate = c("linear"))(n)

#pdf("Outside_Phylbars_PC1.pdf")
add=T
plot.new()
plotTree.wBars(obj$tree,states3, method="plotSimmap", lwd=3,
               tip.labels=T,fsize=0.6,colors=obj$cols,type="fan",scale=0.002)
#dev.off()





# If need to retrieve an ancestral state value at a specific node, 
# do not use 'writeAncestors' function (issue with sorting) but rather
# plot labels with 'nodelabels()' and then check appropriate node value with fit$ace


########################################################
# For testing outside loop
states2 <- data.frame(states[,50]) # change states[XX] to the appropriate variable number; example here: col 50 is normalized LI
rownames(states2)<- states$Species_FullName
#write.csv(states2, 'testStates2.csv') # check line
states3 <- states2[match(Picidae$tip.label,rownames(states2)),]
names(states3) <- Picidae$tip.label
#write.csv(states3, 'testStates3.csv') # check line
fit<-fastAnc(Picidae, states3, vars=TRUE,CI=TRUE)

obj<-contMap(Picidae,states3,plot=F)

#### if want to use greyscale instead of color ramp:
n<-length(obj$cols)
# obj$cols[1:n]<-grey(0:(n-1)/(n-1))
#### if want to use defined colorRamp instead of color ramp:
obj$cols[1:n]<-colorRampPalette(c("gray100", "gray50", "yellow","goldenrod1", "orange","chocolate", "chocolate4"), space="Lab", interpolate = c("linear"))(n)

plot(obj, type="fan", legend=0.7*max(nodeHeights(Picidae)), fsize=c(0.7,0.9))
#nodelabels(cex = 0.5)
ancestralState <- round(fit$ace[names(fit$ace)==93], digits=3)

Lambdaval <- phylosig(Picidae, states3, method="lambda", test=TRUE, nsim=999)
roundedLamdaval <- round(Lambdaval$lambda, digits=3)

# write info into plot:
mtext(paste("Var =",names(states[50])), at=30, line=-1, cex=1, las=1) # change states[XX] to the appropriate variable number
mtext(paste("AncState =",ancestralState), at=-35, line=-1, cex=1, las=1)
mtext(paste("     lambda =",roundedLamdaval), at=30, line=-35, cex=1, las=1)
