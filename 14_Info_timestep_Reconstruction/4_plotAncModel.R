# Plots data and analytical information models

# External Libraries Needed
library(ggplot2)

# Go to the right place, load the data and functions needed.
setwd("../14_Info_timestep_Reconstruction/")
rm(list=ls())#empty the workspace
load(file='infoRes.RData')
source('infoFromP.R')
source('infoFromPZero.R')

dataPlot <- data.frame(tVals, ancInfo, ancInfoSD, ancInfoBi, ancInfoBiSD, ancInfoNull, ancInfoSDNull, ancInfoBiNull, ancInfoBiSDNull)
ntimes <- length(tVals)

# The plot to compare with null: only one strategy. plot.n for null model
plot.n<- ggplot(dataPlot, aes(x=tVals, y=ancInfo/log2(nspecies), colour='Actual') ) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=(ancInfo-ancInfoSD)/log2(nspecies), ymax=(ancInfo+ancInfoSD)/log2(nspecies), colour='Actual'), width=.2,
                position=position_dodge(0.05)) +
  geom_line(aes(x=tVals, y=ancInfoNull/log2(nspecies), colour='One Strat')) +
  geom_point()+
  geom_errorbar(aes(ymin=(ancInfoNull-ancInfoSDNull)/log2(nspecies), ymax=(ancInfoNull+ancInfoSDNull)/log2(nspecies), colour='One Strat'), width=.2,
                position=position_dodge(0.05))+
  scale_colour_manual(values = c("black", "grey"))

plot.n <- plot.n + labs(x="Time (mYrs)", y = "Normalized Info %") + ylim(c(0.0, 1.0)) +  xlim(c(-21, -5)) +
          theme_classic()
plot.n
#ggsave("ancInfoNull.eps", device = "eps", width = 15, height = 10, units = "cm")

# p.start is the probability for correct detection at start of simulation
p.start <- uniroot(infoFromPZero, c(0.3, 0.9), nspecies[1], ancInfo[1])

# This is the model where that probability remains constant - simulation a high pressure for preserving species identity
infoModel.start <-  array(0, dim=ntimes)
for (it in 1:ntimes) {
  infoModel.start[it] <- infoFromP(p.start$root, nspecies[it])
}

# This is a model to show how we could end up at currect with with constant
# probability of being correct.
p.end <- uniroot(infoFromPZero, c(0.3, 0.9), nspecies[ntimes], ancInfo[ntimes])
infoModel.end <-  array(0, dim=ntimes)
for (it in 1:ntimes) {
  infoModel.end[it] <- infoFromP(p.end$root, nspecies[it])
}

# The random model - no pressure 
p.2 <- (p.start$root)^(1/(nspecies[1]-1))
infoModel.rand <-  array(0, dim=ntimes)
for (it in 1:ntimes) {
  infoModel.rand[it] <- infoFromP(p.2^(nspecies[it]-1), nspecies[it])
}

dataPlot2 <- data.frame(tVals, ancInfo, ancInfoSD, infoModel.start, infoModel.rand, nspecies)

# The plot to compare with two analytical models
plot.m<- ggplot(dataPlot2, aes(x=tVals, y=ancInfo/log2(nspecies), colour='Actual') ) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=(ancInfo-ancInfoSD)/log2(nspecies), ymax=(ancInfo+ancInfoSD)/log2(nspecies), colour='Actual'), width=.2,
                position=position_dodge(0.05)) +
  geom_line(aes(x=tVals, y=infoModel.start/log2(nspecies), colour='High P')) +
  geom_line(aes(x=tVals, y=infoModel.rand/log2(nspecies), colour='No P')) +
  scale_colour_manual(values = c("black", "red", "blue"))

plot.m <- plot.m + labs(x="Time (mYrs)", y = "Normalized Info %") + ylim(c(0.0, 1.0)) +  xlim(c(-15, 0.5)) +
          theme_classic()
plot.m
#ggsave("ancInfoModel.eps", device = "eps", width = 15, height = 10, units = "cm")
