# Wrapper Scrit to calculate information by sampling
library(MASS)
library(ggplot2)

# Change this directory to folder with your files as needed
setwd("../14_Info_timestep_Reconstruction/")
rm(list=ls())#empty the workspace
source('localInfo.R')

# Read the data that has the strategies and different points in time
load('stratMapp2.RData')  

# Read the data base of acoustic features for all birds
dataAll=read.table(file="Acoustic_alldata_postDFA.csv", header=T, sep=",",dec=".")

# Global Variables - user specified
nboot <- 30         # Repeat info calculation 30 times at each time point

# Other static variables
Strategies <- data$AcousticClust
Strategies[which(Strategies==1)] = "AC"
Strategies[which(Strategies==2)] = "DK"
Strategies[which(Strategies==3)] = "SF"
Strategies[which(Strategies==4)] = "SS"
Strategies[which(Strategies==5)] = "RS"
Strategies[which(Strategies==6)] = "IS"

nstrat <- length(unique(Strategies))
specInd <- list()
for (is in 1:nstrat) {
  specInd[[is]] <- which(Strategies == rownames(Q)[is]) 
}

# Make space for Info and SE Info
ntimes <- length(tVals)
ancInfo <- array(0, dim=ntimes)
ancInfoSD <- array(0, dim=ntimes)
ancInfoBi <- array(0, dim=ntimes)
ancInfoBiSD <- array(0, dim=ntimes)
nspecies <- array(0, dim=ntimes)

# Breaks for histogram
breaks <- seq(from=0, to=1, by=0.05)
histLI <- array(0, dim=c(ntimes, length(breaks)-1))
nbadboot <- 0

for (it in 1:ntimes) {
  infoVal <- array(0, nboot)
  infoValBi <- array(0, nboot)
  ns <- xVals[[it]][[1]]   # this is the number of "species" at this time point
  nspecies[it] <- ns
  ib <- 1
  while (ib <= nboot) {
    # Generate ns strategies - ss = species strategies sp = species
    ss <- array(0, dim=ns)
    sp <- array(0, dim=ns)
    randunif <- runif(ns) 
    for (is in 1:ns) {
      cumP <- cumsum(xVals[[it]][[is+1]])
      ss[is] <- which(randunif[is] < cumP)[1]
    }
    badboot <-  FALSE
    # for each species strategy find a random species
    for (istrat in 1:nstrat) {
      nspec <- sum(ss == istrat)
      if ( nspec > length(specInd[[istrat]]) ) {
        print(sprintf('At time slice %d, Strat %s: Asking %d sp but on %d available', it, rownames(Q)[istrat], nspec, length(specInd[[istrat]]) ))
        badboot <-  TRUE
        break
      }
      sp[ss==istrat] <-  sample(specInd[[istrat]], size = nspec, replace = FALSE)
    }
    if (badboot) {
      nbadboot <- nbadboot + 1
      if ( nbadboot > 100 ) {
        break
      }
      next
    }

    # My random sample of species
    spind <-  NULL
    for (is in 1:ns) {
      spind <-  c(spind, which(dataAll$Species == data$Species[sp[is]]))
    }
    
    # Calculate Information for subset
    dataSubset <-  dataAll[spind,]
    info <- localInfo(dataSubset)
    
    # Store mean
    infoVal[ib] <- sum(info$Full, na.rm=TRUE)/length(info$Full)
    infoValBi[ib] <- sum(info$Binary, na.rm=TRUE)/length(info$Binary)
    
    # Store histogram
    histval <- hist(info$Full/log2(ns), breaks = breaks, plot=FALSE)
    histLI[it,] <- histLI[it,] + histval$counts
    ib <- ib + 1
  }
  ancInfo[it] <- mean(infoVal, na.rm=TRUE)
  ancInfoSD[it] <- sd(infoVal, na.rm=TRUE)
  ancInfoBi[it] <- mean(infoValBi, na.rm=TRUE)
  ancInfoBiSD[it] <- sd(infoValBi, na.rm=TRUE)
}

# Calculate info at current time.
infoNow <- localInfo(dataAll)
tVals[ntimes+1] <- tVals[ntimes]+tVals[1]
nspecies[ntimes+1] <- length(unique(dataAll$Species))
tVals <- tVals - tVals[ntimes+1]   # Set current time to zero
ancInfo[ntimes+1] <- sum(infoNow$Full, na.rm=TRUE)/length(infoNow$Full)
ancInfoSD[ntimes+1] <- 0
ancInfoBi[ntimes+1] <- sum(infoNow$Binary, na.rm=TRUE)/length(infoNow$Binary)
ancInfoBiSD[ntimes+1] <- 0


dataPlot <- data.frame(tVals, ancInfo, ancInfoSD, ancInfoBi, ancInfoBiSD)

# Print Info
p<- ggplot(dataPlot, aes(x=tVals, y=ancInfo) ) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=ancInfo-ancInfoSD, ymax=ancInfo+ancInfoSD), width=.2,
                position=position_dodge(0.05))
p+labs(x="Time (mYrs)", y = "Info (Bits)")+ylim(c(0.0, 2.5)) +
  theme_classic()

# Normalized Info
p2<- ggplot(dataPlot, aes(x=tVals, y=100.0*ancInfo/log2(nspecies)) ) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=100.0*(ancInfo-ancInfoSD)/log2(nspecies), ymax=100.0*(ancInfo+ancInfoSD)/log2(nspecies)),
                width=.2, position=position_dodge(0.05))
p2+labs(x="Time (mYrs)", y = "Normalized Info (%)")+ylim(c(0.0, 100.0)) +
  theme_classic()


# Dual Plot
#pdf("DualPlot_InfoReconstruction.pdf")
p4 <- ggplot(dataPlot, aes(x=tVals)) + theme(panel.grid.major = element_blank())
  # The normalized info
  p4 <- p4 + geom_line( aes(y=100.0*ancInfo/log2(nspecies), colour='Normalized') ) 
  p4 <- p4 + geom_point(aes(y=100.0*ancInfo/log2(nspecies), colour='Normalized') ) 
  p4 <- p4 + geom_ribbon(aes(ymin=100.0*(ancInfo-ancInfoSD)/log2(nspecies), ymax=100.0*(ancInfo+ancInfoSD)/log2(nspecies), colour="Normalized"),
                alpha=.3, fill = 'blue')
  
  p4 <- p4 + labs(x="Time (mYrs)", y = "Normalized Info (%)")
  
  # The total information on a secondary axis.
  p4 <- p4 + geom_line( aes(y=ancInfo*40, colour = 'Raw'))
  p4 <- p4 + geom_point(aes(y=ancInfo*40, colour = 'Raw')) 
  p4 <- p4 + geom_ribbon( aes(ymin=(ancInfo-ancInfoSD)*40, ymax=(ancInfo+ancInfoSD)*40, colour='Raw'), alpha=.3, fill='red')
  p4 <- p4 + scale_y_continuous(sec.axis = sec_axis(~./40, name =  "Info (Bits)"), limits = c(0,100)) 

  p4 <- p4 + scale_colour_manual(values = c("blue", "red"))
  p4 <- p4 + theme_minimal()
  p4 <- p4 + theme(legend.position = c(0.9, 0.1)) 
  p4
#dev.off()
#ggsave("ancInfo.eps", device = "eps", width = 15, height = 10, units = "cm")
