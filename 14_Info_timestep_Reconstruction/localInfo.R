localInfo <- function(dataSubset) {

  # Calculates the information from the confusion matrix obtained from an LDA
  # dataSubset is a data frame with Species as the grouping factor and PC1 to PC6 as measures for discrimination
  species.names <-  unique(dataSubset$Species)  #'Species' is the variable by which we want to discriminate drums
  nbirds <- length(species.names)
  ndat <- dim(dataSubset)[1]
  species.names <- factor(species.names, species.names)
  dataSubset$Species <- factor(dataSubset$Species, species.names)
  
  # Use a variable number of PCs depending on the data size
  ndf <- ndat - nbirds
  npcMax <- floor((-1+sqrt(1+8*ndf))/2)
  if (npcMax > 6) {
    npcMax <- 6
  }
  eqLDA <- 'Species~PC1'
  for (ipc in 2:npcMax) {
    eqLDA <- sprintf('%s+PC%d', eqLDA, ipc)
  }
  resLDA <- lda(as.formula(eqLDA), prior=rep(1,nbirds)/nbirds, data=dataSubset, CV=T) # LDA computation. Important for information calculation: equal prior probabilities
  
  # The confusion matrix
  cprob.tab.post=matrix(0,nrow=nbirds,ncol=nbirds)
  for (i in 1:nbirds) {
    for (j in 1:nbirds) {
      # Select the rows that correspond to a specific species
      ind.good <- (dataSubset$Species == species.names[i])
      cprob.tab.post[i,j] <- mean(resLDA$posterior[ind.good,j], na.rm=TRUE)
    }
  }
  # Go from conditional probability to joint probability distribution
  normval_matrix <- cprob.tab.post/nbirds 
  localInfoFull <- array(0, dim = nbirds)
  localInfoBinary <- array(0, dim = nbirds)
  
  # Calculate the unconditional probability p(predicted=j)
  p.pred <- array(0, dim=nbirds)
  for (j in 1:nbirds) {
    p.pred[j] = sum(normval_matrix[,j])
  }
  
  for(i in 1:nbirds){ # for each row (actual species)
    li <- array(0, dim=nbirds)
    for (j in 1:nbirds) {
      if ( cprob.tab.post[i,j] != 0.0) {
        li[j] <-  log2(cprob.tab.post[i,j]/p.pred[j])
      } else {
        li[j] <- 0.0
      }
    }
    localInfoFull[i] = sum(cprob.tab.post[i,]*li)
    pcc = cprob.tab.post[i,i]
    # localInfoBinary[i] = pcc*log2(pcc/p.pred[i]) + (1-pcc)*log2((1-pcc)/(1-p.pred[i]))
    localInfoBinary[i] = 1 + pcc*log2(pcc) + (1-pcc)*log2(1-pcc)
  }
  
  info <- list(Full=localInfoFull, Binary=localInfoBinary)
  
  return(info)
}
