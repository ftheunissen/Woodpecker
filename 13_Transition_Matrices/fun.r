## Utils functions for discrete traits analysis
## Julien Clavel - 2015

# Function to extract the aic
aicc <- function(x, n, k=NULL){
  
  if(is.null(k)){
    p <- dim(x$Q)[1]
    k <- length(unique(x$Q[-seq(1,p^2,p+1)]))
  }else{
    p <- dim(k)[1]
    vec<-unique(k[-seq(1,p^2,p+1)])
    k <- length(vec[vec!=0])
  }
  if(inherits(x,"simmap")){
    AIC <- -2*x$logL+2*k
  }else{# assuming rerooting is used
    AIC <- -2*x$loglik+2*k
  }
  
  AICc <- AIC+((2*k*(k+1))/(n-k-1))
  return(data.frame(AIC=AIC, AICc=AICc, k=k))
}

# Function to extract the aic-weights
aicw <- function(x){
  aics <- data.frame(models=1:length(x), AIC=x, diff=x)
  aics <- aics[order(-aics$AIC),]
  for(i in 1:length(x)){aics$diff[i] <- aics$AIC[i]-min(aics$AIC)}
  aics$wi <- exp(-0.5*aics$diff)
  aics$aicweigths <- aics$wi/sum(aics$wi)
  aics <- aics[order(aics$models),]
  res <- aics$aicweigths
  return(res)
}

# Compute the weighted transition matrix
Model_averaging <- function(x, weight){
  p <- length(x)
  Qmat <- lapply(1:p, function(i){x[[i]]$Q * weight[[i]]})
  Qaverage <- t(Reduce('+',Qmat)) # phytools return the transpose of the Q matrix
  diag(Qaverage) <- 0
  diag(Qaverage) <- -rowSums(Qaverage)
  return(t(Qaverage)) # revert to the same order as phytools
}

## Indices for transition matrix
p=6 # number of states to reconstruct

# ER model
ER_mod <- matrix(1,p,p)
diag(ER_mod) <- 0

# SYM model
index <- 1:(p*(p-1)/2)
SYM_mod <- matrix(0,p,p)
SYM_mod[upper.tri(SYM_mod)]<-SYM_mod[lower.tri(SYM_mod)]<-index

# ARD model
index <- 1:(p*(p-1))
ARD_mod <- matrix(0,p,p)
ARD_mod[-seq(1,p^2,p+1)]<-index

## -------------- à modifier car l'ordre des colonnes et lignes des matrices doit être l'ordre alphabétique - juste pour illustrer ici le cas de l'augmentation de l'information dans le papier.
# specific models
level_structure <- c("SF", "SS", "AC", "DK", "RS", "IS")

# Quelques examples de modèles personnalisés
# Model 1: transitions sequentielles entre strategies avec taux unique
transMod1 <- matrix(0, p, p)
row.names(transMod1) <- colnames(transMod1) <- level_structure
for(i in 1:(p-1)) transMod1[i+1,i] <- 1
for(i in 1:(p-1)) transMod1[i,i+1] <- 1
diag(transMod1) <- 0

# Model 2: transitions sequentielles avec différents taux
transMod2 <- matrix(0, p, p)
row.names(transMod2) <- colnames(transMod2) <- level_structure
for(i in 1:(p-1)) transMod2[i+1,i] <- i
for(i in 1:(p-1)) transMod2[i,i+1] <- (p-1)+i
diag(transMod2) <- 0


# Model 3: transitions sequentielles avec différents taux pour chacun des sens
transMod3 <- matrix(0, p, p)
row.names(transMod3) <- colnames(transMod3) <- level_structure
for(i in 1:(p-1)) transMod3[i,i+1] <- 1
for(i in 1:(p-1)) transMod3[i+1,i] <- 2
diag(transMod3) <- 0


