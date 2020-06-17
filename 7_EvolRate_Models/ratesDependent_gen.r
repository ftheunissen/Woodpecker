## ------------------------------------------------------------------------------- ##
##                                                                                 ##
##  Estimating rates dependency on reconstructed traits - J. Clavel 2016           ##
##                                                                                 ##  
## ------------------------------------------------------------------------------- ##

# tree : is a "phylo" object.
# trait : is the trait for which we want to estimate the rates.
# data : is the trait we think is driving "data" rates (default is NULL, a simple BM is estimated).
# phytrans : a tree with branch length transformed (used for the ancestral states reconstructions); if nothing is provided the default tree is used and therefore a BM process is assumed for the ancestral states reconstructions.
# model : the functional relationship between rates and variable. "exponential" is sigma_0*exp(-r*variable(t)) and "linear" is sigma_0 + r*variable(t); where sigma_0 is basal rate, r the the effect strength, t is the time.
# startvalue : starting values for the parameter seach.
# optimization: method used for the parameter search; either "BB" (see spg package), or any methods in optim. Default is "Nelder-Mead".
# control : control list for the optimization method. See ?optim.
# scaling: "none", "sd" = standardize by the sd, "total" = center by ancestral state and then standardize by sd, "normal" = raw mean centering and scaling
# subdivisions (...): number of subdivisions for numerical integration

ratesModel2 <- function(tree, trait, data=NULL, model=c("exponential","linear"), startvalue=NULL,  optimization="Nelder-Mead", control=list(maxit=50000), scaling=c("none","sd","total","normal"), intra=NULL, phytrans=NULL, echo=TRUE, REML=FALSE, ...){

  # Reorder the tree & turn the loglik option to check=FALSE
  tree <- reorder(tree,"postorder")
  model = model[1]
  scaling = scaling[1]
  if(!is.null(data)) data = as.matrix(data)
  if(!is.null(intra) && intra==FALSE) intra = NULL
  args <- list(...)
  if(is.null(args[["subdivisions"]])) subdivisions = 500L else subdivisions = args$subdivisions
  
  # Number of species
  tips <- length(tree$tip.label)
  if(is.null(data)) p <- 0 else p <- ncol(data)
  
  # Identifying tips values
  tipsIndices <- which(tree$edge[, 2] <= tips)
  indE <- p+2
  
  # Parameters
  sig <- sum(pic(trait,tree)^2)/length(trait)
  if(is.null(startvalue)){
    if(is.null(intra)){
      startvalue <- c(rep(0.0001, p), log(sig))
    }else{
      startvalue <- c(rep(0.0001, p), log(sig), sqrt(sig*0.01))
    }
  }
  
  # brownian bridge expectation path (Guindon 2013 - Syst. Bio.)
  bm <- function(s,v0,vt,t) v0 + matrix(((vt - v0)/t), ncol=p)%x%s # to vectorize the calculus
  
  # Make matrix
  if(is.null(phytrans)) phytrans <- tree # assuming brownian motion
  
    if(!is.null(data)){
      ancestral <- apply(data, 2, function(x) fastAnc(phytrans, x))
      DaAnc <- rbind(data, ancestral)
  
    if(scaling=="sd"){
      DaAnc = scale(DaAnc, center = FALSE) # just standardized by sd?
    }else if(scaling=="total"){
        # with centring = using the ancestral state as the central data
        DaAnc = scale(DaAnc, center = ancestral[1,])
    }else if(scaling=="normal"){
      DaAnc = scale(DaAnc)
    }else{
      DaAnc = DaAnc
    }
  }else{
      model = "brownian"

  }

  
  # functions for rates changes
  switch(model,
         "exponential"={
             f <- function(x, param, i, maxt){
                 exp(param[p+1]) * exp( rowSums(param[1:p]*bm(x, DaAnc[tree$edge[i,1], 1:p], DaAnc[tree$edge[i,2], 1:p], maxt)))
             }
         },
         "linear"={
           f <- function(x, param, i, maxt){
             exp(param[p+1]) + rowSums(param[1:p]*bm(x, DaAnc[tree$edge[i,1], 1:p], DaAnc[tree$edge[i,2], 1:p], maxt))
           }
         },
         "brownian"={
             f <- function(x, param, i, maxt) exp(param[p+1]) * (x/(maxt/2))
         })
  
  # Tree branches length transformations
  brTransf <- function(phy, param){
    
    # loop over the tree
    for(i in 1:nrow(phy$edge)){
        bl<-phy$edge.length[i]
        int <- try(integrate(f, lower=0, upper=bl, param=param, i=i, maxt=bl, subdivisions=subdivisions, rel.tol = .Machine$double.eps^0.05), silent = TRUE)
        if(inherits(int ,'try-error')){
            warning(as.vector(int))
            integrated <- NA_real_
        } else {
            integrated <- int$value
        }
        phy$edge.length[i] <- integrated
    }
    
    return(phy)
  }
  
  # Log-likelihood function with contrast in mvMORPH; method="pic"
  loglik<-function(param){
    
    # Compute the rates
    phy <- brTransf(phy=tree, param=param)

    if(any(is.na(phy$edge.length))) return(list(loglik=1000000))
    if(any(phy$edge.length<=0)) return(list(loglik=1000000))
    if(!is.null(intra)){
      error = param[indE]*param[indE]
      phy$edge.length[tipsIndices] <- phy$edge.length[tipsIndices] + error
    }
    # compute the loglik with the fast algorithm in mvMORPH
    if(REML==FALSE){
      loglik <- mvLL(phy, trait, method="pic", param=list(estim=FALSE, sigma=1, check=FALSE))
    }else{
      loglik <- ll_reml(trait, phy)
    }
    
    if(is.na(loglik$logl)){return(list(loglik=1000000))}
    return(list(loglik=-loglik$logl, theta=loglik$theta, sigma=loglik$sigma))
  }
  
  if(optimization=="fixed"){
    
    # transformed tree
    trans_phy <- brTransf(phy=tree, param=startvalue)
    
    # relative rates
    rates_est <- trans_phy$edge.length/tree$edge.length
    
    param <- startvalue
    
    results <- list(loglik=loglik(param)$loglik, tree=trans_phy, rates=rates_est, model=model, par=param)
    return(results)
  }else{
    # Optimization; default to Nelder-Mead but we can use L-BFGS-B with bounds or to speed up the computations
    if(optimization=="BB"){
      require(BB)
      estim<-spg(par=startvalue, fn=function(par){loglik(par)$loglik}, method=3)
    }else{
      estim<-optim(par=startvalue, fn=function(par){loglik(par)$loglik}, method=optimization, hessian=TRUE, control=control)
    }
    
    # estimated parameters
    if(!is.null(intra)) intra <- estim$par[indE]*estim$par[indE]
    # theta <- loglik(estim$par)$theta # to remove
    
    # number of parameters
    nparam <- length(estim$par) + 1
   
    LL <- estim$value
    AIC <- 2*LL+(2*nparam)
    AICc <- AIC+((2*nparam*(nparam+1))/(tips-nparam-1))
    
    
    # Check ML
    if(optimization=="BB"){
      require(numDeriv)
      hmat<-hessian(x=estim$par, func=function(par){loglik(par)$loglik})
      hess<-eigen(hmat)$value
    }else{
      hess<-eigen(estim$hessian)$values
    }
    
    if(any(hess<0)){
      hess.value<-1
    }else{
      hess.value<-0
    }
    
    # transformed tree
    trans_phy <- brTransf(phy=tree, param=estim$par)
   
    # ancestral state
    theta <- mvLL(trans_phy, trait, method="pic", param=list(estim=FALSE, sigma=1, check=FALSE))$theta
    
    # relative rates
    rates_est <- trans_phy$edge.length/tree$edge.length
    
    # transform the tree while accounting for error
     if(!is.null(intra)){
      trans_phy$edge.length[tipsIndices] <- trans_phy$edge.length[tipsIndices] + intra
     }
    
    # parameters
    if(is.null(intra)){
        if(!is.null(data)){
            param = c(exp(estim$par[p+1]),estim$par[1:p])
            coeffNames = paste("coef",1:p)
        }else{
            param = c(exp(estim$par[p+1]))
            coeffNames = NULL
        }
      names(param) = c("sigma_0", coeffNames)
    }else{
         if(!is.null(data)){
            param = c(exp(estim$par[p+1]),estim$par[1:p],intra)
            coeffNames = paste("coef",1:p)
         }else{
            param = c(exp(estim$par[p+1]),intra)
            coeffNames = NULL
         }
      names(param) = c("sigma_0", coeffNames,"error")
    }
   

    results<-list(par=param,  ancestral.state=theta,
        logl=-LL, AIC=AIC, AICc=AICc, convergence=estim$convergence, hess.values=hess.value, tree=trans_phy, rates=rates_est, model=model, intra=intra, optim_par=estim$par, REML=REML, optimization=optimization)
    
    ## ---- plot the results
    if(echo==TRUE){
      if(estim$convergence==0){  
       cat("successful convergence of the optimizer","\n")
       }else if(estim$convergence==1){  
        cat("\n","maximum limit iteration has been reached, please consider increase maxit","\n")
       }else{  
        cat("\n","convergence of the optimizer has not been reached, try simpler model","\n") 
     }
    
     if(hess.value!=0){
       cat("unreliable solution has been reached, check hessian eigenvectors or try simpler model","\n")
     }else{
       cat("R think that a reliable solution has been reached","\n")
      }
    
     if(!is.null(data)) cat("-- Summary results for",model,"dependent rate model --","\n") else cat("-- Summary results for Brownian Motion --","\n")
     cat("LogLikelihood:","\t",-LL,"\n")
     cat("AIC:","\t",AIC,"\n")
     cat("AICc:","\t",AICc,"\n")
     cat(nparam,"parameters","\n")
     cat("\n")
     cat("Estimated parameters","\n")
     cat("______________________","\n")
     print(param)
     cat("\n")
     cat("Estimated root state","\n")
     cat("______________________","\n")
     print(theta)
     cat("\n")
    }
    # Return the results
    invisible(results)

  }
}


# library(mvMORPH)
# tree <- pbtree(n=150)
# 
# trait <- rTraitCont(tree)
# data <-  cbind(rTraitCont(tree),  rTraitCont(tree))
# 
# ratesModel2(tree, trait, data, model=c("linear"))
# ratesModel2(tree, trait, data, model=c("linear"), startvalue=runif(3))
# ratesModel2(tree, trait, data, model=c("linear"), startvalue=runif(3), optimization = "BB")
# ratesModel2(tree, trait, data, model=c("exponential"))
# ratesModel2(tree, trait, data, model=c("exponential"), startvalue=runif(3))
# ratesModel2(tree, trait, data, model=c("exponential"), startvalue=runif(3), optimization = "BB")
# # plot traits
# par(mfrow=c(2,2))
# contMap(tree,trait, main="reconstructed traits values")
# contMap(tree,variable, main="reconstructed driving variable values")
# contMap(fit_1$tree,trait, main="tree with branch length scaled by the rates inferred")
# 
# # rates shifts:
# plotShifts(tree, fit_1$rates , log=F, width=5, horizontal=TRUE, main="Inferred Rates")
# fit_1$tree$edge.length
# 
# # Compare the AIC
# fit_1$AIC # Higher than BM when slope tend to zero
# mvBM(tree,trait,method="pic", echo=F, diagnostic=F)$AIC
# 
# # Estimated parameters (r is the value used for the trend):
# print(fit_1)




## Plot shift function for the mcmcTrait package

plotShifts <-
  function(phylo, chain, burnin=1000, ...){
    require(fields)
    args <- list(...)
    # options
    if(is.null(args[["fun"]])) args$fun <- mean
    if(is.null(args[["show.tip.label"]])) args$show.tip.label <- TRUE
    if(is.null(args[["horizontal"]])) args$horizontal <- TRUE
    if(is.null(args[["color"]])) args$color <- c("blue", "red")
    if(is.null(args[["scale"]])) args$scale <- FALSE
    if(is.null(args[["log"]])) args$log <- FALSE
    if(is.null(args[["palette"]])) args$palette <- FALSE
    if(is.null(args[["main"]])) args$main <- NULL
    if(is.null(args[["cex"]])) args$cex <- 0.8
    if(is.null(args[["width"]])) args$width <- 1
    
    if(inherits(chain, "mcmc")){
      tot <- nrow(chain)
      if(burnin>tot) stop("Error! the burnin value is higher than the chain length")
      chain <- chain[c(burnin:tot),-1]
      meanRate <- apply(chain, 2, args$fun)
      if(args$log==TRUE) meanRate <- log(meanRate)
    }else{
      meanRate <- chain
      if(args$log==TRUE) meanRate <- log(chain)
    }
    
    # Check the order of the tree 
    if(attr(phylo,"order")!="postorder") phylo <- reorder.phylo(phylo, "postorder")
    
    # colors mapping
    if(any(args$palette==FALSE)){
      Colors = colorRampPalette(args$color)( 100 ) 
    }else{
      Colors = args$palette
    }
    # 0 index induce error I scale it between 1 and 100
    linScale <- function(x, from, to) round( (x - min(x)) / max(x - min(x)) * (to - from) + from)
    col <- linScale(meanRate, from=1, to=100)
    
    if(args$scale==TRUE){
      phylo$edge.length <- phylo$edge.length*meanRate
    }
    plot(phylo, edge.color = Colors[col], show.tip.label = args$show.tip.label, main = args$main, cex = args$cex, edge.width= args$width)
    
    image.plot(z = as.matrix(meanRate),col = Colors,
               legend.only = T, horizontal = args$horizontal)
    
  }


## REML loglik (to be improved)
ll_reml <- function(data, tree){
  picX<-pic(data,tree,scaled=FALSE,var.contrasts=TRUE)
  logl<-sum(dnorm(picX[,1],sd=sqrt(picX[,2]),log=TRUE))
  results <- list(logl=logl, sigma=1, theta="Not estimated")
  return(results)
}
