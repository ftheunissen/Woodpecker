treeCutVals2 <- function(tree, tvals, simMod, leafGroup) {
  # tree is a phylogenetic tree - an object of class phylo
  # tvals are the time points corresponding to the temporal cut in the tree
  # simMod is the ouput of the describe.simmap() summarizing the simulations of ancestral states performed by make.simmap
  # treeCutVals is a list of the same length as tvals that has number of branches and
  # the interporlated value for the discrete states modeled in ace. 
  library(utils)
  # Generate a distance for each beg and end of segment.
  # These are distances from root (0) to tip.
  segLengths <- array(0, dim=c(length(tree$edge.length),2))
  branchLengths <- array(0, dim=tree$Nnode)
  
  # Loop 1
  # Loop through all the edges until you find leaf and fill in segLengths
  rootNode <-  Ntip(tree)+1
  branchLengths[1] <-  0
  tempLength = 0
  for (i in 1:length(tree$edge.length) ) {
    segLengths[i,1] <- tempLength
    tempLength <-  tempLength + tree$edge.length[i]
    segLengths[i,2] <- tempLength
    startNode <-  tree$edge[i,1]
    endNode <- tree$edge[i,2]
    if (endNode < rootNode) {
      # Found a leaf
      if (i < length(tree$edge.length) ) {
        # Peak at the next edge to set a starting point
        nextNode <- tree$edge[i+1,1]
        tempLength <- branchLengths[nextNode-rootNode+1] 
      }
    } 
    else {
        branchLengths[endNode-rootNode+1] <- tempLength
    }
  }
  
  ## Loop2: Loop through tvals to generate probabilities for all states
  nstates <- dim(simMod$ace)[2]   # Number of states

  tXvals = list()

  for (i in 1:length(tvals)) {
    nbranches <- 0
    xvalinterp <- list()
    for (j in 1:length(tree$edge.length)) {
      if (tvals[i] > segLengths[j,1] && tvals[i] <= segLengths[j,2]) {
        nbranches <- nbranches + 1
        
        # startNode must be an node
        startNode <-  tree$edge[j,1]-rootNode+1
        xStart <- simMod$ace[startNode,]
        
        endNode <- tree$edge[j,2]
        if (endNode < rootNode) {
          # Found a leaf
          xEnd <- array(0, dim=c(1,nstates))
          colnames(xEnd) <- colnames(simMod$ace)
          xEnd[1,leafGroup[endNode]] = 1.0
        } else {
          endNode <- endNode - rootNode + 1
          xEnd <- simMod$ac[endNode,]
        }
        
        # Linear interpolation between xStart and xEnd
        dt <- tvals[i] - segLengths[j,1]
        xvalinterp[[nbranches]] <- xStart + (dt/tree$edge.length[j])*(xEnd-xStart)
        
        print(sprintf('tval= %f', tvals[i]))
        print(sprintf('\tBranch %d', nbranches))
        print(xStart)
        print(xvalinterp[[nbranches]])
        print(xEnd)
        print(' ')
      }
    }
    tXvals[[i]] <- c(nbranches, xvalinterp)
  }
  return(tXvals)
}
  
  
  