treeLengths <- function(tree) {
  # tree is a phylogenetic tree - an object of class phylo
  # treelents is a 1-d array of all the lengths from the root node to the leaf
  
  # Assign space for treelengths and branchlengths
  # These are the distances to the tips (treelengths) and nodes from root
  tLengths <- array(0, dim=Ntip(tree))
  branchLengths <- array(0, dim=tree$Nnode)
  
  # Loop through all the edges until you find leaf
  rootNode <-  Ntip(tree)+1
  branchLengths[1] <-  0
  tempLength = 0
  for (i in 1:length(tree$edge.length) ) {
    tempLength <-  tempLength + tree$edge.length[i]
    startNode <-  tree$edge[i,1]
    endNode <- tree$edge[i,2]
    if (endNode < rootNode) {
      # Found a leaf
      tLengths[endNode] <- tempLength
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
  return(tLengths)
}
  
  
  