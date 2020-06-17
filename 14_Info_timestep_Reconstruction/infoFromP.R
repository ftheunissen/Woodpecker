infoFromP <- function(p, n) {
  # Calculates the information for n groups with p as the probability of correct classification
  
  I <- log2(n) - (1-p)*log2(n-1) + p*log2(p) + (1-p)*log2(1-p)
  
  return(I)
}