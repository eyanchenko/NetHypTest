generateDCBM <- function(theta, p11, p12, p22, prop=0.50){
  #theta: different degrees
  #pin: prob. of connection in comm.
  #pout: prob. of connectionb between comm.
  
  n <- length(theta)
  m <- as.integer(round(n*prop, digits=0))
  
  Omega11 <- matrix(p11,  ncol = m,   nrow = m)
  Omega12 <- matrix(p12, ncol = n-m, nrow = m)
  Omega22 <- matrix(p22,  ncol = n-m, nrow = n-m) 
  
  Omega <- rbind(cbind(Omega11, Omega12), cbind(t(Omega12), Omega22) )
  
  P = as.vector(theta %*% t(theta) * Omega)
  
  A <- matrix(rbinom(n^2, 1, P), ncol=n, nrow=n)
  
  A[lower.tri(A, diag=TRUE)] <- 0
  A = A + base::t(A)
  
  return(A) 
  
}


generateDCBMK <- function(theta, pin, pout, props){
  
  n <- length(theta)
  K <- length(props)
  
  Omega <- matrix(0, nrow=n, ncol=n)
  
  idx <- 1
  # Block probabilities
  for(k in 1:K){
    intra <- idx:(idx + n*props[k]-1)
    inter <- (1:n)[-intra]
    Omega[intra, intra] <- pin  # intra-community probs.
    Omega[intra, inter] <- pout # inter-community probs.
    Omega[inter, intra] <- pout
    
    idx =  idx + n*props[k]
  }
  
  P = as.vector(theta %*% t(theta) * Omega)
  
  A <- matrix(rbinom(n^2, 1, P), ncol=n, nrow=n)
  
  A[lower.tri(A, diag=TRUE)] <- 0
  A = A + base::t(A)
  
  return(A)
}
