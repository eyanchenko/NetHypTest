generateER <- function(n, p){
  
  A <- matrix(0, nrow = n, ncol = n)
  
  A[col(A) > row(A)] <- rbinom(n*(n-1)/2, 1, p)
  
  A = A+t(A)
  
  return(A)
  
}

generateCLa <- function(A){
  
  n = dim(A)[1]
  m = sum(A)/2
  d = rowSums(A)
  
  #NOTE THAT IF A HAS LOTS OF EDGES, THEN SOME ENTRIES >1
  P.hat  = as.vector(d %*% base::t(d) / (2*m))
  P.hat[P.hat>1] <- 0.999
  
  CL <- matrix(rbinom(n^2, 1, P.hat), ncol=n, nrow=n)
  
  CL[lower.tri(CL, diag=TRUE)] <- 0
  
  CL = CL + t(CL)
  
  return(CL)
}


generateLSM <- function(n, d, beta, sigma2=1){
  
  expit <- function(x){1/(1+exp(-x))}
  expit <- Vectorize(expit)
  
  require(mvtnorm)
  
  #Latent space position
  z <- rmvnorm(n, mean=rep(0,d), sigma=sigma2 * diag(d))
  
  P = as.vector(expit(beta - as.matrix(dist(z, method="euclidean"))))
  
  A <- matrix(rbinom(n^2, 1, P), ncol=n, nrow=n)
  
  A[lower.tri(A, diag=TRUE)] <- 0
  A = A + t(A)
  
  return(A)
  
}
