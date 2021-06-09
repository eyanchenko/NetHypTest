asym_cut <- function(n, K, p, alpha=0.05){
  
  n.seq = seq(1,n,1)
  
  logNnk <- sum(log(n.seq)) - log(factorial(K)/2) - K * sum(log(seq(1,n/K,1)))
  
  return( sqrt(8*(logNnk - log(alpha))/(n*(n-2))) / p )
  
}

emp.pval2 <- function(ts, A, nsim=200, null=c("ER", "CL", "LSM"), com_detect_alg = cluster_walktrap, ...){
  
  n = dim(A)[1]
  p.hat <- sum(A)/(n*(n-1))
  
  #Only fit MLE once and then reuse Z and beta to generate new graphs
  if(null=="LSM"){
    #Estimate beta intercept
    require(mvtnorm)
    require(latentnet)
    A.net <- as.network(A, directed=FALSE)
    A.fit <- ergmm(A.net ~ euclidean(d=1), tofit="mle")
    beta <- A.fit$mle$beta
    sigma2 <- var(A.fit$mle$Z)
  }
  
  pval = 0
  
  for(sim in 1:nsim){
    if(null=="ER"){#Generate ER graph
      A.hat = matrix(rbinom(n^2,1,p.hat), ncol=n)
      A.hat[lower.tri(A.hat, diag = TRUE)] <- 0
      A.hat = A.hat + base::t(A.hat)
    }else if(null=="CL"){#Generate CL graph from degree nodes from A
      A.hat <- generateCLa(A)
    }else if(null=="LSM"){
      A.hat <- generateLSM(n, d=1, beta, sigma2)
    }
    
    A.graph <- as.undirected(graph.adjacency(A.hat, mode="undirected"))
    
    comm.det <- com_detect_alg(A.graph, ...)
    
    temp.stat <- test.stat(A.hat, comm.det)
    
    pval = pval + as.numeric(((temp.stat-ts)>=0))/nsim
    
  }
  
  return(pval)
}

spectral.pval <- function(A){
  
  require(RMTstat)

  n=dim(A)[1]
  p.hat <- sum(A)/(n*(n-1))
  
  P.hat <- p.hat - p.hat*diag(1,n)
  A.prime <- (A-P.hat)/sqrt((n-1)*p.hat*(1-p.hat))
  
  princ.eigen <- eigs(A.prime,1)[[1]]
  
  obs.stat <- n^(2/3)*(princ.eigen-2)
  return(ptw(obs.stat, beta=1, lower.tail = FALSE))
}

spectral.adj.pval <- function(A){

  require(RMTstat)
  
  n=dim(A)[1]
  p.hat <- sum(A)/(n*(n-1))
  
  P.hat <- p.hat - p.hat*diag(1,n)
  A.prime <- (A-P.hat)/sqrt((n-1)*p.hat*(1-p.hat))
  
  princ.eigen <- eigs(A.prime,1)[[1]]
  
  obs.stat <- n^(2/3)*(princ.eigen-2)
  
  mu.tw <- -1.2065335745820 #from wikipedia
  sigma.tw <- sqrt(1.607781034581) #from wikipedia
  
  emp.stats <- numeric(50)
  
  for(i in 1:50){
    A.i <- generateER(n, p.hat)
    A.i.prime <- (A.i-P.hat)/sqrt((n-1)*p.hat*(1-p.hat))
    princ.eigen.i <- eigs(A.i.prime,1)[[1]]
    emp.stats[i] <- n^(2/3)*(princ.eigen.i-2)
  }
  
  mu.theta <- mean(emp.stats)
  sigma.theta <- sqrt(var(emp.stats))
  
  theta.prime <- mu.tw + (obs.stat-mu.theta)/sigma.theta * sigma.tw
  
  return(ptw(theta.prime, beta=1, lower.tail = FALSE))
}

test.stat2 <- function(adj, comm.det){
  n = dim(adj)[1]
  
  K = length(comm.det)
  
  if (K == 1 || K == n) {
    
    return(-1)
    
  }
  
  
  P = choose(n,2)
  P.in <- 0
  for(i in 1:K){P.in = P.in + choose(length(comm.det[[i]]),2)} 
  P.out = P-P.in
  
  m.in <- 0
  for(i in 1:K){m.in = m.in + sum(adj[comm.det[[i]],comm.det[[i]]])/2}
  m.out <- sum(adj)/2-m.in
  
  p.in <- m.in/P.in 
  p.out <- m.out/P.out
  
  p.bar <- sum(adj)/(n*(n-1))
    
  return((p.in - p.out)/p.bar)
  
}
