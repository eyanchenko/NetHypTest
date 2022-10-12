asymp.cut <- function(A, K=2, g0=0){
  n = dim(A)[1]
  p.hat <- sum(A)/(n*(n-1))
  eps = 0.0001
  
  cutoff = (g0 + (log(K)/n)^(1/2)/p.hat/K) * (1+eps)
  
  return(cutoff)
  
}

emp.pval2 <- function(ts, A, nsim, null=c("ER", "CL"), rn=1){
  #ts to compute p-value for, from true fit
  #A: observed adjacency matrix
  #nsim: # iterations
  #Null: What null hypothesis do you want to compare against?
  #ER: Erdos-Renyi, CL: Chung-Lu
  
  if(ts <= -1){
    warning("The observed test statistic corresponds to no communities.")
    return(1)
  }
  
  n = dim(A)[1]
  p.hat <- sum(A)/(n*(n-1))
  
  fit <- eigs(A, 1, which="LM")
  theta.hat <- abs(as.numeric(fit$vectors * sqrt(fit$values)))
  theta.hat[theta.hat > 1] <- 1
  # theta.hat <- colSums(A)/sqrt(sum(A))
  # theta.hat[theta.hat > 1] <- 1
  
  
  apply.fun <- function(i, nl){
    if(nl=="ER"){#Generate ER graph
      A.hat = generateDCBM(rep(1,n), p.hat, p.hat, p.hat)
    }else if(nl=="CL"){#Generate CL graph from degree nodes from A
      theta.boot <- sample(theta.hat, n, replace=T)
      A.hat <- generateDCBM(theta.boot, 1, 1, 1)
    }
    
    # Only keep largest connected component
    A.graph <- as.undirected(graph.adjacency(A.hat, mode="undirected"))
    comps <- components(A.graph)
    big_id <- which.max(comps$csize)
    v_ids <- V(A.graph)[comps$membership == big_id]
    A.graph <- induced_subgraph(A.graph, v_ids)
    A.hat <- as.matrix(as_adjacency_matrix(A.graph))   
    K <- length(cluster_fast_greedy(A.graph))
    
    greedy(A.hat, K, runs=rn)$e2d2
  }
  
  emp.samp <- unlist(mclapply(1:nsim, apply.fun, nl=null, mc.cores = 6))
  
  
  pval <- mean(emp.samp > ts)
  
  return(pval)
}

greedy <- function(A, K, runs=2, max.iter=100){
  
  n = dim(A)[1]
  
  # Find max e2d2 and corresponding labels
  apply.fun <- function(i){
    # Initialize communities
    comm = sample(1:K,n,replace=T)
    comm_old = comm
    
    # Community sizes
    nc = numeric(K)
    for(i in 1:K){nc[i] <- sum(comm==i)}
    
    
    # Initial values
    M   = sum(A)/2
    Xin = 0
    for(i in 1:K){Xin = Xin + sum(A[comm==i, comm==i])/2}
    Xout = M - Xin
    
    m   = 0.5*n*(n-1)
    m_in = 0
    for(i in 1:K){m_in = m_in + 0.5*nc[i]*(nc[i]-1)}
    m_out = m - m_in
    
    e2d2 = Xin/m_in - Xout/m_out
    
    run = T
    iter = 1
    while(run && iter < max.iter){
      run = F
      # Randomize node order
      node.order = sample(1:n)
      for(i in node.order){
        
        #Find community neighbors of node i
        neighs = sample(unique(comm[as.logical(A[i,])]))
        e2d2_hold = numeric(length(neighs))
        
        # Compute change in e2d2 for all neighboring communities
        Xlost = sum(A[i, comm==comm[i]])
        m_in_lost = nc[comm[i]] - 1
        
        apply.fun2 <- function(j){
          # Swap community of node i
          Xin_i = Xin - Xlost + sum(A[i, comm==neighs[j]])
          Xout_i = M - Xin_i
          
          m_in_i = m_in - m_in_lost + nc[neighs[j]] - 1*(neighs[j]==comm[i])
          m_out_i = m - m_in_i
          
          return(Xin_i/m_in_i - Xout_i/m_out_i)
        }
        
        e2d2_hold <- unlist(lapply(1:length(neighs), apply.fun2))
        
        # Swap node to whichever community yields largest e2d2 and update values
        ck = neighs[which.max(e2d2_hold)]
        
        if(length(ck)==0){
          ck = comm[i]
        }
        
        Xin = Xin - Xlost + sum(A[i, comm==ck])
        Xout = M - Xin
        
        m_in = m_in - m_in_lost + nc[ck] - 1*(ck==comm[i])
        m_out = m - m_in
        
        e2d2 <-  Xin/m_in - Xout/m_out
        nc[comm[i]] = nc[comm[i]] - 1
        nc[ck]      = nc[ck] + 1
        comm[i] <- ck
        
        
        
      }
      
      
      
      # Stop if no labels updated
      if(sum(comm==comm_old) < n){
        run = T
        comm_old = comm
      }
      
      iter = iter + 1
      
    }
    
    return(c(comm, e2d2/(M/(n*(n-1)/2))/K ))
  }
  
  out <- lapply(1:runs, apply.fun)
  
  # Keep run with max e2d2
  e2d2_vec <- numeric(runs)
  for(i in 1:runs){
    e2d2_vec[i] <- out[[i]][n+1]
  }
  
  keeper = which.max(e2d2_vec)
  
  return(list(e2d2=out[[keeper]][n+1], comm=out[[keeper]][1:n]))
  
}

spectral.pval <- function(A){
  
  require(RMTstat)

  n=dim(A)[1]
  p.hat <- sum(A)/(n*(n-1))
  
  P.hat <- p.hat - p.hat*diag(1,n)
  A.prime <- (A-P.hat)/sqrt((n-1)*p.hat*(1-p.hat))
  
  princ.eigen <- eigs_sym(A.prime,1,which="LA")[[1]]
  
  obs.stat <- n^(2/3)*(princ.eigen-2)
  return(ptw(obs.stat, beta=1, lower.tail = FALSE))
}

spectral.adj.pval <- function(A){

  require(RMTstat)
  
  n=dim(A)[1]
  p.hat <- sum(A)/(n*(n-1))
  
  P.hat <- p.hat - p.hat*diag(1,n)
  A.prime <- (A-P.hat)/sqrt((n-1)*p.hat*(1-p.hat))
  
  princ.eigen <- eigs_sym(A.prime,1,which="LA")[[1]]
  
  obs.stat <- n^(2/3)*(princ.eigen-2)
  
  mu.tw <- -1.2065335745820 #from wikipedia
  sigma.tw <- sqrt(1.607781034581) #from wikipedia
  
  emp.stats <- numeric(50)
  
  for(i in 1:50){
    A.i <- generateER(n, p.hat)
    A.i.prime <- (A.i-P.hat)/sqrt((n-1)*p.hat*(1-p.hat))
    princ.eigen.i <- eigs_sym(A.i.prime,1,which="LA")[[1]]
    emp.stats[i] <- n^(2/3)*(princ.eigen.i-2)
  }
  
  mu.theta <- mean(emp.stats)
  sigma.theta <- sqrt(var(emp.stats))
  
  theta.prime <- mu.tw + (obs.stat-mu.theta)/sigma.theta * sigma.tw
  
  return(ptw(theta.prime, beta=1, lower.tail = FALSE))
}
