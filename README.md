# NetHypTest

R code for implementing network homophily hypothesis testing based on methods proposed in the paper "A model-agnostic hypothesis test for community structure and homophily in networks" by Eric Yanchenko and Srijan Sengupta [ArXiv link]. Codes include functions to generate random graph models, asymptotic testing procedure, bootstrap testing procedure, Spectral and Spectral adjusted methods from Bickel and Sarkar (2016) and real-data examples. Please cite the paper if you use these codes.

######################################################

asym.cut

Description

Asymptotic cutoff for the method in Yanchenko and Sengupta (2021+)

Usage

asym.cut(n, K, p, alpha=0.05)

Arguments

n: number of nodes

K: number of communities

p: overall probability of edge in graph

alpha: level of test


######################################################

emp.pval

Description

Computes the p-value for the method in Yanchenko and Sengupta (2021+)

Usage

emp.pval(ts, A, nsim=200, null=c("ER", "CL", "LSM"), com_detect_alg = cluster_walktrap, ...)

Arguments

ts: test statistic

A: adjacency matrix

nsim: number of bootstrap simulations

null: null distribution for bootstrap iterations (ER=Erdos-Renyi, CL=Chung-Lu, LSM=latent space model)

com_detect_alg: community detection algorithm to be used during bootstrap iterations

...: extra parameters for the community detection algorithm

######################################################

generateER

Description

Generates a random graph from an Erdos-Renyi model 

Usage

generateER(n,p)

Arguments

n: number of nodes in the network

p: probability of an edge between nodes

######################################################

generateCL

Description

Generates a random graph from a Chung-Lu model based on an existing adjacency matrix

Usage

generateCL(A)

Arguments

A: adjacency matrix from which to extract the degree distribution

######################################################

generateLSM

Description

Generates a random graph from a Latent Space Model

Usage

generateLSM(n, d=1, beta, sigma2=1)

Arguments

n: number of nodes in the network

d: dimension of latent space

beta: parameter that controls the overall sparsity. Larger values means more sparse

sigma2: variance parameter for generating latent positions


######################################################

spectral.pval

Description

Computes the p-value for the unadjusted method (Algorithm 1) from Bickel and Sarkar (2016)

Usage

spectral.pval(A)

Arguments

A: adjacency matrix

######################################################

spectral.adj.pval

Description

Computes the p-value for the adjusted method (Algorithm 2) from Bickel and Sarkar (2016)

Usage

spectral.adj.pval(A)

Arguments

A: adjacency matrix

######################################################

test.stat

Description

Computes the test statistic for the method in Yanchenko and Sengupta (2021+)

Usage

test.stat(adj, comm.det)

Arguments

adj: adjacency matrix

comm.det: community labels

######################################################

Data sets

Game of Thrones 

Co-appearance data for Season 4 of Game of Thrones

Hospital Encounters (Vanhems et. al 2013)

Interactions between patients and workers in geriatric unit of French hospital

Jazz (Gleiser and Danon 2003)

Common musicians between jazz bands

Karate (Zachary 1977)

Social connection between members of a karate club

Law Friends (Snijder et. al 2006)

Friendships among attorneys at a law firm

Political Blogs (Adamic and Glance 2005)

Political blogs citation network










