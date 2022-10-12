# NetHypTest

R code for implementing network community detection hypothesis testing based on methods proposed in the paper "A generalized hypothesis test for community structure in networks" by Eric Yanchenko and Srijan Sengupta [https://arxiv.org/abs/2107.06093]. Codes include functions to generate random graph models, asymptotic testing procedure, bootstrap testing procedure, Spectral and Spectral adjusted methods from Bickel and Sarkar (2016) and real-data examples. Please cite the paper if you use these codes.

######################################################

asymp.cut

Description

Asymptotic cutoff for the method in Yanchenko and Sengupta (2021+)

Usage

asym.cut(A, K, g0=0)

Arguments

A: network

K: number of communities

g0: null value of E2D2 parameter


######################################################

emp.pval

Description

Computes the p-value for the method in Yanchenko and Sengupta (2021+)

Usage

emp.pval(ts, A, nsim=200, null=c("ER", "CL"), rn=1)

Arguments

ts: test statistic

A: adjacency matrix

nsim: number of bootstrap simulations

null: null distribution for bootstrap iterations (ER=Erdos-Renyi, CL=Chung-Lu)

rn: number of runs of greedy maximization algorithm

######################################################

greedy

Description

Finds labels and value of maximum of E2D2 parameter in Yanchenko and Sengupta (2021+)

Usage

greedy(A, K, runs=2, max.iter=100)

A: adjacency matrix

K: number of communities 

runs: number of runs 

max.iter: maximum number of iterations

######################################################

generateDCBM

Description

Generates a random graph from an Degree corrected block model with 2 blocks

Usage

generateDCBM(theta, p11, p12, p22, prop=0.50)

Arguments

theta: degree parameters (length n=number of nodes)

p11: probability of an edge between nodes in group 1

p12: probabilty of an edge between nodes in group 1 and 2

p22: probablity of an edge between nodes in group 2

prop: proportion of nodes in group 1

######################################################

generateDCBMK

Description

Generates a random graph from an Degree corrected block model with K blocks

Usage

generateDCBMK(theta, pin, pout, props)

Arguments

theta: degree parameters (length n=number of nodes)

pin: probability of an edge within a group

pout: probabilty of an edge between groups

prop: proportion of nodes in each group (length K=number of groups)
 the network

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

Data sets


Hospital Encounters (Vanhems et. al 2013)

Interactions between patients and workers in geriatric unit of French hospital

DBLP (Gao et al., 2009 and Ji et al., 2010)

Interactions between researchers at conferences for computer science










