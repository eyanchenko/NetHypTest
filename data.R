##############################################

Game of Thrones

library(networkdata)
n=vcount(got[[4]])

A = matrix(0, ncol=n, nrow=n)

for(i in 1:n){
  
  ni = length(got[[4]][[i]][[1]])
  
  if(ni >0 ){
    for(j in 1:ni){
      
      index = as.numeric(got[[4]][[i]][[1]][j])
      
      A[i,index] <- 1
      
    }
  }
}
  

##############################################

Hospital

data(package="igraphdata")
data(rfid)
n=vcount(rfid)

A = matrix(0, ncol=n, nrow=n)

for(i in 1:n){
  
  ni = length(rfid[[i]][[1]])
  
  if(ni >0 ){
    for(j in 1:ni){
      
      index = as.numeric(rfid[[i]][[1]][j])
      
      A[i,index] <- 1
      A[index,i] <- 1
    }
  }
  

  
}

##############################################

Jazz

drat::addRepo("schochastics")
library(networkdata)
data(jazz)

n=vcount(jazz)

A = matrix(NA, ncol=n, nrow=n)

for(i in 1:n){
  A[,i] <- jazz[i]
}

##############################################

Karate

library(igraphdata)
data(karate)
n=vcount(karate)

A = matrix(0, ncol=n, nrow=n)

for(i in 1:n){

  ni = length(karate[[i]][[1]])
  
  for(j in 1:ni){
    
    index = as.numeric(karate[[i]][[1]][j])
    
    A[i,index] <- 1
  }
  
}

##############################################

Law Friends

library(networkdata)
n=vcount(law_friends)

A = matrix(0, ncol=n, nrow=n)

for(i in 1:n){
  
  ni = length(law_friends[[i]][[1]])
  
  if(ni >0 ){
    for(j in 1:ni){
      
      index = as.numeric(law_friends[[i]][[1]][j])
      
      A[i,index] <- 1
      A[index, i] <- 1
      
    }
  }
  
  
  
}

##############################################

Poltical Blogs

install.packages("remotes")
remotes::install_github("barrpet/gclust")


data("polblogs")

## Just look at one community to see if Bickel rejects and we don't


n=vcount(polblogs)

A = matrix(0, ncol=n, nrow=n)

for(i in 1:n){
  
  ni = length(polblogs[[i]][[1]])
  
  if(ni > 0){
    for(j in 1:ni){
      
      index = as.numeric(polblogs[[i]][[1]][j])
      
      A[i,index] <- 1
      A[index, i] <- 1
    }
  }
  
}

diag(A) <- 0

#plot_matrix(A*4)

##Drop unconnected components

B <- A[colSums(A)>0, colSums(A)>0]
n=dim(B)[1]

G <- as.undirected(graph.adjacency(B, mode="undirected"))

#Unconnected components

c = seq(1,n)[components(G)$membership==2]

B <- B[-c, -c]
