#Source for the Stochastic Block Model
#Written by James D. Wilson; University of North Carolina at Chapel Hill
require(Rlab)
require(Matrix)
stochastic.block = function(n, k = 2, p.diff, p.between, ind = c(0, floor(n/2), n), background = FALSE, bn = 0){
  ##Generate a Stochastic block model
  #Function will create a realization of the stochastic block model where each block does not overlap. 
  #Here, the probability of an edge between blocks is given by p and the probability of connection within blocks is given by r. 
  #Input: 
  #   k = number of blocks (default = 2)
  #   p.diff = additional probability of connection within blocks (p.within - p.between)
  #   p.between = probability of connection between blocks
  #   ind = index of which communities belong to which community
  #         ex: ind = c(0, 100, 200) contains two blocks: from 1 to 100 and from 101 to 200
  #         default = c(0, floor(n/2), n) gives a stochastic 2 block model with equal sized blocks
  #   Requires the use of "Rlab" and "Matrix"
  Y = matrix(rep(0,n*k),ncol = k) 
  for(i in 1:k){
    Y[(ind[i]+1):ind[i+1],i] = 1
  }
  B = diag(rep(1,k))*p.diff + p.between * rep(1,k) %*% t(rep(1,k))
  if(background == TRUE){
    M = matrix(rep(bn,(k+1)^2),ncol = k+1)
    M[1:k,1:k] = B
    B = M
    indx = which(rowSums(Y) < 1)
    L = rep(0,n)
    L[indx] = 1
    Y = cbind(Y,L)
  } 
  #ind = c(0,cumsum(size))
  #Left Sender community
  
  #Z = 1 - Y
  expected.A = Y%*%B%*%t(Y)/2
  Adj = matrix(rbern(n^2, matrix(expected.A, ncol = 1)), ncol = n)
  #Make this symmetric
  Adj = Adj + t(Adj)
  Adj[which(Adj > 0)] = 1
  diag(Adj) = 0
  Adj = Matrix(Adj, sparse = TRUE)
  return(Adj)
}

## Examples:

# Example 1: a block model with 3 blocks, 1000 nodes, inner block probability = 0.10, between block probability = 0.05

Adjacency <- stochastic.block(1000, 3, 0.10, 0.05, ind = c(0,300, 600, 1000))

image(Adjacency)

# Example 2: a block model with 3 blocks as above but with background vertices belonging to no community.
# Here, background vertices are connected at random with probability bn = 0.05
Adjacency.background <- stochastic.block(1000, 3, 0.05, 0.05, ind = c(0,300, 600, 800), background = TRUE, bn = 0.05)

image(Adjacency.background)
