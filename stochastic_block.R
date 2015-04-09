install.packages("Rlab")
install.packages("Matrix")
require(Rlab)
require(Matrix)
stochastic.block = function(n,k,p,r,ind,background = FALSE, bn = 0){
  ##Generate a Stochastic block model
  #Function will create a realization of the stochastic block model where each block does not overlap. Here, the probability of an edge between blocks is given by p and the probability of connection within blocks is given by r. (These can be replaced by vectors as well for different probabilities)
  #Written by James Wilson
  #Input: k = number of blocks
  #p = additional probability of connection within blocks (i.e amount over r)
  #r = probability of connection between blocks
  #size = k vector of sizes for each block
  #Requires the use of "Rlab" and "Matrix"
  Y = matrix(rep(0,n*k),ncol = k) 
  for(i in 1:k){
    Y[(ind[i]+1):ind[i+1],i] = 1
  }
  B = diag(rep(1,k))*p + r* rep(1,k) %*% t(rep(1,k))
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
  Adj = matrix(rbern(n^2,matrix(expected.A,ncol = 1)),ncol = n)
  #Make this symmetric
  Adj = Adj + t(Adj)
  Adj[which(Adj > 0)] = 1
  diag(Adj) = 0
  Adj = Matrix(Adj, sparse = TRUE)
  return(Adj)
}