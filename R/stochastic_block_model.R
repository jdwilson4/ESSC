#' stochastic.block
#' 
#' Generate a realization of the stochastic block model with non-overlapping blocks
#' @param n Number of nodes
#' @param k Number of blocks. Default is 2
#' @param P k x k matrix where the (i,j)th entry specifies the probability of connection between community i and j
#' @param sizes an array of length k where the ith entry specifies the size of the ith block. This must sum to n
#' @param random.community.assignment logical specifying whether or not the labels are assigned randomly. If FALSE, labels are assigned in order starting with the upper left corner of the matrix. Default is FALSE
#' @param vector logical specifying whether or not the returned adjacency matrix should be returned as a vector of edge variables. Default is FALSE
#' @import Matrix
#' @import Rlab
#' @keywords block model, community detection, SBM
#' @return a list containing the objects 
#' \itemize{
#'    \item Adjacency: either an adjacency matrix of the generated realization, or if vector == TRUE, an edgelist of the generated realization
#'    \item Membership: a numeric vector specifying the community label of each vertex
#' }
#' @author James D. Wilson
#' @examples
#' net <- stochastic.block(n = 300, k = 2, P = cbind(c(0.1, 0.01), c(0.01, 0.1)), sizes = c(100, 200))
#' image(net$Adjacency)
#' @export

#require(Rlab)
#require(Matrix)
stochastic.block = function(n, k = 2, P, sizes, random.community.assignment = c(FALSE, TRUE), vector.output = c(FALSE, TRUE)){
  #Input: 
  #   n = number of nodes
  #   k = number of blocks (default = 2)
  #   P = k x k matrix where the (i,j)th entry specifies the probability of connection between community i and j
  #   sizes = an array of length k where the ith entry specifies the size of the ith block. This must sum to n
  #   random.community.assignment = logical specifying whether or not the labels are assigned randomly. If FALSE, labels
  #                       are assigned in order starting with the upper left corner of the matrix.
  #   vector = logical specifying whether or not the returned adjacency matrix should be returned as a vector
  #            of edge variables. 
  #Output: 
  #   A list containing three components:
  #     Adjacency: the adjacency matrix of the generated network
  #     Membership: an n x 1 vector specifying community membership of each node
  
  if(sum(sizes) != n){
    stop("argument sizes must sum to n")
  }
  if(length(sizes) != k){
    stop("argument sizes must be of length k")
  }
  
  #Generate Membership vector
  Membership <- rep(1, n)
  for(i in 2:k){
    z <- cumsum(sizes)
    Membership[(z[i - 1] + 1):z[i]] <- i
  }
  
  vector.output <- vector.output[1]
  
  
  Y <- matrix(rep(0, n*k), ncol = k)
  index <- list()
  possible <- 1:n
  
  random.community.assignment <- random.community.assignment[1]
  #assign vertex labels randomly if random.assignment == TRUE
  if(random.community.assignment == TRUE){
    for(i in 1:k){
      index[[i]] <- sample(possible, sizes[i])
      possible <- setdiff(possible, index[[i]])
    }
  }
  
  #assign vertex labels in order if random.assignment == FALSE
  if(random.community.assignment == FALSE){
    for(i in 1:k){
      index[[i]] <- possible[1:sizes[i]]
      possible <- setdiff(possible, index[[i]])
    }
  }
  
  for(i in 1:k){
    Y[index[[i]], i] = 1
  }
  
  expected.A <- Y%*%P%*%t(Y)
  
  #In case any of these probabilities are greater than 1
  expected.A[which(expected.A > 1)] <- 1
  
  Adj <- Matrix(rbinom(n^2, size = 1, matrix(expected.A, ncol = 1)), ncol = n, sparse = TRUE)
  diag(Adj) <- 0
  #Adj <- Matrix(Adj, sparse = TRUE)
  
  #Return results
  if(vector.output == TRUE){
    return(list(Adjacency = as.vector(Adj), Membership = Membership))
  }
  if(vector.output == FALSE){
    return(list(Adjacency = Adj, Membership = Membership))
  }
}

