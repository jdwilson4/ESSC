#' essc
#'
#' Identify statistically significant communities in undirected networks
#' @param Adj.Matrix Adjacency matrix of the network for which you'd like to find communities
#' @param alpha False discovery rate. Value between 0 and 1 that sets an upper limit on the false discovery rate of the hypothesis testing procedure in ESSC
#' @param Null The null distribution used for comparing the observed number of connections between a single vertex and a community of vertices. Can be set to either "Binomial" or "Poisson." Default is "Binomial"
#' @param Num.Samples The number of randomly selected neighborhoods used to initiate the ESSC algorithm. The maximum number is the number of vertices in Adj.Matrix. Defaults to the number of vertices in Adj.Matrix
#' @keywords community detection, community extraction
#' @return a list containing the objects 
#' \itemize{
#'    \item Communities: a list of identified communities
#'    \item Background: a numeric vector identifying vertices that did not belong to any statistically significant community
#'    \item PValues: a n x k matrix whose (i,j)th entry is the Benjamini-Hochberg adjusted p-value expressing connectivity of node i to community j
#' }
#'@references
#'\itemize{
#'     \item Wilson, James D.,Wang, Simi, Mucha, Peter J., Bhamidi, Shankar, and Nobel, Andrew B. (2014). “A testing based extraction algorithm for identifying significant communities in networks.” 
#'     The Annals of Applied Statistics Vol. 8, No. 3, 1853-1891 
#' } 
#' @author James D. Wilson
#' @examples
#' net <- stochastic.block(n = 1000, k = 3, P = cbind(c(0.1, 0.01, 0.01), c(0.01, 0.1, 0.01), c(0.01, 0.01, 0.1)), sizes = c(300, 300, 400))
#' results <- essc(net$Adjacency, alpha = 0.10, Null = "Poisson")
#' print(results$Communities)
#' @export 

###essc Function. Run this function to extract communities
essc = function(Adj.Matrix, alpha, Null = c("Binomial", "Poisson"), Num.Samples = nrow(Adj.Matrix)){
    symdiff = function(X,Y){setdiff(union(X,Y),intersect(X,Y))}
    degrees <- rowSums(Adj.Matrix)
    n <- length(degrees)
    index <- 1:dim(Adj.Matrix)[1]
    extractFrom <- sample(index, Num.Samples, replace = FALSE)
    
    #initialize
    Community <- rep(list(NULL), Num.Samples)
    pvals <- rep(list(NULL), Num.Samples)
    unq <- TRUE
    nodes <- matrix(0, 2, Num.Samples)
    nodes[1, ] <- extractFrom
    which_match <- 0

    
    i <- 1
    count <- 0
    #Running long loop
    for(node in extractFrom){
        count <- count+1
        #cat("######################################################\n")
        #cat(paste0("extraction ",count,"\n"))
        B0 <- c(node,which(Adj.Matrix[node,]>0))
        temp <- Main.Search(Adj.Matrix, alpha, B0, Null)
        
        if(length(temp$Community) > 1){
            if(i > 1){
                unq <- TRUE
                which_match <- 0
                for(j in 1: i){
                    if(length(symdiff(Community[[j]], temp$Community)) == 0){
                        unq <- FALSE
                        which_match <- j
                    }
                }
            }
            if(unq){
                Community[[i]] <- temp$Community
                pvals[[i]] <- temp$PValues
                nodes[2, count] <- i
                i = i+1
            } else{
                nodes[2, count] <- which_match
            }
        } else{
            nodes[2, count] <- 0
        }
    }
    
    Community <- Community[unlist(lapply(Community, function(comm) length(comm) > 0))]
    Background <- setdiff(1:n, unlist(Community))
    pvals <- pvals[unlist(lapply(pvals, function(comm) length(comm) > 0))]
    pvals <- matrix(unlist(pvals), nrow = n, byrow = T)
    
    #divide each of these by the sum of the p-values for each community
    tot.values <- rowSums(pvals)
    normalized_pvals <- pvals / tot.values
    
    if(length(Community) == 0){
      binary <- NULL
      Community <- list(NULL)
      pvals <- NULL
    }
    if(length(Community) > 0){
      #Creating a binary matrix indicating communities
      binary <- matrix(rep(0, n*length(Community)), nrow = n)
      for(i in 1:length(Community)){
        binary[Community[[i]], i] <- 1
      }
    }
    
    return(list(Communities = Community, Background = Background, 
                PValues = as.matrix(pvals), Indicator_Matrix = binary))
}

###Single Search Function

Main.Search = function(Adj.Matrix,alpha,B0, Null){
    j <- 0
    degrees <- rowSums(Adj.Matrix)
    n <- length(degrees)
    #Ensuring B1 and B0 are not equal to start the loop
    B1 <- integer(0)
    #Function for checking equality of matrices 
    vectorequal = function(x,y){
        is.numeric(x) && is.numeric(y) && length(x) == length(y) && all(x==y)
    }
    #cat("###### Starting Main.Search Loop ######\n")
    while(vectorequal(B0,B1) == FALSE & j <= 30){
        
        j <- j + 1
        if(j < 5||j%%5==0)
           #cat(paste("iteration",j,"\n"))
        if(j > 1)
            B0 <- B1
	  if(length(B0) > 1)
	        duBs <- rowSums(Adj.Matrix[, B0])        
	  if(length(B0)==1)
		  duBs <- Adj.Matrix[, B0]
        pB <- sum(degrees[B0])/sum(degrees) #probability of connection to B
        if(Null == "Binomial")
            pvals <- pbinom(duBs, degrees, pB, lower.tail = FALSE)
        if(Null == "Poisson")
            pvals <- ppois(duBs, degrees*pB, lower.tail = FALSE)
        pvals[degrees == 0] <- 1
        #cat(paste0("length pvals is ",length(pvals),"; mean is ",round(mean(pvals),2),"\n"))
        pvals_bh <- pvals*n/rank(pvals)
        if(sum(pvals_bh <= alpha) == 0)
        {
          B1 <- integer(0)
          #cat("B1 reached size 0 in Main.Search\n")
          break
        }
        threshold <- max(pvals[pvals_bh <= alpha])
        B1 <- which(pvals <= threshold)
        if(length(B1) < 1 ){break}
    }
    Community <- B1
    if(j > 30){
        #cat("Main.Search failed to converge before 30 iterations\n")
        Community <- integer(0)
    }
    
    return(list(Community = Community, PValues = pvals_bh))
}


