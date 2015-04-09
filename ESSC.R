##Extraction of Statistically Significant Communities (ESSC)
##Original Code by James David Wilson; modified code by John Palowich
##Function for extracting statistically significant communities in undirected networks
##For more information, see the paper "A Testing Based Extraction Algorithm for Identifying 
##Significant Communities in Networks", Wilson, et al. 2014 in the Annals of Applied Statistics

symdiff = function(X,Y){setdiff(union(X,Y),intersect(X,Y))}

###ESSC Function. Run this function to extract communities
ESSC = function(Adj.Matrix,alpha, Null = "Binomial", Num.Samples = nrow(Adj.Matrix)){
    degrees = rowSums(Adj.Matrix)
    n = length(degrees)
    index = 1:dim(Adj.Matrix)[1]
    extractFrom = sample(index,Num.Samples,replace = FALSE)
    #initialize
    Community = rep(list(NULL),Num.Samples)
    unq = TRUE
    nodes = matrix(0,2,Num.Samples)
    nodes[1,] = extractFrom
    which_match = 0

    
    i = 1
    count = 0
    #Running long loop
    for(node in extractFrom){
        count = count+1
        cat("######################################################\n")
        cat(paste0("extraction ",count,"\n"))
        B0 = c(node,which(Adj.Matrix[node,]>0))
        temp = Main.Search(Adj.Matrix, alpha, B0, Null)
        
        if(length(temp$Community)>1){
            if(i>1){
                unq = TRUE
                which_match = 0
                for(j in 1:i){
                    if(length(symdiff(Community[[j]],temp$Community))==0){
                        unq = FALSE
                        which_match = j
                    }
                }
            }
            if(unq){
                Community[[i]] = temp$Community
                nodes[2,count] = i
                i = i+1
            } else{
                nodes[2,count] = which_match
            }
        } else{
            nodes[2,count] = 0
        }
    }
    
    Community = Community[unlist(lapply(Community,function(comm)length(comm)>0))]
    Background = setdiff(1:n,unlist(Community))
    
    if(length(Community) == 0)
        Community = list(NULL)
    
    return(list(Communities = Community,Background = Background,nodes = nodes))
}

###Main Search Function

Main.Search = function(Adj.Matrix,alpha,B0, Null){
    j = 0
    degrees = rowSums(Adj.Matrix)
    n = length(degrees)
    #Ensuring B1 and B0 are not equal to start the loop
    B1 = integer(0)
    #Function for checking equality of matrices 
    vectorequal = function(x,y){
        is.numeric(x) && is.numeric(y) && length(x) == length(y) && all(x==y)
    }
    cat("###### Starting Main.Search Loop ######\n")
    while(vectorequal(B0,B1) == FALSE & j <= 30){
        
        j = j+1
        if(j<5||j%%5==0)
            cat(paste("iteration",j,"\n"))
        if(j>1)
            B0 = B1
	  if(length(B0)>1)
	        duBs = rowSums(Adj.Matrix[,B0])        
	  if(length(B0)==1)
		  duBs = Adj.Matrix[,B0]
        pB = sum(degrees[B0])/sum(degrees) #probability of connection to B
        if(Null == "Binomial")
            pvals = pbinom(duBs,degrees,pB,lower.tail = FALSE)
        if(Null == "Poisson")
            pvals = ppois(duBs,degrees*pB,lower.tail = FALSE)
        pvals[degrees==0] = 1
        #cat(paste0("length pvals is ",length(pvals),"; mean is ",round(mean(pvals),2),"\n"))
        pvals_bh = pvals*n/rank(pvals)
        if(sum(pvals_bh<=alpha)==0){B1 = integer(0);cat("B1 reached size 0 in Main.Search\n");break}
        threshold = max(pvals[pvals_bh<=alpha])
        B1 = which(pvals<=threshold)
        if(length(B1)<1){break}
    }
    Community = B1
    if(j>30){
        cat("Main.Search failed to converge before 30 iterations\n")
        Community = integer(0)
    }
    return(list(Community = Community))
}


