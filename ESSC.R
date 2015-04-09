##Extraction of Statistically Significant Communities (ESSC)
##Written by James David Wilson, 9/6/2013
##Function for extracting statistically significant communities in undirected networks
##For more information, see the paper "A Testing Based Extraction Algorithm for Identifying 
##Significant Communities in Networks", Wilson, et al. 2014 in the Annals of Applied Statistics


###ESSC Function. Run this function to extract communities
ESSC = function(Adj.Matrix,alpha, Fast = TRUE){
	degrees = rowSums(Adj.Matrix)
	n = length(degrees)
	index = 1:dim(Adj.Matrix)[1]
	#initialize
	Community = list()
	p.values = 1
	#Running Greedy Search
	if(Fast == TRUE){
	#Run Loop across entire network starting at highest degree vertex on each iteration
	
	#Begin with community of largest degree
	vertex.max.degree = which.max(degrees)[1]
	B0 = c(vertex.max.degree,which(Adj.Matrix[,vertex.max.degree]>0))
	i = 0
	#Run First Extraction loop
		i = i+1
		temp = Main.Search(Adj.Matrix,alpha,B0) #run main search
		#Update
		#index = temp$index
		Community[[i]] = temp$Community
		index = setdiff(index,temp$Community)
		
		#Update B0 and the vertex of max degree
		temp.adj = Adj.Matrix[index,index]
		temp.degs = degrees[index]
		new.max = which.max(temp.degs)
		vertex.max.degree = index[new.max]
		B0 = index[which(temp.adj[,new.max]>0)]
		print(i)
		#Function to check the overlap of the communities
		match.check = function(x){
			match = rep(0,(length(x)-1))
			for(i in 1:(length(x)-1)){
				match[i] = length(intersect(x[[length(x)]],x[[i]]))/length(x[[length(x)]])
			}
			max.match = max(match)
			return(max.match)
		}
	j = 0
	while(length(Community[[i]]) > 0 & length(index) > 0 & j < 30){
        j = j + 1
				i = i+1
				temp = Main.Search(Adj.Matrix,alpha,B0) #run main search
				if(length(temp$Community) == 0 | length(temp$Community) == n ){
					break}
					max.match = match.check(list(Community,temp$Community))
				if(max.match > 0.75){break}
				#Update
				index = setdiff(index,temp$Community)
				Community[[i]] = temp$Community
        
				if(length(index) == 1){break}
				#Update B0 and the vertex of max degree
				temp.adj = Adj.Matrix[index,index]
				temp.degs = degrees[index]
				new.max = which.max(temp.degs)
				vertex.max.degree = index[new.max]
				#B0 = index[which(temp.adj[,new.max]>0)]
				B0 = which(Adj.Matrix[,vertex.max.degree]>0)
			}
	#Keep only unique communities
	Community = unique(Community)
	}
	
	
	#Run Search across all neighborhoods
	#Here, we use a parallel for loop to run quickly (automatically uses available cores)
	if(Fast == FALSE){
		Results = foreach(i = 1:n, .export = c("Main.Search")) %dopar% 
		    Main.Search(Adj.Matrix,alpha,which(Adj.Matrix[,i]>0))
		
		#Keep only unique communities
		Results = unique(Results)
		k = length(Results)
		for(i in 1:k){
			Community[[i]] = Results[[i]]$Community
			Final.p.values[[i]] = Results[[i]]$Final.p.values
			Mean.p.value[i] = Results[[i]]$Mean.p.value
		}
		all.vertices.sign = unlist(Community)
		all.vertices.sign = all.vertices.sign[unique(all.vertices.sign)]
		}

		
v = length(Community)
Background = setdiff(1:n, matrix(unlist(Community[1:v])))
return(list(Communities = Community[1:v], Background = Background))
}


###Main Search Function

Main.Search = function(Adj.Matrix,alpha,B0){
  j = 0
  
  degrees = rowSums(Adj.Matrix)
  n = length(degrees)
  B1 = B0 - 1 #Arbitrary, ensuring B1 and B0 are not equal
  #Function for checking equality of matrices 
  vectorequal = function(x,y){
    is.numeric(x) && is.numeric(y) && length(x) == length(y) && all(x==y)
  }
  while(vectorequal(B0,B1) == FALSE & j < 200){
    j = j+1
    if(j > 1){B0 = B1}
    if(length(B0)<2){num.edges.incident = Adj.Matrix[B0,]} 
    if(length(B0) >=2){
      num.edges.incident = colSums(Adj.Matrix[B0,])
    }
    
    prob.B = sum(degrees[B0])/sum(degrees) #probability of connection to B
    
    #Handle vertices that are not connected to the remainder of the graph
    #Modify p-values of these to be 1 since they are not connected
    indx = which(degrees == 0)
    p.values = pbinom(num.edges.incident,degrees,prob.B,lower.tail = FALSE)
    p.values[indx] = 1 
    
    #Adjust the p-values according to Benjamini-Hochberg procedure
    adj.p.values = p.adjust(p.values,method="BH")
    B1 = which(adj.p.values < alpha)
    Mean.p.value = mean(adj.p.values[B1])
    if(length(B1) == 0){break}
    if(is.na(vectorequal(B0,B1))==TRUE){
      Mean.p.value = mean(p.values[B1])
      break
    }
  }
  Community = B1
  return(list(Community = Community))
}