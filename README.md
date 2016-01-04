# ESSC Package
An R package that implements the ESSC method to identify statistically significant communities in undirected networks. The key references for ESSC include the following papers:

- Wilson, J.D., Wang, S.  Mucha, P.J., Bhamidi, S., and Nobel, A.B. (2014) **A testing based extraction algorithm for identifying significant communities in networks**. *The Annals of Applied Statistics* Vol. 8, No. 3, 1853-1891. 
- Wilson, J.D., Bhamidi, S., and Nobel, A.B. (2013) **Measuring the statistical significance of local connections in directed networks**. *Neural Information Processing Systems Workshop on Frontiers of Network Analysis: Methods, Models and Applications*.

## Synopsis

The Extraction of Statistically Significant Communities (ESSC) package is used to identify statistically significant communities in undirected networks. The ESSC procedure 
is based on p-values or the strength of connection between a single vertex and a set of vertices under a reference distribution derived from the conditional configuration network model. 
The procedure automatically selects both the number of communities in the network and their size. Moreover, ESSC can handle overlapping communities and, unlike the majority of existing methods, 
identifies “background” vertices that do not belong to a well-defined community. The method has only one parameter, which controls the stringency of the hypothesis tests. 


## Installation

To install ESSC.R, use the following commands. Be sure to include the required packages *Matrix*, *Rlab*, and *devtools* from R version 3.1.2 or higher.

``` 

#install the latest version of devtools
install.packages("devtools")

#install and load ESSC
devtools::install_github("jdwilson4/ESSC")
library(ESSC, quietly = TRUE)

#load other required packages
library(Matrix, quietly = TRUE)
library(Rlab, quietly = TRUE)
library(devtools, quietly = TRUE)

```

## Examples
This package contains two functions

- essc(): the function used to identify communities via the ESSC procedure
- stochastic.block(): a function to generate a realization of the stochastic block model with non-overlapping communities

First, we will generate a stochastic block model with 3 non-overlapping communities.

```
net <- stochastic.block(n = 1000, k = 3, P = cbind(c(0.1, 0.01, 0.01), c(0.01, 0.1, 0.01), c(0.01, 0.01, 0.1)), sizes = c(300, 300, 400))

#view the network
image(net$Adjacency)
```

Alternatively, we can generate the same stochastic block model with randomly assigned community labels using the argument *random.community.assignment = TRUE*

```
net2 <- stochastic.block(n = 1000, k = 3, P = cbind(c(0.1, 0.01, 0.01), c(0.01, 0.1, 0.01), c(0.01, 0.01, 0.1)), sizes = c(300, 300, 400), random.community.assignment = TRUE)

#view the network. notice that the community labels are assigned at random
image(net2$Adjacency)

#now let's rearrange the network according to the true community labels and view
image(net2$Adjacency[order(net2$Membership), order(net2$Membership)])
```

Now, let's identify the communities in the generated network using the ESSC method.

First, we will use the Poisson null distribution

```
results.pois <- essc(net$Adjacency, alpha = 0.10, Null = "Poisson")
print(results.pois$Communities)

```

Then, we will try using the Binomial null distribution

```
results.bin <- essc(net$Adjacency, alpha = 0.10, Null = "Binomial")
print(results.bin$Communities)
```


## Contributors

Please send any comments or questions to the developer James D. Wilson at jdwilson4@usfca.edu. 
