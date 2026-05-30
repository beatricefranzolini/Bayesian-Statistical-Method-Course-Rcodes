library(igraph)
library(latex2exp)
#Simulate from an Erdos Renyi random graph likelihood 
set.seed(5)
n = 100 
p = 0.05
A = matrix(0, nrow = n, ncol = n)
for(u in 1:(n-1)){
  for(v in (u+1):n){
    A[u,v] = A[v,u] = rbinom(1, 1, p) 
  }
}
# create the network object
network_ER <- graph_from_adjacency_matrix(A, mode = "undirected")

# plot it
plot(network_ER, vertex.size=4, vertex.label=NA, vertex.color='blue')
#graph density 
sum(A)/(n*(n-1)) 
#node degree distribution
hist(degree(network_ER),xlim = c(0, 20), 
     main= TeX("Observed node degree distribution $V = 100, \\psi = .05$"),
     xlab = "node degree", ylab = "num. of nodes",
     col = 'blue')

plot(1:99, dbinom(1:99, size=100, prob=p),type='h')

#Simulate from a preferential attachment model likelihood 
set.seed(2)
n = 100 
m = 3 
m0 = 3
#define f up to a notmalizing constant
f <-function(A){
  network =  graph_from_adjacency_matrix(A, mode = "undirected")
  log(2*degree(network))
}
A = matrix(0, nrow = n, ncol = n)
#set the starting network with three nodes and links between 1 and 2, 2 and 3 
A[1,2] = A[2,1] = 1
A[3,2] = A[2,3] = 1
for(u in 4:n){
  prob = f(A[1:(u-1),1:(u-1)])
  node_to_connect_with_u = sample(c(1:(u-1)), m0, prob = prob)
  A[u,node_to_connect_with_u] = A[node_to_connect_with_u,u] = 1
}
# create the network object
network_BA <- graph_from_adjacency_matrix(A, mode = "undirected")

# plot it
plot(network_BA, vertex.size=4, vertex.label=NA)
plot(network_BA, vertex.size=4)
plot(network_BA, vertex.size=4, vertex.label=NA, vertex.color='blue')
#graph density 
sum(A)/(n*(n-1)) 
#node degree distribution
hist(degree(network_BA))
hist(degree(network_BA),xlim = c(0, 20), 
     main= TeX("Observed node degree distribution $V = 100$"),
     xlab = "node degree", ylab = "num. of nodes",
     col = 'blue')
mean(degree(network_BA))

#Simulate from a SBM model likelihood 
set.seed(1)
n = 100 
H = 3 
Psi = matrix(c(.09, .03, .02, .03, .09, .03, .02, .03, .09 ), nrow = 3)
A = matrix(0, nrow = n, ncol = n)
z = sample(c(1,2,3), size = 100, replace = TRUE, prob =c(0.5, 0.25, 0.25))
for(u in 1:(n-1)){
  for(v in (u+1):n){
    A[u,v] = A[v,u] = rbinom(1, 1, Psi[z[u],z[v]]) 
  }
}
# create the network object
network_SBM <- graph_from_adjacency_matrix(A, mode = "undirected")

# plot it
plot(network_SBM, vertex.size=4, vertex.label=NA, vertex.color='blue')
#graph density 
sum(A)/(n*(n-1)) 
#node degree distribution
hist(degree(network_SBM))
hist(degree(network_SBM),xlim = c(0, 20), 
     main= TeX("Observed node degree distribution $V = 100$"),
     xlab = "node degree", ylab = "num. of nodes",
     col = 'blue')
mean(degree(network_SBM))

#Simulate from a SBM model likelihood 
set.seed(2)
n = 100 
H = 3 
Psi = matrix(c(.05, .02, .1, .02, .02, .08, .1, .08, .5 ), nrow = 3)
A = matrix(0, nrow = n, ncol = n)
z = sample(c(1,2,3), size = 100, replace = TRUE, prob =c(0.45, 0.4, 0.15))
for(u in 1:(n-1)){
  for(v in (u+1):n){
    A[u,v] = A[v,u] = rbinom(1, 1, Psi[z[u],z[v]]) 
  }
}
# create the network object
network_SBM <- graph_from_adjacency_matrix(A, mode = "undirected")

# plot it
plot(network_SBM, vertex.size=4, vertex.label=NA, vertex.color='blue')
#graph density 
sum(A)/(n*(n-1)) 
#node degree distribution
hist(degree(network_SBM))
hist(degree(network_SBM),xlim = c(0, 22), 
     main= TeX("Observed node degree distribution $V = 100$"),
     xlab = "node degree", ylab = "num. of nodes",
     col = 'blue')
mean(degree(network_SBM))
