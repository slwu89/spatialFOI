# -------------------------------------------------------------------------------- #
#
#   Equilibrium calculation of lambda under movement
#   January 2020
#   Sean Wu (slwu89@berkeley.edu)
#
# -------------------------------------------------------------------------------- #

rm(list=ls());gc()
set.seed(694892397L)

nodes <- c("a","a","a","b","b","b")
a_set <- which(nodes=="a")
b_set <- which(nodes=="b")

n <- length(nodes)
gamma <- 1/rgamma(n = n,shape = 20,rate = 20*(1/3))
g <- 1/rgamma(n = n,shape = 30,rate = 30*(1/10))

Fmat <- matrix(data = runif(n = n*n),nrow = n,ncol = n)
diag(Fmat) <- 0
Fmat <- Fmat/rowSums(Fmat)

M <- rpois(n = n,lambda = 50)

# solve for lambda
lambda <- rep(0,length(a_set))
for(i in a_set){
  lambda[i] <- ((g[i] + gamma[i])*M[i]) - sum(gamma[b_set]*Fmat[b_set,i]*M[b_set])
}