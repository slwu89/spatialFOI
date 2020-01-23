# -------------------------------------------------------------------------------- #
#
#   when is taking averages misleading? mosquito to human bit
#   January 2020
#   Sean Wu (slwu89@berkeley.edu)
#
# -------------------------------------------------------------------------------- #

rm(list=ls());gc()
set.seed(4369576L)

# this many habitats
A <- 10

# this many humans
H <- 1e3


# -------------------------------------------------------------------------------- #
#   Human bits
# -------------------------------------------------------------------------------- #

# b param for each
b_mean <- 0.55
b_sd <- 0.15
b_shape1 <- (((1-b_mean)/(b_sd^2)) - (1/b_mean)) * (b_mean^2)
b_shape2 <- b_shape1 * ((1/b_mean)-1)
b <- rbeta(n = H,shape1 = b_shape1,shape2 = b_shape2)

# biting heterogeneities
Z_mean <- 1
Z_sd <- 1
Z_shape <- (Z_mean^2)/(Z_sd^2)
Z_rate <- Z_mean/(Z_sd^2)
Z <- rgamma(n = H, shape = Z_shape, rate = Z_rate)


# -------------------------------------------------------------------------------- #
#   Mosquito bits
# -------------------------------------------------------------------------------- #

# just need the d'th row of the K matrix since that's enough to evaluate misleading or not
K_A2d <- runif(n=A)

# how many bites they are giving
EIR_A <- rep(10,A)


# -------------------------------------------------------------------------------- #
#   full sampling
# -------------------------------------------------------------------------------- #

sample_inf <- function(){

  # sample how many bites arise from each a \in A
  beta_A_hat <- rpois(n=A,lambda=EIR_A)

  # how many bites came here
  beta_d_hat <- sum(rbinom(n=A, size=beta_A_hat, prob=K_A2d))

  # how many bites each person got
  beta_h_Hat <- as.vector(rmultinom(n = 1,size = beta_d_hat,prob = Z))

  # infected?
  p_inf <- (1 - (1 - b)^beta_h_Hat)
  rbinom(n = H,size = 1,prob = p_inf)
}


# -------------------------------------------------------------------------------- #
#   marginal sampling
# -------------------------------------------------------------------------------- #

marginal_inf <- function(){
  lambda_h <- (Z/sum(Z)) * sum(EIR_A*K_A2d) * b
  p_inf <- 1-exp(-lambda_h)
  rbinom(n = H,size = 1,prob = p_inf)
}


# -------------------------------------------------------------------------------- #
#   simulate a large number of times in parallel
# -------------------------------------------------------------------------------- #

library(parallel)
parseed <- 19306148L
N_mc <- 1e4
cl <- makeCluster(4)
clusterSetRNGStream(cl,parseed)
clusterExport(cl,varlist = ls())

sample_inf_mc <- parSapply(cl, 1:N_mc,function(x){sample_inf()})
marginal_inf_mc <- parSapply(cl,1:N_mc,function(x){marginal_inf()})

stopCluster(cl)

mean(colMeans(x = sample_inf_mc))
mean(colMeans(x = marginal_inf_mc))

sd(colMeans(x = sample_inf_mc))
sd(colMeans(x = marginal_inf_mc))

library(matrixStats)

mean(colSds(x = sample_inf_mc))
mean(colSds(x = marginal_inf_mc))

sd(colSds(x = sample_inf_mc))
sd(colSds(x = marginal_inf_mc))
