rm(list=ls());gc()

# params
D <- 25 # this many dwellings
A <- 20 # this many aquatic habitats
H <- rpois(n = D,lambda = 5) # vector of humans at each house
alpha <- 1/3 # feeding rate on humans
pfpr <- 0.25 # 25% prevalence

# vector for dwellings, each element is d(p)
d_p <- lapply(X = H,FUN = function(x){
  p <- rexp(n = x,rate = 1)
  p/sum(p)
})

# vector for dwellings, each element is d(c)
d_c <- lapply(X = H,FUN = function(x){
  sapply(X = 1:x,FUN = function(h){
    if(runif(1) < pfpr){
      rbeta(n = 1,shape1 = 5,shape2 = 10)
    } else {
      0
    }
  })
})

# vector X "net infectiousness to mosquitos"
X <- mapply(function(p,c){
  c %*% p
},p = d_p,c = d_c)

# simulate some points
library(spatstat)

xlim <- c(0,1); ylim <- c(0,1)
pts_spatstat <- rmpoint(n = c(D,A),win = owin(xrange = xlim,yrange = ylim))
pts <- data.frame(x=as.numeric(pts_spatstat$x),y=as.numeric(pts_spatstat$y),marks=as.integer(pts_spatstat$marks))
xy_d <- as.matrix(pts[which(pts$marks==1),c("x","y")]) # dwellings
xy_a <- as.matrix(pts[which(pts$marks==2),c("x","y")]) # habitats

# calculate sigma for all a
sigma_a <- rep(0,A)

pb <- txtProgressBar(1,A)
for(a in 1:A){

  dist <- as.matrix(dist(x=rbind(xy_a[a,],xy_d)))[1,2:(nrow(xy_d)+1)] # vector of distance from a to all D
  sigma_a[a] <- min(dist) # sigma(a)

  setTxtProgressBar(pb,a)
}

# make the bite dispersal matrix K (each row will be kappa for aquatic habitat a)
K <- matrix(0,nrow=A,ncol=D)

for(a in 1:A){

  dist <- as.matrix(dist(x=rbind(xy_a[a,],xy_d)))[1,2:(nrow(xy_d)+1)] # vector of distance from a to all D
  kappa <- sapply(dist,function(x){
    dnorm(x,mean=0,sd=sigma_a[a])
  })
  K[a,] <- kappa/sum(kappa)
}

# force of infection on mosquitos
lambdaV <- alpha * (K %*% X)

# plot FOI on mosquitos against Psi
psi_grid <- expand.grid(
  x=seq(xlim[1]-0.05,xlim[2]+0.05,by=0.01),
  y=seq(ylim[1]-0.05,ylim[2]+0.05,by=0.01)
)

library(foreach)
library(iterators)
library(doParallel)

cl <- makeCluster(4)
registerDoParallel(cl)

habitats <- data.frame(x=xy_a[,"x"],y=xy_a[,"y"],sigma=sigma_a,lambda=lambdaV)
dwellings <- data.frame(x=xy_d[,"x"],y=xy_d[,"y"])

psi <- foreach(xy = iter(psi_grid,by="row"),.combine = "rbind",.inorder = TRUE,.verbose = F) %:%
  foreach(hab = iter(habitats,by = "row"),.combine = "+") %dopar% {
    dist <- as.matrix(dist(x = rbind(as.vector(xy),c(hab$x,hab$y))))[1,2]
    psi <- dnorm(dist,mean=0,sd=hab$sigma)
    psi
  }

stopCluster(cl)

psi_grid$psi <- as.vector(psi)

library(ggplot2)
library(viridis)

library(tikzDevice)

tikz(file = "/Users/slwu89/Desktop/git/mosymodel/tikzplots/rf_LambdaV.tex",width = 5,height = 4)

ggplot() +
  geom_raster(aes(x=x,y=y,fill=psi),data=psi_grid) +
  stat_contour(aes(x=x,y=y,z=psi),colour=grey(0.75,0.75),size=0.25,data = psi_grid,geom = "contour") +
  scale_fill_viridis() +
  geom_point(aes(x=x,y=y,colour=lambda),data=habitats) +
  scale_color_viridis(begin = 0.5,end = 1,option = "A") +
  theme_bw() +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank())

dev.off()
