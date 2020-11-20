# 4-Prediction.R

# Script: Spatial prediction of Rattiness
# Software authors: Max T. Eyre, Emanuele Giorgi, Peter Diggle (CHICAS, Lancaster University Medical School)

# Article title: A multivariate geostatistical framework for combining multiple indices of abundance for disease 
# vectors and reservoirs: A case study of rattiness in a low-income urban Brazilian community

# Description: In this script we make predictions at unobserved locations within the study area at the points in
# the dataset 2-PdLPredictionGrid.csv

###-----------------------------------------------

rm(list=ls())
library(dplyr)
library(PrevMap)

#### DATA IN
rat <- read.csv("Data/1-PdLRattinessData.csv")
rat <- na.omit(rat)
ID <- create.ID.coords(rat,~X + Y)
coords <- unique(rat[,c("X","Y")])
U <- dist(coords)

rat$valley <- as.factor(rat$valley)
sp.trash <- function(x) max(0,x-90)
sp.trash <- Vectorize(sp.trash)

# choose covariates
D.aux <- as.matrix(model.matrix(~-1+ elevation + dist_trash + sp.trash(dist_trash) + lc30_prop_veg + valley, data=rat))

p <- ncol(D.aux)
N <- nrow(coords)
D <- matrix(NA,nrow=N,ncol=p)
for(i in 1:p) {
  D[,i] <- tapply(D.aux[,i],ID,max)
}

# standardise
D.orig <- D
mean.D <- apply(D.orig,2,mean)
sd.D <- apply(D.orig,2,sd)
D <- sapply(1:ncol(D.orig),function(i) (D.orig[,i]-mean.D[i])/sd.D[i])

ind3 <- which(rat$data_type=="plates")
ind1 <- which(rat$data_type=="signs")
ind2 <- which(rat$data_type=="traps")

# read in the estimated parameter values from 3-FitFullModel.R here
estim.par <- readRDS("estim.par.RDS")
alpha1.0 <- estim.par$par[1]
alpha2.0 <- estim.par$par[2]
alpha3.0 <- estim.par$par[3]
beta0 <- estim.par$par[7:(p+6)]
sigma1.0 <- exp(estim.par$par[4])
sigma2.0 <- exp(estim.par$par[5])
sigma3.0 <- exp(estim.par$par[6])
phi0 <- exp(estim.par$par[p+7])
psi0 <- exp(estim.par$par[p+8])/(1+exp(estim.par$par[p+8]))

Sigma0 <- as.matrix(psi0*exp(-U/phi0))
Sigma0.inv <- solve(Sigma0)

mu0 <- as.numeric(D%*%beta0)
y <- rat$outcome
offset <- rat$offset

integrand <- function(R) {
  
  eta1 <- alpha1.0+sigma1.0*R[ID[ind1]]
  llik1 <- sum(y[ind1]*eta1-offset[ind1]*log(1+exp(eta1)))
  
  lambda2 <- exp(alpha2.0+sigma2.0*R[ID[ind2]])
  prob2 <- 1-exp(-offset[ind2]*lambda2)
  llik2 <- sum(y[ind2]*log(prob2/(1-prob2))+log(1-prob2))
  
  eta3 <- alpha3.0+sigma3.0*R[ID[ind3]]
  llik3 <- sum(y[ind3]*eta3-offset[ind3]*log(1+exp(eta3)))
  
  diff.R <- R-mu0
  out <- as.numeric(-0.5*t(diff.R)%*%Sigma0.inv%*%diff.R)+
    llik1+llik2+llik3
  as.numeric(out)
}

ID1 <- sort(unique(ID[ind1]))
ID2 <- sort(unique(ID[ind2]))
ID3 <- sort(unique(ID[ind3]))
grad.integrand <- function(R) {
  der.tot <- rep(0,N)
  
  eta1 <- alpha1.0+sigma1.0*R[ID[ind1]]
  der.tot[ID1] <- 
    tapply((y[ind1]-offset[ind1]*exp(eta1)/(1+exp(eta1)))*sigma1.0,ID[ind1],sum)
  
  lambda2 <- exp(alpha2.0+sigma2.0*R[ID[ind2]])
  prob2 <- 1-exp(-offset[ind2]*lambda2)
  der.prob2 <- offset[ind2]*exp(-offset[ind2]*lambda2)*lambda2*sigma2.0
  der.tot[ID2] <- der.tot[ID2]+
    tapply((y[ind2]/(prob2*(1-prob2))-1/(1-prob2))*der.prob2,ID[ind2],sum)
  
  eta3 <- alpha3.0+sigma3.0*R[ID[ind3]]
  der.tot[ID3] <- der.tot[ID3]+
    tapply((y[ind3]-offset[ind3]*exp(eta3)/(1+exp(eta3)))*sigma3.0,ID[ind3],sum)
  
  diff.R <- R-mu0
  out <- -Sigma0.inv%*%diff.R+der.tot
  as.numeric(out)
}

hessian.integrand <- function(R) {
  hess.tot <- rep(0,N)
  
  eta1 <- alpha1.0+sigma1.0*R[ID[ind1]]
  hess.tot[ID1] <- hess.tot[ID1]+
    tapply(-(offset[ind1]*exp(eta1)/((1+exp(eta1))^2))*sigma1.0^2,ID[ind1],sum)
  
  lambda2 <- exp(alpha2.0+sigma2.0*R[ID[ind2]])
  prob2 <- 1-exp(-offset[ind2]*lambda2)
  der.prob2 <- offset[ind2]*exp(-offset[ind2]*lambda2)*lambda2*sigma2.0
  der2.prob2 <- -((offset[ind2])^2)*exp(-offset[ind2]*lambda2)*(lambda2*sigma2.0)^2+
    offset[ind2]*exp(-offset[ind2]*lambda2)*lambda2*(sigma2.0)^2
  hess.tot[ID2] <- hess.tot[ID2]+
    tapply((y[ind2]/(prob2*(1-prob2))-1/(1-prob2))*der2.prob2+
             (y[ind2]*((2*prob2-1)/((prob2*(1-prob2))^2))-1/(1-prob2)^2)*(der.prob2^2),ID[ind2],sum)
  
  eta3 <- alpha3.0+sigma3.0*R[ID[ind3]]
  hess.tot[ID3] <- hess.tot[ID3]+
    tapply(-(offset[ind3]*exp(eta3)/((1+exp(eta3))^2))*sigma3.0^2,ID[ind3],sum)
  
  out <- -Sigma0.inv
  diag(out) <- diag(out)+hess.tot
  out
}

estim <- nlminb(start=rep(0,N),
                function(x) -integrand(x),
                function(x) -grad.integrand(x),
                function(x) -hessian.integrand(x),
                control=list(trace=1))
H <- hessian.integrand(estim$par)

Sigma.sroot <- t(chol(solve(-H)))
A <- solve(Sigma.sroot)
Sigma.W.inv <- solve(A%*%Sigma0%*%t(A))
mu.W <- as.numeric(A%*%(mu0-estim$par))

cond.dens.W <- function(W,R) {
  eta1 <- alpha1.0+sigma1.0*R[ID[ind1]]
  llik1 <- sum(y[ind1]*eta1-offset[ind1]*log(1+exp(eta1)))
  
  lambda2 <- exp(alpha2.0+sigma2.0*R[ID[ind2]])
  prob2 <- 1-exp(-offset[ind2]*lambda2)
  llik2 <- sum(y[ind2]*log(prob2/(1-prob2))+log(1-prob2))
  
  eta3 <- alpha3.0+sigma3.0*R[ID[ind3]]
  llik3 <- sum(y[ind3]*eta3-offset[ind3]*log(1+exp(eta3)))
  
  diff.W <- W-mu.W
  -0.5*as.numeric(t(diff.W)%*%Sigma.W.inv%*%diff.W)+
    llik1+llik2+llik3
}

lang.grad <- function(W,R) {
  der.tot <- rep(0,N)
  
  eta1 <- alpha1.0+sigma1.0*R[ID[ind1]]
  der.tot[ID1] <- 
    tapply((y[ind1]-offset[ind1]*exp(eta1)/(1+exp(eta1)))*sigma1.0,ID[ind1],sum)
  
  lambda2 <- exp(alpha2.0+sigma2.0*R[ID[ind2]])
  prob2 <- 1-exp(-offset[ind2]*lambda2)
  der.prob2 <- offset[ind2]*exp(-offset[ind2]*lambda2)*lambda2*sigma2.0
  der.tot[ID2] <- der.tot[ID2]+
    tapply((y[ind2]/(prob2*(1-prob2))-1/(1-prob2))*der.prob2,ID[ind2],sum)
  
  eta3 <- alpha3.0+sigma3.0*R[ID[ind3]]
  der.tot[ID3] <- der.tot[ID3]+
    tapply((y[ind3]-offset[ind3]*exp(eta3)/(1+exp(eta3)))*sigma3.0,ID[ind3],sum)
  
  diff.W <- W-mu.W
  
  as.numeric(-Sigma.W.inv%*%diff.W+
               t(Sigma.sroot)%*%der.tot)
}


h <- 1.65/(N^(1/6))
n.sim <- 110000
burnin <- 10000
thin <- 10
c1.h <- 0.001
c2.h <- 0.0001
W.curr <- rep(0,N)

R.curr <- as.numeric(Sigma.sroot%*%W.curr+estim$par)
mean.curr <- as.numeric(W.curr + (h^2/2)*lang.grad(W.curr,R.curr))
lp.curr <- cond.dens.W(W.curr,R.curr)
acc <- 0
n.samples <- (n.sim-burnin)/thin
sim <- matrix(NA,nrow=n.samples,ncol=N)

h.vec <- rep(NA,n.sim)
for(i in 1:n.sim) {
  W.prop <- mean.curr+h*rnorm(N)
  R.prop <-  as.numeric(Sigma.sroot%*%W.prop+estim$par)
  mean.prop <- as.numeric(W.prop + (h^2/2)*lang.grad(W.prop,R.prop))
  lp.prop <- cond.dens.W(W.prop,R.prop)
  
  dprop.curr <- -sum((W.prop-mean.curr)^2)/(2*(h^2))
  dprop.prop <- -sum((W.curr-mean.prop)^2)/(2*(h^2))
  
  log.prob <- lp.prop+dprop.prop-lp.curr-dprop.curr
  
  if(log(runif(1)) < log.prob) {
    acc <- acc+1
    W.curr <- W.prop
    R.curr <- R.prop
    lp.curr <- lp.prop
    mean.curr <- mean.prop
  }
  
  if( i > burnin & (i-burnin)%%thin==0) {
    sim[(i-burnin)/thin,] <- R.curr
  }
  
  h.vec[i] <- h <- max(0,h + c1.h*i^(-c2.h)*(acc/i-0.57))
  cat("Iteration",i,"out of",n.sim,"\r")
  flush.console()
}

acf.plot <- acf(sim[,1],plot=FALSE)
plot(acf.plot$lag,acf.plot$acf,type="l",xlab="lag",ylab="autocorrelation",
     ylim=c(-0.1,1),main="Autocorrelogram of the simulated samples")
for(i in 2:ncol(sim)) {
  acf.plot <- acf(sim[,i],plot=FALSE)
  lines(acf.plot$lag,acf.plot$acf)
}
abline(h=0,lty="dashed",col=2)

# --- Prediction ---
library(sf)
library(tmap)
library(pdist)
library(splancs)
library(raster)

points <- read.csv("Data/2-PdLPredictionGrid.csv")
points$valley <- as.factor(points$valley)
grid.pred <- points[,c("X","Y")]
n.pred <- nrow(grid.pred)

# For rattiness
D.pred <- as.matrix(model.matrix(~-1+ Z + dist_trash + sp.trash(dist_trash) + lc30_prop_veg + valley, data=points)) # change covariates as necessary
D.pred <- sapply(1:ncol(D.pred),function(i) (D.pred[,i]-mean.D[i])/sd.D[i])
mu.pred <- as.numeric(D.pred%*%beta0)
U.pred <- as.matrix(pdist(grid.pred,coords))
C <- psi0*exp(-U.pred/phi0)
A <- C%*%Sigma0.inv
R.pred.cond.mean <- sapply(1:n.samples,function(i) mu.pred+A%*%(sim[i,]-mu0))
R.pred.hat <- apply(R.pred.cond.mean,1,mean)
R.pred.sd <- sqrt(psi0-apply(A*C,1,sum))

# For S(x)
D.pred <- matrix(0,ncol=p,nrow=n.pred)
mu.pred <- as.numeric(D.pred%*%beta0)
U.pred <- as.matrix(pdist(grid.pred,coords))
C <- psi0*exp(-U.pred/phi0)
A <- C%*%Sigma0.inv
S.pred.cond.mean <- sapply(1:n.samples,function(i) mu.pred+A%*%(sim[i,]-mu0))
S.pred.hat <- apply(S.pred.cond.mean,1,mean)

# -- Make into a raster --
points1 <- read_csv("Data/Outlines/outline1_points.csv")
poly1 <- data.frame(x=points1$X, y=points1$Y)
points2 <- read_csv("Data/Outlines/outline2_points.csv")
poly2 <- data.frame(x=points2$X, y=points2$Y)
res <-5
grid.pred1 <- gridpts(as.matrix(poly1),xs=res,ys=res)
grid.pred2 <- gridpts(as.matrix(poly2),xs=res,ys=res)

# Rattiness
r.R.hat_1 <- rasterFromXYZ(cbind(grid.pred1,R.pred.hat[1:nrow(grid.pred1)]))
r.R.hat_2 <- rasterFromXYZ(cbind(grid.pred2,R.pred.hat[(nrow(grid.pred1)+1):nrow(grid.pred)]))
crs(r.R.hat_1) <- CRS("+proj=utm +zone=24 +south +datum=WGS84 +units=m +no_defs")
crs(r.R.hat_2) <- CRS("+proj=utm +zone=24 +south +datum=WGS84 +units=m +no_defs")
plot(r.R.hat_1)
plot(r.R.hat_2)

r.R.sd_1 <- rasterFromXYZ(cbind(grid.pred1,R.pred.sd[1:nrow(grid.pred1)]))
r.R.sd_2 <- rasterFromXYZ(cbind(grid.pred2,R.pred.sd[(nrow(grid.pred1)+1):nrow(grid.pred)]))
crs(r.R.sd_1) <- CRS("+proj=utm +zone=24 +south +datum=WGS84 +units=m +no_defs")
crs(r.R.sd_2) <- CRS("+proj=utm +zone=24 +south +datum=WGS84 +units=m +no_defs")
plot(r.R.sd_1)
plot(r.R.sd_2)

# S(x)
r.S.hat_1 <- rasterFromXYZ(cbind(grid.pred1,S.pred.hat[1:nrow(grid.pred1)]))
r.S.hat_2 <- rasterFromXYZ(cbind(grid.pred2,S.pred.hat[(nrow(grid.pred1)+1):nrow(grid.pred)]))
crs(r.S.hat_1) <- CRS("+proj=utm +zone=24 +south +datum=WGS84 +units=m +no_defs")
crs(r.S.hat_2) <- CRS("+proj=utm +zone=24 +south +datum=WGS84 +units=m +no_defs")
plot(r.S.hat_1)
plot(r.S.hat_2)

