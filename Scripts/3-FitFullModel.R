# 3-FitFullModel.R

# Script: Fitting the full model
# Software authors: Max T. Eyre, Emanuele Giorgi, Peter Diggle (CHICAS, Lancaster University Medical School)

# Article title: A multivariate geostatistical framework for combining multiple indices of abundance for disease 
# vectors and reservoirs: A case study of rattiness in a low-income urban Brazilian community

# Description: In this script we fit the full geostatistical model. Please see the published article for a 
# detailed explanation of the methodology.

###-----------------------------------------------

rm(list=ls())
library(dplyr)
library(PrevMap)
set.seed(590)

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
D.aux <- as.matrix(model.matrix(~-1 + elevation + dist_trash + sp.trash(dist_trash) + lc30_prop_veg + valley, data=rat))

p <- ncol(D.aux)
N <- nrow(coords)
D <- matrix(NA,nrow=N,ncol=p)
for(i in 1:p) {
  D[,i] <- tapply(D.aux[,i],ID,max)
}
D <- scale(D)

ind3 <- which(rat$data_type=="plates")
ind1 <- which(rat$data_type=="signs")
ind2 <- which(rat$data_type=="traps")

# initial parameter guesses
estim.par <- readRDS("estim.full.RDS")
par0 <- c(estim.par$par, 3, 0.5)

k <- 0
while(k <= 2) {
  
  k <- k + 1

alpha1.0 <- par0[1]
alpha2.0 <- par0[2]
alpha3.0 <- par0[3]
sigma1.0 <- exp(par0[4])
sigma2.0 <- exp(par0[5])
sigma3.0 <- exp(par0[6])
beta0 <- par0[7:(p+6)]
phi0 <- exp(par0[p+7])
psi0 <- exp(par0[p+8])/(1+exp(par0[p+8]))

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

if(k==1) {
  n.sim <- 10000
  burnin <- 2000
  thin <- 8  
} else {
  n.sim <- 110000
  burnin <- 10000
  thin <- 10
}

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

log.integrand <- function(R,val) {
  eta1 <- val$alpha1+val$sigma1*R[ID[ind1]]
  llik1 <- sum(y[ind1]*eta1-offset[ind1]*log(1+exp(eta1)))
  
  lambda2 <- exp(val$alpha2+val$sigma2*R[ID[ind2]])
  prob2 <- 1-exp(-offset[ind2]*lambda2)
  llik2 <- sum(y[ind2]*log(prob2/(1-prob2))+log(1-prob2))
  
  eta3 <- val$alpha3+val$sigma3*R[ID[ind3]]
  llik3 <- sum(y[ind3]*eta3-offset[ind3]*log(1+exp(eta3)))
  
  diff.R <- R-val$mu
  out <- as.numeric(-0.5*(val$log.det.Sigma+t(diff.R)%*%val$Sigma.inv%*%diff.R))+
    llik1+llik2+llik3
}

compute.log.f <- function(par,ldetR=NA,R.inv=NA) {
  val <- list()
  val$alpha1 <- par[1]
  val$alpha2 <- par[2]
  val$alpha3 <- par[3]
  
  val$sigma1 <- exp(par[4])
  val$sigma2 <- exp(par[5])
  val$sigma3 <- exp(par[6])
  
  beta <- par[7:(p+6)]
  val$mu <- as.numeric(D%*%beta)
  
  phi <- exp(par[p+7])
  psi <- exp(par[p+8])/(1+exp(par[p+8]))
  
  Sigma <- as.matrix(psi*exp(-U/phi))
  val$Sigma.inv <- solve(Sigma)
  val$log.det.Sigma <- determinant(Sigma)$modulus
  
  sapply(1:(dim(sim)[1]),function(i) log.integrand(sim[i,],val))
}

par0 <- c(alpha1.0,alpha2.0,alpha3.0,log(sigma1.0),log(sigma2.0),log(sigma3.0),
          beta0,log(phi0),log(psi0/(1-psi0)))
log.f.tilde <- compute.log.f(par0)

MC.log.lik <- function(par) {
  log(mean(exp(compute.log.f(par)-log.f.tilde)))
}

MC.log.lik(par0)

estim.par <- nlminb(par0,
                    function(x) -MC.log.lik(x),
                    control=list(trace=1))
par0 <- estim.par$par

}
# your parameter estimates are now in estim.par

# We recommend repeated re-fitting model with new parameter estimate plug-in until 
# euclid norm of relative difference between two consecutive parameter estimates falls below 
# a chosen value. Parameter estimates reported in the published article were estimated 
# following this method.

saveRDS(estim.par,"estim.par.RDS")
