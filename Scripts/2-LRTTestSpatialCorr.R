# 2-LRTTestSpatialCorr.R

# Script: LRT tests for indices and assess evidence for residual spatial autocorrelation
# Software authors: Max T. Eyre, Emanuele Giorgi, Peter Diggle (CHICAS, Lancaster University Medical School)

# Article title: A multivariate geostatistical framework for combining multiple indices of abundance for disease 
# vectors and reservoirs: A case study of rattiness in a low-income urban Brazilian community

# Description: In this script we fit four non-spatial models (all three indices, and three models each excluding a different index) 
# using the quasi-Monte Carlo method with the covariate selected for inclusion in 1-ExploreCovariates.R. We then conduct three likelihood 
# ratio tests to test for evidence against the null hypotheses that each index does not contribute to Rattiness. Finally we predict Uj 
# at each location and use a semivariogram to test for evidence of residual spatial correlation after controlling for covariate effects.

###-----------------------------------------------

rm(list=ls())
library(dplyr)
library(PrevMap)
library(geoR)
library(randtoolbox)
library(numDeriv)

#### DATA IN
rat <- read.csv("Data/1-PdLRattinessData.csv")
rat <- na.omit(rat)
ID <- create.ID.coords(rat,~X+Y)
U.halton <- qnorm(halton(1000,1))

rat$valley <- as.factor(rat$valley)
sp.trash <- function(x) max(0,x-90)
sp.trash <- Vectorize(sp.trash)

inst <- levels(rat$data_type)
n.inst <- 3 # 3 indices
coords <- unique(rat[,c("X","Y")])

# choose covariates
D.aux <- as.matrix(model.matrix(~-1 + elevation + dist_trash + sp.trash(dist_trash) + lc30_prop_veg + valley, data=rat))

# standardise covariate values
p <- ncol(D.aux)
N <- nrow(coords)
D <- matrix(NA,nrow=N,ncol=p)
for(i in 1:p) {
  D[,i] <- tapply(D.aux[,i],ID,max)
}
D <- scale(D)

p <- ncol(D)
n.par <- 2*n.inst+p
par.start <- rep(0,n.par)
par <- par.start

llik <- function(par,which="none") {
  if(all(which!=c("none","traps","signs","plates"))) {
    stop("'which' must be set to either 'none', 'signs', 'plates' or 'traps'")
  }
  # Set up parameters
  if(which!="none") {
    sigma <- exp(par[((n.inst+1):(n.inst+2))])
    beta <- par[(n.inst+3):(n.inst+2+p)]  
  } else {
    sigma <- exp(par[((n.inst+1):(2*n.inst))]) s
    beta <- par[(2*n.inst+1):(2*n.inst+p)]    
  }
  
  alpha <- par[1:n.inst]

  n.loc <- length(unique(ID))
  
  eta.R <- as.numeric(D%*%beta) 
  
  out <- 0
  
  for(j in 1:n.loc) { # please not that in the code we use j to denote location
    R.halton <- eta.R[j]+U.halton
    rat.j <- rat[ID==j,] 
    data_type.j <- rat.j$data_type 
    ind.plates.j <- which(data_type.j=="plates")
    ind.traps.j <- which(data_type.j=="traps")
    ind.signs.j <- which(data_type.j=="signs")
    
    out.signs.j <- 1
    out.traps.j <- 1
    out.plates.j <- 1
    
    if(length(ind.signs.j)>0) { 
      if(which=="signs") {
        eta.signs.j <- sapply(R.halton,function(R.j) alpha[1]) 
      } else {
        eta.signs.j <- sapply(R.halton,function(R.j) alpha[1]+sigma[1]*R.j) 
      }
      p.signs.j <- 1/(1+exp(- eta.signs.j)) 
      
      llik.signs.j <- 
        rat.j$outcome[ind.signs.j]*log(p.signs.j)+
        (1-rat.j$outcome[ind.signs.j])*log(1-p.signs.j)
      out.signs.j <- exp(llik.signs.j)
    }
    
    if(length(ind.traps.j) > 0) {
      if(which=="traps") {
        eta.traps.j <- sapply(R.halton,function(R.j) alpha[2]) 
      } else if(which!="none" & which!="traps") {
        eta.traps.j <- sapply(R.halton,function(R.j) alpha[2]+sigma[1]*R.j) 
      }  else {
        eta.traps.j <- sapply(R.halton,function(R.j) alpha[2]+sigma[2]*R.j) 
      }

      lambda.traps.j <- exp(eta.traps.j) 
      t.j <- rep(1,length(ind.traps.j)) 
      if(any(rat.j$offset_req[ind.traps.j]==1)) t.j[rat.j$offset_req[ind.traps.j]==1] <- 0.5
      
      llik.traps.j <- sapply(1:length(ind.traps.j), 
                             function(h) dbinom(rat.j$outcome[ind.traps.j][h],1, 
                                                prob = 1-exp(-t.j[h]*lambda.traps.j),
                                                log=TRUE))
      if(length(ind.traps.j) > 1) {
        out.traps.j <- apply(llik.traps.j,1,function(row) exp(sum(row))) 
      } else {
        out.traps.j <- exp(llik.traps.j)
      }
    }
    
    if(length(ind.plates.j) > 0) {
      if(which=="plates") {
        eta.plates.j <- sapply(R.halton,function(R.j) alpha[3]) 
      } else if(which!="none" & which!="plates") {
        eta.plates.j <- sapply(R.halton,function(R.j) alpha[3]+sigma[2]*R.j)
      } else {
        eta.plates.j <- sapply(R.halton,function(R.j) alpha[3]+sigma[3]*R.j)
      }
      
      q.plates.j <- 1/(1+exp(-eta.plates.j)) 
      llik.plates.j <- sapply(1:length(ind.plates.j), 
                              function(h) dbinom(rat.j$outcome[ind.plates.j][h],
                                                 rat.j$offset[ind.plates.j][h],
                                                 prob = q.plates.j,
                                                 log=TRUE))
      if(length(ind.plates.j) > 1) {
        out.plates.j <- apply(llik.plates.j,1,function(row) exp(sum(row)))
      } else {
        out.plates.j <- exp(llik.plates.j)
      }
    }
    
    out <- out+log(mean(exp(log(out.signs.j)+log(out.traps.j)+log(out.plates.j)))) 
  }
  return(out)
}

# Fit model
estim.full <- nlminb(par.start, 
                function(x) -llik(x,which="none"),
                control=list(trace=1))

estim.no.rat.for.traps <- nlminb(estim.full$par[-(n.inst+2)], 
                     function(x) -llik(x,which="traps"),
                     control=list(trace=1))
estim.no.rat.for.signs <- nlminb(estim.full$par[-(n.inst+1)], 
                                 function(x) -llik(x,which="signs"),
                                 control=list(trace=1))
estim.no.rat.for.plates <- nlminb(estim.full$par[-(2*n.inst)], 
                                 function(x) -llik(x,which="plates"),
                                 control=list(trace=1))

# Likelihood ratio test H0: an index does not contribute to Rattiness
(1-pchisq(-2*(estim.full$objective-estim.no.rat.for.traps$objective),1))/2
(1-pchisq(-2*(estim.full$objective-estim.no.rat.for.signs$objective),1))/2
(1-pchisq(-2*(estim.full$objective-estim.no.rat.for.plates$objective),1))/2

# save parameter estimates 
saveRDS(estim.full, "estim.full.RDS")
saveRDS(estim.no.rat.for.traps, "estim.no.rat.for.traps.RDS")
saveRDS(estim.no.rat.for.signs, "estim.no.rat.for.signs.RDS")
saveRDS(estim.no.rat.for.plates, "estim.no.rat.for.plates.RDS")

# --- Make spatial variogram from the full model ---
# Calculate mean Uj
compute.pred.mean <- function(par,j) {
  n.loc <- length(unique(ID))
  
  alpha <- par[1:n.inst]
  sigma <- exp(par[((n.inst+1):(2*n.inst))])
  beta <- par[(2*n.inst+1):(2*n.inst+p)]
  
  n.loc <- length(unique(ID))
  
  eta.R.j <- sum(D[j,]*beta)
  
  R.halton <- eta.R.j+U.halton
  R.halton <- U.halton
  rat.j <- rat[ID==j,]
  data_type.j <- rat.j$data_type
  ind.plates.j <- which(data_type.j=="plates")
  ind.traps.j <- which(data_type.j=="traps")
  ind.signs.j <- which(data_type.j=="signs")
  
  out.signs.j <- 1
  out.traps.j <- 1
  out.plates.j <- 1
  
  if(length(ind.signs.j)>0) {
    eta.signs.j <- sapply(R.halton,function(R.j) alpha[1]+sigma[1]*R.j)
    p.signs.j <- 1/(1+exp(- eta.signs.j))
    
    llik.signs.j <- 
      rat.j$outcome[ind.signs.j]*log(p.signs.j)+
      (1-rat.j$outcome[ind.signs.j])*log(1-p.signs.j)
    out.signs.j <- exp(llik.signs.j)
  }
  
  if(length(ind.traps.j) > 0) {
    eta.traps.j <- sapply(R.halton,function(R.j) alpha[2]+sigma[2]*R.j)
    lambda.traps.j <- exp(eta.traps.j)
    t.j <- rep(1,length(ind.traps.j))
    if(any(rat.j$offset_req[ind.traps.j]==1)) t.j[rat.j$offset_req[ind.traps.j]==1] <- 0.5
    llik.traps.j <- sapply(1:length(ind.traps.j), 
                           function(h) dbinom(rat.j$outcome[ind.traps.j][h],1,
                                              prob = 1-exp(-t.j[h]*lambda.traps.j),
                                              log=TRUE))
    if(length(ind.traps.j) > 1) {
      out.traps.j <- apply(llik.traps.j,1,function(row) exp(sum(row)))
    } else {
      out.traps.j <- exp(llik.traps.j)
    }
  }
  
  if(length(ind.plates.j) > 0) {
    eta.plates.j <- sapply(R.halton,function(R.j) alpha[3]+sigma[3]*R.j)
    q.plates.j <- 1/(1+exp(-eta.plates.j))
    llik.plates.j <- sapply(1:length(ind.plates.j), 
                            function(h) dbinom(rat.j$outcome[ind.plates.j][h],
                                               rat.j$offset[ind.plates.j][h],
                                               prob = q.plates.j,
                                               log=TRUE))
    if(length(ind.plates.j) > 1) {
      out.plates.j <- apply(llik.plates.j,1,function(row) exp(sum(row)))
    } else {
      out.plates.j <- exp(llik.plates.j)
    }
  }
  den <- mean(exp(log(out.signs.j)+log(out.traps.j)+log(out.plates.j)))
  num <- mean(U.halton*exp(log(out.signs.j)+log(out.traps.j)+log(out.plates.j)))
  num/den
}
n.loc <- length(unique(ID))
estim <- estim.full
U.pred.mean <- sapply(1:n.loc,function(j) compute.pred.mean(estim$par,j))

# Merge predicted Uj values with main dataframe, then get unique locations only
U.data <- data.frame(X=unique(rat[,c("X","Y")])[,1],
                     Y=unique(rat[,c("X","Y")])[,2],
                     U=U.pred.mean)

ratUj <- left_join(rat,U.data, by=c("X","Y"))
ratUj <- ratUj[,c("valley","X","Y","elevation","U","dist_sewer","dist_trash",
                  "lc10_prop_soil","lc10_prop_veg","lc10_prop_pave","lc_perv_10m",
                  "lc30_prop_soil","lc30_prop_veg","lc30_prop_pave","lc_perv_30m","rel_elev","dist_pri_sewer")]
ratUj <- unique(ratUj)

#### Exploratory analysis - test for residual spatial correlation

# Histogram of Uj
hist(ratUj$U)

# Spatial correlation in data
U.data <- as.data.frame(U.data)
ratUj <- as.data.frame(ratUj)
vari.no.cov <- variogram(ratUj,~U,~X+Y,uvec=seq(0,200,length=15))
env <- variog.mc.env(geodata=as.geodata(U.data),obj.variog=vari.no.cov,nsim=1000)
plot(vari.no.cov,envelope=env, xlab="Distance (m)", ylab = "Semivariance")
