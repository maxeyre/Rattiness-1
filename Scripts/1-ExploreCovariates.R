# 1-ExploreCovariates.R

# Script: Exploring covariate relationships with estimated rattiness
# Software authors: Max T. Eyre, Emanuele Giorgi, Peter Diggle (CHICAS, Lancaster University Medical School)

# Article title: A multivariate geostatistical framework for combining multiple indices of abundance for disease 
# vectors and reservoirs: A case study of rattiness in a low-income urban Brazilian community

# Description: In this script we explore covariate relationships by fitting a model (using the quasi-Monte Carlo method) 
# without covariates which assumes no spatial correlation. We then predict Uj (=Rj) at each location for which data was observed and
# assess the data for evidence of spatial correlation. Relationships of each covariate can then be assess by plotting these values 
# against covariate values

###-----------------------------------------------

rm(list=ls())
library(dplyr)
library(PrevMap)
library(randtoolbox)
library(numDeriv)

#### DATA IN
rat <- read.csv("Data/1-PdLRattinessData.csv")
rat <- na.omit(rat)
ID <- create.ID.coords(rat,~X+Y)
U.halton <- qnorm(halton(1000,1))

# standardise covariate values
inst <- levels(rat$data_type)

n.inst <- 3 # 3 indices
n.par <- 2*n.inst
par.start <- rep(0,n.par) # first guess at parameters


llik <- function(par) {
  # Set up parameters
  alpha <- par[1:n.inst]
  sigma <- exp(par[((n.inst+1):(2*n.inst))])
  n.loc <- length(unique(ID))
  
  out <- 0
  for(j in 1:n.loc) { # please not that in the code we use j to denote location
    
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
    
    out <- out+log(mean(exp(log(out.signs.j)+log(out.traps.j)+log(out.plates.j)))) 
  }
  return(out)
}

# Fit model
estim <- nlminb(par.start, 
                function(x) -llik(x),
                control=list(trace=1))

# Calculate mean Uj
compute.pred.mean <- function(par,j) {
  n.loc <- length(unique(ID))
  
  alpha <- par[1:n.inst]
  sigma <- exp(par[((n.inst+1):(2*n.inst))])
  n.loc <- length(unique(ID))

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

#### Exploratory analysis

# Histogram of Uj
hist(ratUj$U)

# You can now plot Uj against covariates of interest (using ratUj) to explore relationships and decide how they 
# should be included in the model. For example:

# elevation
x.elevation <- cut(ratUj$elevation,quantile(ratUj$elevation,seq(0,1,0.12)),
                   include.lowest=TRUE)
elevation.bin <- x.elevation
x.elevation.mean <- tapply(ratUj$elevation,elevation.bin,mean)
prop.elevation <- tapply(ratUj$U,elevation.bin,
                         function(x) mean(x,na.rm=TRUE))
n.obs.elevation <- tapply(elevation.bin,elevation.bin,length)
plot(x.elevation.mean,prop.elevation,cex=n.obs.elevation*0.06, main="",xlab="Mean elevation (m)",
     ylab="Estimated Rj", xlim=c(25,50)) 
