
#####################################################
### Fit NB, POLN, POWB, POGG to data from POISSON ###
#####################################################

############
## FIT = NB

rm(list=ls())

args <- commandArgs(trailingOnly = TRUE)

out = args[1] #working directory
cat(",out=",out)
cores = as.numeric(args[2]) #number of cores to run in parallel
cat(",cores=",cores)

library(fitdistrplus)
library(pracma)
library(MASS)
library(flexsurv)
library(foreach)
library(doParallel)

load("data_POIS.RData")
source("distfunctions.R")
B.gam.pois = 1000


registerDoParallel(cores=35)
res.sim.gam.pois <- list()
res.sim.gam.pois <- foreach(i=1:B.gam.pois) %dopar%{
    fit.gampois = fitdistrplus::fitdist(mat[i,],"nbinom")
    res.sim.gam.pois[[i]] = c(fit.gampois$estimate[2],fit.gampois$estimate[1],fit.gampois$sd[2],fit.gampois$sd[1],fit.gampois$aic,fit.gampois$loglik)
}

name<-"Simulation_pois_gampois.RData"
setwd(out)
save(res.sim.gam.pois,file = name)

############
## FIT = LN

rm(list=ls())

args <- commandArgs(trailingOnly = TRUE)

out = args[1] #working directory
cat(",out=",out)
cores = as.numeric(args[2]) #number of cores to run in parallel
cat(",cores=",cores)

library(fitdistrplus)
library(pracma)
library(MASS)
library(flexsurv)
library(foreach)
library(doParallel)

load("data_POIS.RData")
source("distfunctions.R")
B.gam.pois = 1000

registerDoParallel(cores=35)
res.sim.log.pois <- list()
res.sim.log.pois <- foreach(i=1:B.gam.pois) %dopar%{
    mlog = log(mean(mat[i,])^2/sqrt(var(mat[i,]) + mean(mat[i,])^2))
    slog = sqrt(log(1+(var(mat[i,])/mean(mat[i,])^2)))
    if(slog<0) slog=1
    fit.logpois=try(fitdistrplus::fitdist(mat[i,],"logpois",start=list(sigma=slog,mu=mlog),lower=c(0,-Inf),optim.method="L-BFGS-B",discrete=T),silent=T)
    while(class(fit.logpois)=="try-error"){
      mlog=mlog+runif(1,-0.1,0.1); slog=slog+runif(1,-0.1,0.1)
      if(slog<0) slog=1
      fit.logpois=try(fitdistrplus::fitdist(mat[i,],"logpois",start=list(sigma=slog,mu=mlog),lower=c(0,-Inf),optim.method="L-BFGS-B",discrete=T),silent=T)
    }    
    res.sim.log.pois[[i]] = c(fit.logpois$estimate[2],fit.logpois$estimate[1],fit.logpois$sd[2],fit.logpois$sd[1],fit.logpois$aic,fit.logpois$loglik)
}

name<-"Simulation_pois_logpois.RData"
setwd(out)
save(res.sim.log.pois,file = name)

############
## FIT = WB

rm(list=ls())

args <- commandArgs(trailingOnly = TRUE)

out = args[1] #working directory
cat(",out=",out)
cores = as.numeric(args[2]) #number of cores to run in parallel
cat(",cores=",cores)

library(fitdistrplus)
library(pracma)
library(MASS)
library(flexsurv)
library(foreach)
library(doParallel)

load("data_POIS.RData")
source("distfunctions.R")
B.gam.pois = 1000

fun.k = function(x) v + mu^2 - ((mu^2*(gamma(1+2/x)))/(gamma(1+1/x)^2))

registerDoParallel(cores=35)
res.sim.bull.pois <- list()
res.sim.bull.pois <- foreach(i=1:B.gam.pois) %dopar%{
    mu = mean(mat[i,]); v = var(mat[i,])
    k_w = uniroot(fun.k,c(0.1,100))$root
    l_w = mean(mat[i,])/(gamma(1+1/k_w))
    if(k_w<0) k_w=1; if(l_w<0) l_w=1
    fit.bullpois = try(fitdistrplus::fitdist(mat[i,],"bullpois",start=list(p=k_w,l=l_w),lower=c(0,0),optim.method="L-BFGS-B",discrete=T),silent=T)
    while(class(fit.bullpois)=="try-error"){
      k_w=k_w+runif(1,-0.1,0.1); l_w=l_w+runif(1,-0.1,0.1)
      if(k_w<0) k_w=1; if(l_w<0) l_w=1
      fit.bullpois = try(fitdistrplus::fitdist(mat[i,],"bullpois",start=list(p=k_w,l=l_w),lower=c(0,0),optim.method="L-BFGS-B",discrete=T),silent=T)
    }   
    res.sim.bull.pois[[i]] = c(fit.bullpois$estimate[1],fit.bullpois$estimate[2],fit.bullpois$sd[1],fit.bullpois$sd[2],fit.bullpois$aic,fit.bullpois$loglik)
}

name<-"Simulation_pois_bullpois.RData"
setwd(out)
save(res.sim.bull.pois,file = name)

############
## FIT = GG

rm(list=ls())

args <- commandArgs(trailingOnly = TRUE)

out = args[1] #working directory
cat(",out=",out)
cores = as.numeric(args[2]) #number of cores to run in parallel
cat(",cores=",cores)

library(fitdistrplus)
library(pracma)
library(MASS)
library(flexsurv)
library(foreach)
library(doParallel)

load("data_POIS.RData")
source("distfunctions.R")
B.gam.pois = 1000

registerDoParallel(cores=35)
res.sim.gengam.pois <- list()
res.sim.gengam.pois <- foreach(i=1:B.gam.pois) %dopar%{
      a.g = var(mat[i,])/mean(mat[i,])
      p.g = mean(mat[i,])^2/var(mat[i,])
      d.g = 1
      if(a.g<0) a.g=0; if(p.g<0) p.g=0
      fit.gengampois = try(fitdistrplus::fitdist(mat[i,], "gengampois", start = list(a = a.g, p = p.g, d = d.g),lower=c(0,0,0),optim.method="L-BFGS-B",discrete=T),silent=T)
      while(class(fit.gengampois)=="try-error"){
        a.g=a.g+runif(1,-0.1,0.1); p.g=p.g+runif(1,-0.1,0.1); d.g = d.g+runif(1,-0.1,0.1)
        if(a.g<0) a.g=0; if(p.g<0) p.g=0; if(d.g<0) d.g=1
        fit.gengampois = try(fitdistrplus::fitdist(mat[i,], "gengampois", start = list(a = a.g, p = p.g, d = d.g),discrete=T, lower=c(0,0,0),optim.method="L-BFGS-B"),silent=T)
      }   
      res.sim.gengam.pois[[i]] = c(fit.gengampois$estimate[1], fit.gengampois$estimate[2],fit.gengampois$estimate[3],fit.gengampois$sd[1],fit.gengampois$sd[2],
                                fit.gengampois$sd[3], fit.gengampois$aic,fit.gengampois$loglik)
}

name<-"Simulation_pois_gengampois.RData"
setwd(out)
save(res.sim.gengam.pois,file = name)
