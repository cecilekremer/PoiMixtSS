
##################################################
### Fit POWB to data from NB, POLN, POWB, POGG ###
##################################################

############
## DATA = NB

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

load("data_NB.RData")
source("distfunctions.R")
B.gam.pois = 1000

fun.k = function(x) v + mu^2 - ((mu^2*(gamma(1+2/x)))/(gamma(1+1/x)^2))

registerDoParallel(cores=35)
res.sim.bull.pois <- list()

res.sim.bull.pois <- foreach(i=1:B.gam.pois) %dopar%{
    mu = mean(mat[i,]); v = var(mat[i,])
    k_w = uniroot(fun.k,c(0.1,100))
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

name<-"Simulation_gampois_bullpois.RData"
setwd(out)
save(res.sim.bull.pois,file = name)

############
## DATA = LN

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

load("data_LN.RData")
source("distfunctions.R")
B.gam.pois = 1000

fun.k = function(x) v + mu^2 - ((mu^2*(gamma(1+2/x)))/(gamma(1+1/x)^2))

registerDoParallel(cores=35)
res.sim.bull.pois <- list()
res.sim.bull.pois <- foreach(i=1:B.gam.pois) %dopar%{
    mu = mean(mat[i,]); v = var(mat[i,])
    k_w = uniroot(fun.k,c(0.1,100))
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

name<-"Simulation_logpois_bullpois.RData"
setwd(out)
save(res.sim.bull.pois,file = name)

############
## DATA = WB

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

load("data_WB.RData")
source("distfunctions.R")
B.gam.pois = 1000

fun.k = function(x) v + mu^2 - ((mu^2*(gamma(1+2/x)))/(gamma(1+1/x)^2))

registerDoParallel(cores=35)
res.sim.bull.pois <- list()
res.sim.bull.pois <- foreach(i=1:B.gam.pois) %dopar%{
    mu = mean(mat[i,]); v = var(mat[i,])
    k_w = uniroot(fun.k,c(0.1,100))
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

name<-"Simulation_bullpois_bullpois.RData"
setwd(out)
save(res.sim.bull.pois,file = name)

############
## DATA = GG

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

load("data_GG.RData")
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

name<-"Simulation_GGpois_bullpois.RData"
setwd(out)
save(res.sim.bull.pois,file = name)
