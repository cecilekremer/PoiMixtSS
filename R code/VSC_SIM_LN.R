
##################################################
### Fit POLN to data from NB, POLN, POWB, POGG ###
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

name<-"Simulation_gampois_logpois.RData"
setwd(out)
save(res.sim.log.pois,file = name)

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

name<-"Simulation_logpois_logpois.RData"
setwd(out)
save(res.sim.log.pois,file = name)

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
library(parallel)

load("data_WB.RData")
source("distfunctions.R")

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

name<-"Simulation_bullpois_logpois.RData"
setwd(out)
save(res.sim.log.pois,file = name)

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

name<-"Simulation_GGpois_logpois.RData"
setwd(out)
save(res.sim.log.pois,file = name)

