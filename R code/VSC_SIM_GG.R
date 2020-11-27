
##################################################
### Fit POGG to data from NB, POLN, POWB, POGG ###
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

name<-"Simulation_gampois_gengampois.RData"
setwd(out)
save(res.sim.gengam.pois,file = name)

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

registerDoParallel(cores=70)
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

name<-"Simulation_logpois_gengampois.RData"
setwd(out)
save(res.sim.gengam.pois,file = name)

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

registerDoParallel(cores=70)
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

name<-"Simulation_bullpois_gengampois.RData"
setwd(out)
save(res.sim.gengam.pois,file = name)

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

name<-"Simulation_GGpois_gengampois.RData"
setwd(out)
save(res.sim.gengam.pois,file = name)

