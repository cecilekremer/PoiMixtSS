
################################################
### Fit NB to data from NB, POLN, POWB, POGG ###
################################################


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
res.sim.gam.pois <- list()
res.sim.gam.pois <- foreach(i=1:B.gam.pois) %dopar%{
    fit.gampois = fitdistrplus::fitdist(mat[i,],"nbinom")
    res.sim.gam.pois[[i]] = c(fit.gampois$estimate[2],fit.gampois$estimate[1],fit.gampois$sd[2],fit.gampois$sd[1],fit.gampois$aic,fit.gampois$loglik)
}

name<-"Simulation_gampois_gampois.RData"
setwd(out)
save(res.sim.gam.pois,file = name)

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
res.sim.gam.pois<-list()
res.sim.gam.pois <- foreach(i=1:B.gam.pois) %dopar%{
    fit.gampois = fitdistrplus::fitdist(mat[i,],"nbinom")
    res.sim.gam.pois[[i]] = c(fit.gampois$estimate[2],fit.gampois$estimate[1],fit.gampois$sd[2],fit.gampois$sd[1],fit.gampois$aic,fit.gampois$loglik)
}

name<-"Simulation_logpois_gampois.RData"
setwd(out)
save(res.sim.gam.pois,file = name)

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

registerDoParallel(cores=35)
res.sim.gam.pois<-list()
res.sim.gam.pois <- foreach(i=1:B.gam.pois) %dopar%{
    fit.gampois = fitdistrplus::fitdist(mat[i,],"nbinom")
    res.sim.gam.pois[[i]] = c(fit.gampois$estimate[2],fit.gampois$estimate[1],fit.gampois$sd[2],fit.gampois$sd[1],fit.gampois$aic,fit.gampois$loglik)
}

name<-"Simulation_bullpois_gampois.RData"
setwd(out)
save(res.sim.gam.pois,file = name)


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
res.sim.gam.pois <- list()
res.sim.gam.pois <- foreach(i=1:B.gam.pois) %dopar%{
    fit.gampois = fitdistrplus::fitdist(mat[i,],"nbinom")
    res.sim.gam.pois[[i]] = c(fit.gampois$estimate[2],fit.gampois$estimate[1],fit.gampois$sd[2],fit.gampois$sd[1],fit.gampois$aic,fit.gampois$loglik)
}

name<-"Simulation_GGpois_gampois.RData"
setwd(out)
save(res.sim.gam.pois,file = name)

