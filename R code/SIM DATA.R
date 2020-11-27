
##########################################
### Simulate data for simulation study ###
##########################################

library(fitdistrplus)
library(pracma)
library(MASS)
library(flexsurv)
library(foreach)
library(doParallel)
source("distfunctions.R")

## INITIALISE PARAMETERS
n = 10000
mu = 0.8
v = (3^2)-0.8 # variance of mixing distribution = variance of mixture - mean
# v = (1.5^2)-0.8
# v = (1^2)-0.8
s = sqrt(v)
k = mu^2/((mu+v)-mu)

# Set parameters 
# gamma
a = mu^2/v # gamma shape = k
b = v/mu # gamma scale
# lognormal
mlog = log(mu^2/sqrt(v + mu^2))
slog = sqrt(log(1+(v/mu^2)))
mean.LN = as.numeric(exp(mlog + 1/2 * slog^2))
var.LN = as.numeric(exp(mlog + 1/2 * slog^2) + (exp(2*mlog + slog^2)*(exp(slog^2)-1)))
# weibull
fun.k = function(x) v + mu^2 - ((mu^2*(gamma(1+2/x)))/(gamma(1+1/x)^2))
k_w = uniroot(fun.k,c(0.1,100))$root
l_w = mu/(gamma(1+1/k_w))
mean.WB = as.numeric(l_w*gamma(1+1/k_w))
var.WB =  as.numeric(l_w*gamma(1+1/k_w) + (l_w^2*(gamma(1+2/k_w)-(gamma(1+1/k_w))^2)))

## SIMULATE DATA POISSON
B.gam.pois=1000
mat = zeros(B.gam.pois,n)
set.seed(251020)

for(sim in 1:B.gam.pois){
  Y.pois = rpois(n,lambda=mu)
  mat[sim,] = Y.pois
}
# check mean and sd
mean(unlist(mat))
sd(unlist(mat))

## SIMULATE DATA NB
B.gam.pois=1000
mat = zeros(B.gam.pois,n)
set.seed(251020)

for(sim in 1:B.gam.pois){
  Y.gam.pois = rnbinom(n,mu=mu,size=k)
  mat[sim,] = Y.gam.pois
}
mean(unlist(mat))
sd(unlist(mat))

## SIMULATE DATA POLN
B.gam.pois = 1000
mat = zeros(B.gam.pois,n)
set.seed(251020)

for(sim in 1:B.gam.pois){
  Y.logn.pois = rpois(n,lambda=rlnorm(n,meanlog=mlog,sdlog=slog))
  mat[sim,] = Y.logn.pois
}
mean(unlist(mat))
sd(unlist(mat))

## SIMULATE DATA POWB
B.gam.pois = 1000
mat = zeros(B.gam.pois,n)
set.seed(251020)

for(sim in 1:B.gam.pois){
  Y.bull.pois = rpois(n,lambda=rweibull(n,shape=k_w,scale=l_w)) 
  mat[sim,] = Y.bull.pois
}
mean(unlist(mat))
sd(unlist(mat))

## SIMULATE DATA POGG
B.gam.pois = 1000
mat = zeros(B.gam.pois,n)
set.seed(251020)

# GG sd = 3
a.g = 1.355
p.g = 0.5 
d.g = 0.21

# GG sd = 1.5 
# a.g = 0.69
# p.g = 0.71 
# d.g = 0.68

# GG sd = 1
# a.g = 0.3
# p.g = 1.1
# d.g = 3.3

as.numeric(a.g*gamma((d.g+1)/p.g)/gamma(d.g/p.g)) #mean
sqrt(as.numeric(a.g*gamma((d.g+1)/p.g)/gamma(d.g/p.g)) + (a.g^2*(gamma((d.g+2)/p.g)/gamma(d.g/p.g) -
    (gamma((d.g+1)/p.g)/gamma(d.g/p.g))^2))) #sd

Q.g = 1/sqrt(d.g/p.g)
mu.g = log(a.g)+(1/p.g)*log(1/Q.g^2)
sigma.g = Q.g/p.g

for(sim in 1:B.gam.pois){
    Y.gam.pois = rpois(n,lambda=flexsurv::rgengamma(n,mu=mu.g,sigma=sigma.g,Q=Q.g)) # check parameter specification !!
    mat[sim,] = Y.gam.pois
}
mean(unlist(mat))
sd(unlist(mat))



