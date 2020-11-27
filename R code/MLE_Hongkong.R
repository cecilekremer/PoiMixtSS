
######################################################################################################
### DATA EXAMPLE 1                                                                                 ###
### Adam et al. "Clustering and superspreading potential of SARS-C-oV-2 in-fections in Hong Kong". ###
### In:Nature Medicine(2020).doi:https://doi.org/10.1038/s41591-020-1092-0.                        ###
######################################################################################################


library(magrittr)
library(tidyverse)
library(fitdistrplus)
source("distfunctions.R")
library(flexsurv)
library(mvtnorm)
library(tmvtnorm)



## Read in data
transmission_pairs <- read.table(file = "transmission_pairs.txt",header=T,sep=",")

####OFFSPRING DISTRIBUTION ANALYSIS
#count number of offspring per individual infector
offspring <- transmission_pairs %>%
  dplyr::select(infector.case) %>%
  group_by(infector.case) %>%
  count() %>%
  arrange(desc(n))

#count number of terminal infectees including sporadic local cases
infectee <- transmission_pairs %>%
  dplyr::select(infector.case, infectee.case) %>%
  gather() %>%
  filter(key == 'infectee.case')

infector <-  transmission_pairs %>%
  dplyr::select(infector.case, infectee.case) %>%
  gather() %>%
  filter(key == 'infector.case')

duplicate <- infector %>%
  left_join(., infectee, by = 'value') %>%
  filter(key.y != 'NA') %>%
  dplyr::select(value) %>%
  distinct()

nterminal_infectees <- infectee %>% 
  dplyr::select(value) %>%
  filter(!value %in% duplicate$value) %>%
  transmute(case.no = as.numeric(value)) %>%
  nrow() + 46 #46 Sporadic Local cases without links additional transmission

#create vector of complete offspring distribution with terminal cases having zero secondary cases
complete_offspringd <- enframe(c(offspring$n, rep(0,nterminal_infectees)))

## Empirical offspring distributions
nsec = complete_offspringd$value
hist(nsec)
summary(nsec)
# observed mean & variance
obs.var = var(nsec)
obs.mu = mean(nsec)

#######################################
#### ESTIMATE PARAMETERS USING MLE ####

## NB
log_lik = function(theta){
  m = theta[1]
  k = theta[2]
  logL = sum(log(dnbinom(nsec,size=k,mu=m)))
  return(-logL)
}
opt.NB = optim(c(1,1),log_lik,hessian=T)
AIC.NB = 2*2 + 2*opt.NB$value
AIC.NB
mu = opt.NB$par[1]; k.NB = opt.NB$par[2]
mu
var = mu + mu^2/k
sqrt(var)
mean.NB=mu; var.NB=var

# 95% CI
covmat = solve(opt.NB$hessian)
sample = rmvnorm(100000, mean=c(mu,k), sigma=covmat)
sd.est = c()
for(i in 1:dim(sample)[1]){
  sd.est[i] = sqrt(sample[i,1] + sample[i,1]^2/sample[i,2])
}
quantile(sd.est,c(0.025,0.975))


## POLN
log_lik = function(theta){
  s = theta[1]
  m = theta[2]
  logL = sum(log(dlogpois(nsec,s,m)))
  return(-logL)
}
opt.LN = optim(c(1,-1),log_lik,hessian=T)
mlog = opt.LN$par[2]; slog=opt.LN$par[1]
mean.logpois = as.numeric(exp(mlog + 1/2 * slog^2))
var.logpois = as.numeric(exp(mlog + 1/2 * slog^2) + (exp(2*mlog + slog^2)*(exp(slog^2)-1)))
sqrt(var.logpois)
AIC.LN = 2*2 + 2*opt.LN$value

# 95% CI
covmat = solve(opt.LN$hessian)
sample = rmvnorm(100000, mean=c(slog,mlog), sigma=covmat)
R.est = c(); sd.est = c()
for(i in 1:dim(sample)[1]){
  R.est[i] = exp(sample[i,2] + 0.5*sample[i,1]^2)
  sd.est[i] = sqrt(exp(sample[i,2] + 0.5*sample[i,1]^2) + (exp(2*sample[i,2] + sample[i,1]^2)*(exp(sample[i,1]^2)-1)))
}
quantile(R.est,c(0.025,0.975))
quantile(sd.est,c(0.025,0.975))

## POWB
log_lik = function(theta){
  p = theta[1]
  l = theta[2]
  logL = sum(log(dbullpois(nsec,p,l)))
  return(-logL)
}
fun.k = function(x) obs.var + obs.mu^2 - ((obs.mu^2*(gamma(1+2/x)))/(gamma(1+1/x)^2))
k_w = uniroot(fun.k,c(0.1,100))$root
l_w = obs.mu/(gamma(1+1/k_w))
opt.WB = optim(c(k_w,l_w),log_lik,hessian=T)
k = opt.WB$par[1]; l = opt.WB$par[2]; k; l
mean.bullpois = as.numeric(l*gamma(1+1/k))
var.bullpois = as.numeric(l*gamma(1+1/k) + (l^2*(gamma(1+2/k)-(gamma(1+1/k))^2)))
sqrt(var.bullpois)
AIC.WB = 2*2 + 2*opt.WB$value

# 95% CI
covmat = solve(opt.WB$hessian)
sample = rmvnorm(100000, mean=c(k,l), sigma=covmat)
R.est = c(); sd.est = c()
for(i in 1:dim(sample)[1]){
  R.est[i] = sample[i,2]*gamma(1+1/sample[i,1])
  sd.est[i] = sqrt(sample[i,2]*gamma(1+1/sample[i,1]) + (sample[i,2]^2*(gamma(1+2/sample[i,1])-(gamma(1+1/sample[i,1])^2))))
}
quantile(R.est,c(0.025,0.975))
quantile(sd.est,c(0.025,0.975))

## POGG
log_lik = function(theta){
  a = theta[1]
  p = theta[2]
  d = theta[3]
  logL = sum(log(dgengampois(nsec,a,p,d)))
  return(-logL)
}
opt.GG = optim(c(1,1,1),log_lik,hessian=T)
a.est = opt.GG$par[1]
p.est = opt.GG$par[2]
d.est = opt.GG$par[3]
mean.gengam = as.numeric(a.est*gamma((d.est+1)/p.est)/gamma(d.est/p.est))
var.gengam = as.numeric(a.est*gamma((d.est+1)/p.est)/gamma(d.est/p.est) + (a.est^2*(gamma((d.est+2)/p.est)/gamma(d.est/p.est)-
                                   (gamma((d.est+1)/p.est)/gamma(d.est/p.est))^2)))
sqrt(var.gengam)
AIC.GG = 2*3 + 2*opt.GG$value

# 95% CI
covmat = solve(opt.GG$hessian)
sample = rtmvnorm(100000,mean=c(a.est,p.est,d.est), sigma=covmat, lower=rep(0,3))
R.est = c(); sd.est = c()
for(i in 1:dim(sample)[1]){
  a = sample[i,1]; p = sample[i,2]; d = sample[i,3]
  R.est[i] = a*gamma((d+1)/p)/gamma(d/p)
  sd.est[i] = sqrt(a*gamma((d+1)/p)/gamma(d/p) + (a^2*(gamma((d+2)/p)/gamma(d/p)-(gamma((d+1)/p)/gamma(d/p))^2)))
}
quantile(R.est,c(0.025,0.975))
quantile(sd.est,c(0.025,0.975))

## Plot to observed data
n.sec = nsec
mean.NB = mu; var.NB = var; k.NB = mean.NB^2 / (var.NB-mean.NB)
par(mfrow=c(2,2))
plot(unique(sort(n.sec)),as.numeric(table(n.sec))/length(n.sec), xlim=c(0,max(n.sec)), pch=17, xlab="Number of secondary cases", main="Negative binomial", ylab="Density")
lines(0:max(n.sec), dnbinom(0:max(n.sec), mu=mean.NB, size=k.NB), col="blue", type="b",lty=3)
plot(unique(sort(n.sec)),as.numeric(table(n.sec))/length(n.sec), xlim=c(0,max(n.sec)), pch=17, xlab="Number of secondary cases", main="Poisson-lognormal", ylab="Density")
lines(0:max(n.sec), dlogpois(0:max(n.sec),sigma=slog,mu=mlog), col="purple", type="b",lty=3)
plot(unique(sort(n.sec)),as.numeric(table(n.sec))/length(n.sec), xlim=c(0,max(n.sec)), pch=17, xlab="Number of secondary cases", main="Poisson-Weibull", ylab="Density")
lines(0:max(n.sec), dbullpois(0:max(n.sec), p=k, l=l), col="orange", type="b",lty=3)
plot(unique(sort(n.sec)),as.numeric(table(n.sec))/length(n.sec), xlim=c(0,max(n.sec)), pch=17, xlab="Number of secondary cases", main="Poisson-GG", ylab="Density")
lines(0:max(n.sec), dgengampois(0:max(n.sec),a=a.est,p=p.est,d=d.est), col="red", type="b",lty=3)

########################
### PROP RESPONSIBLE ###

source("funs_cdf.R")
source("distfunctions.R")
library(flexsurv)

### Proportion of cases responsible for 80% of transmission

# Based on Endo et al
# propresponsible(mean.NB,k.NB,0.8)

# Based on Lloyd-Smith et al
prop80.NB = function(s,R){
  fun = function(x) prop.resp.NB(x,s,R) - 0.8
  p = uniroot(fun,c(0.0001,1))$root
  return(p)
}
prop80.NB(sqrt(var.NB),mean.NB)

prop80.LN = function(s,R){
  fun = function(x) prop.resp.LN(x,s,R) - 0.8
  p = uniroot(fun,c(0.0001,1))$root
  return(p)
}
prop80.LN(sqrt(var.logpois),mean.logpois)

prop80.WB = function(s,R){
  fun = function(x) prop.resp.WB(x,s,R) - 0.8
  p = uniroot(fun,c(0.0001,1))$root
  return(p)
}
prop80.WB(sqrt(var.bullpois),mean.bullpois)

prop80.GG = function(R,a,p,d){
  fun = function(x) prop.resp.GG(x,R,a,p,d) - 0.8
  p = uniroot(fun,c(0.0001,1))$root
  return(p)
}
prop80.GG(mean.gengam,a.est,p.est,d.est)

## Figure 2a
par(mfrow=c(1,1))
a = seq(0.0001,1,0.001)
prop.NB <-c()
s=sqrt(var.NB); R=mean.NB
for(i in 1:length(a)){prop.NB[i] = prop.resp.NB(a[i],s,R)}
# jpeg("propAdam.jpeg", width = 18, height = 15, units = 'cm', res = 1200)
plot(a,prop.NB,type="l",col="blue", xlab="Proportion of infectious cases",ylab="Expected proportion of transmission", main="(a)")
prop.LN <-c()
s=sqrt(var.logpois); R=mean.logpois
for(i in 1:length(a)){prop.LN[i] = prop.resp.LN(a[i],s,R)}                                                              
lines(a,prop.LN,col="purple")
prop.WB <-c()
s=sqrt(var.bullpois); R=mean.bullpois
for(i in 1:length(a)){prop.WB[i] = prop.resp.WB(a[i],s,R)}                                                              
lines(a,prop.WB,col="orange")
prop.GG <-c()
s=sqrt(var.gengam); R=mean.gengam
for(i in 1:length(a)){prop.GG[i] = prop.resp.GG(a[i],R,a.est,p.est,d.est)}
lines(a,prop.GG,col="red")
legend(x=0.75,y=0.2,legend=c("Gamma","Lognoral","Weibull","GenGamma"),lty=c(1,1,1,1),col=c("blue","purple","orange","red"),cex=0.8)
abline(h=0.8,lty=2)
# p80%
abline(v=a[289],lty=3,col="blue")
abline(v=a[333],lty=3,col="purple")
abline(v=a[295],lty=3,col="orange")
abline(v=a[305],lty=3,col="red")
# dev.off()

## Proportions based on Poisson mixture
mean.NB=mu; var.NB=var
prop.resp.POGA(mean.NB,sqrt(var.NB),0.8)
prop.resp.POLN(mean.logpois,sqrt(var.logpois),0.8)
prop.resp.POWB(mean.bullpois,sqrt(var.bullpois),0.8)
prop.resp.POGG(mean.gengam,a.est,p.est,d.est,0.8)

prop=seq(0,1,0.01)
p.NB = c()
for(i in 1:length(prop)) p.NB[i] = prop.resp.POGA(mean.NB,sqrt(var.NB),prop[i])
p.LN = c()
for(i in 1:length(prop)) p.LN[i] = prop.resp.POLN(mean.logpois,sqrt(var.logpois),prop[i])
p.WB = c()
for(i in 1:length(prop)) p.WB[i] = prop.resp.POWB(mean.bullpois,sqrt(var.bullpois),prop[i])
p.GG = c()
for(i in 1:length(prop)) p.GG[i] = prop.resp.POGG(mean.gengam,a.est,p.est,d.est,prop[i])

# uncertainty
prop.NB = prop.resp.POGA.discrete(mean.NB,sqrt(var.NB))
prop.LN = prop.resp.POLN.discrete(mean.logpois,sqrt(var.logpois))
prop.WB = prop.resp.POWB.discrete(mean.bullpois,sqrt(var.bullpois))
prop.GG = prop.resp.POGG.discrete(mean.gengam,a.est,p.est,d.est)
max.prop = c(max(prop.NB$x2), max(prop.LN$x2), max(prop.WB$x2), max(prop.GG$x2))
max.p=max(max.prop)

## Figure 2b
# jpeg("propAdam3.jpeg", width = 18, height = 15, units = 'cm', res = 1200)
plot(prop.NB$x1,prop.NB$Prop,type="n",pch=15,col="dimgray",main="(b)",xlab="Proportion of infectious cases",ylab="Realized proportion of transmission",ylim=c(0,1),xlim=c(0,max.p+0.05))
polygon(c(prop.NB$x1,rev(prop.NB$x2)),c(prop.NB$Prop,rev(prop.NB$Prop)),col=rgb(0.2,0.6,0.8,0.2),border=NA)
lines(p.NB,prop,col="blue")
polygon(c(prop.LN$x1,rev(prop.LN$x2)),c(prop.LN$Prop,rev(prop.LN$Prop)),col=rgb(0.6,0.4,0.8,0.2),border=NA)
lines(p.LN,prop,col="purple")
polygon(c(prop.WB$x1,rev(prop.WB$x2)),c(prop.WB$Prop,rev(prop.WB$Prop)),col=rgb(1,0.8,0.4,0.3),border=NA)
lines(p.WB,prop,col="orange")
polygon(c(prop.GG$x1,rev(prop.GG$x2)),c(prop.GG$Prop,rev(prop.GG$Prop)),col=rgb(1,0,0,0.15),border=NA)
lines(p.GG,prop,col="red")
legend(x=0.25,y=0.3,c("NB","POLN","POWB","POGG"),lty=c(1,1,1,1),col=c("blue","purple","orange","red"))
# dev.off()


#######################
### GOODNESS-OF-FIT ###

# p0
sum(n.sec==0)/length(n.sec)
dnbinom(0,size=k.NB,mu=mean.NB)
dbullpois(0,k,l)
dlogpois(0,slog,mlog)
dgengampois(0,a.est,p.est,d.est)

## GOF plot (QQ)
summary(n.sec)
observed = c()
for(i in unique(n.sec)){
  observed[i+1] = sum(n.sec==i)/length(n.sec)
}
observed = observed[!is.na(observed)]

predicted.NB = c()
for(i in unique(n.sec)){
  predicted.NB[i+1] = dnbinom(i,size=k.NB,mu=mean.NB)
}
predicted.NB = predicted.NB[!is.na(predicted.NB)]
par(mfrow=c(2,2))
plot(observed,predicted.NB, xlab="Observed", ylab="Predicted", main="Negative binomial"); abline(coef=c(0,1),lty=3)
text(x=0.1,y=0.6,labels=c(paste0("r = ",round(cor(observed,predicted.NB),5))),cex=0.7)

predicted.LN = c()
for(i in unique(n.sec)){
  predicted.LN[i+1] = dlogpois(i,slog,mlog)
}
predicted.LN = predicted.LN[!is.na(predicted.LN)]
sum(predicted.LN)
plot(observed,predicted.LN,xlab="Observed", ylab="Predicted", main="Poisson-lognormal"); abline(coef=c(0,1),lty=3)
text(x=0.1,y=0.6,labels=c(paste0("r = ",round(cor(observed,predicted.LN),5))),cex=0.7)

predicted.WB = c()
for(i in unique(n.sec)){
  predicted.WB[i+1] = dbullpois(i,k,l)
}
predicted.WB = predicted.WB[!is.na(predicted.WB)]
sum(predicted.WB)
plot(observed,predicted.WB, xlab="Observed", ylab="Predicted", main="Poisson-Weibull"); abline(coef=c(0,1),lty=3)
text(x=0.1,y=0.6,labels=c(paste0("r = ",round(cor(observed,predicted.WB),5))),cex=0.7)

predicted.GG = c()
for(i in unique(n.sec)){
  predicted.GG[i+1] = dgengampois(i,a.est,p.est,d.est)
}
predicted.GG = predicted.GG[!is.na(predicted.GG)]
sum(predicted.GG)
plot(observed,predicted.GG,  xlab="Observed", ylab="Predicted", main="Poisson-GG"); abline(coef=c(0,1),lty=3)
text(x=0.1,y=0.6,labels=c(paste0("r = ",round(cor(observed,predicted.GG),5))),cex=0.7)


####################################
### DISCRETE PARETO DISTRIBUTION ###

library(tolerance)

log_lik = function(theta){
  shape = theta[1]
  logL = sum(log(ddpareto(nsec,shape)))
  return(-logL)
}
opt.GPD = optim(c(0.5),log_lik, method="Brent", lower=0, upper=1,hessian=T)

AIC.GPD = 2 + 2*opt.GPD$value
alpha = opt.GPD$par[1]
sd.alpha = solve(opt.GPD$hessian)

n.sec = nsec
# jpeg("fitPareto.jpeg", width = 22, height = 12, units = 'cm', res = 1200)
par(mfrow=c(1,2))
plot(unique(sort(n.sec)),as.numeric(table(n.sec))/length(n.sec), xlim=c(0,max(n.sec)), ylim=c(0,0.75), pch=17, 
     xlab="Number of secondary cases", ylab="Density")
lines(0:max(n.sec), ddpareto(0:max(n.sec), alpha),type="b",lty=3)

observed = c()
for(i in unique(n.sec)){
  observed[i+1] = sum(n.sec==i)/length(n.sec)
}
observed = observed[!is.na(observed)]

predicted.PR = c()
for(i in unique(n.sec)){
  predicted.PR[i+1] = ddpareto(i,alpha)
}
predicted.PR = predicted.PR[!is.na(predicted.PR)]
sum(predicted.PR)
plot(observed,predicted.PR, xlab="Observed", ylab="Predicted"); abline(coef=c(0,1),lty=3)
text(x=0.1,y=0.6,labels=c(paste0("r = ",round(cor(observed,predicted.PR),5))),cex=0.7)
# dev.off()





