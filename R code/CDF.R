
####################################################################
### Code for Figure 1 (based on simulation setting) & Figure A.1 ###
####################################################################

source("funs_cdf.R")
source("distfunctions.R")
library(flexsurv)

###################################################################
### Calculate prop of transmission due to 20% most infectious cases 

### FIGURE 1a 
R = 0.8; m = R
s = seq(1,3,0.1)
t.20.NB = c(); k = c()
# Gamma
for(i in 1:length(s)){
  v = (s[i]^2) - m
  shape = m^2/v
  scale = v/m
  fun.gamma = function(x) 1 - pgamma(x,shape=shape,scale=scale) - 0.2#prop responsible
  x.a = uniroot(fun.gamma,c(-100,100))$root
  fun.int = function(u) u * dgamma(u,shape=shape,scale=scale) 
  F.trans = function(x) 1/R * integrate(fun.int,0,x)$value
  t.20.NB[i] = 1 - F.trans(x.a) # transmission due to 20% of most infectious cases
}
# jpeg("Fig1a.jpeg", width = 18, height = 15, units = 'cm', res = 1200)
plot(s,t.20.NB,type="l",xlab=expression(paste("Standard deviation (",sigma,")")),ylab="Expected proportion of transmission due to the 20% most infectious cases",
     col="blue")
# lognormal
t.20.LN = c(); k = c()
for(i in 1:length(s)){
  v = (s[i]^2) - m
  mlog = log(m^2/sqrt(v + m^2))
  slog = sqrt(log(1+(v/m^2)))
  fun.logn = function(x) 1 - plnorm(x,mlog,slog) - 0.2#prop responsible
  x.a = uniroot(fun.logn,c(-100,100))$root
  fun.int = function(u) u * dlnorm(u,mlog,slog) 
  F.trans = function(x) 1/R * integrate(fun.int,0,x)$value
  t.20.LN[i] = 1 - F.trans(x.a) # transmission due to 20% of most infectious cases
}
lines(s,t.20.LN,col="purple")
# Weibull
t.20.WB = c(); k = c()
for(i in 1:length(s)){
  v = (s[i]^2) - m
  fun.k = function(x) v + m^2 - ((m^2*(gamma(1+2/x)))/(gamma(1+1/x)^2))
  k_w = uniroot(fun.k,c(0.1,100))$root
  l_w = m/(gamma(1+1/k_w))
  fun.weib = function(x) 1 - pweibull(x,shape=k_w, scale=l_w) - 0.2#prop responsible
  x.a = uniroot(fun.weib,c(-100,100))$root
  fun.int = function(u) u * dweibull(u,shape=k_w, scale=l_w) 
  F.trans = function(x) 1/R * integrate(fun.int,0,x)$value
  t.20.WB[i] = 1 - F.trans(x.a) # transmission due to 20% of most infectious cases
}
lines(s,t.20.WB,col="orange")
# Generalized Gamma
t.20.GG = c(); k = c()
a.g = c(0.3,0.69,1.355)
p.g = c(1.1,0.71,0.5)
d.g = c(3.3,0.68,0.21)
s = c(1,1.5,3)
for(i in 1:length(s)){
  v = (s[i]^2) - m
  a = a.g[i]
  p = p.g[i]  
  d = d.g[i]
  fun.gengam = function(x) 1 - pgengamma.orig(x,shape=p,scale=a,k=d/p) - 0.2#prop responsible
  x.a = uniroot(fun.gengam,c(-100,100))$root
  fun.int = function(u) u * dgengamma.orig(u,shape=p,scale=a,k=d/p)
  F.trans = function(x) 1/R * integrate(fun.int,0,x)$value
  t.20.GG[i] = 1 - F.trans(x.a) # transmission due to 20% of most infectious cases
}
points(s,t.20.GG,col="red",pch=18)
legend(x=1,y=0.95,legend=c("Gamma","Lognoral","Weibull","GenGamma"), lty=c(1,1,1,NA), pch=c(NA,NA,NA,18), col=c("blue","purple","orange","red"))
# dev.off()

##################
### FIGURE 1c 
R = 0.8; m = R
s = seq(1,3,0.1)
# NB
t.20.NB = c()
for(i in 1:length(s)){
  prop.20.fun = function(x) prop.resp.POGA(m,s[i],x) - 0.2
  prop.20 = uniroot(prop.20.fun,c(0,1))$root
  t.20.NB[i] = prop.20
}
# jpeg("Fig1c.jpeg", width = 18, height = 15, units = 'cm', res = 1200)
plot(s,t.20.NB,type="l",xlab=expression(paste("Standard deviation (",sigma,")")),
     ylab="Realized proportion of transmission due to the 20% most infectious cases",
     col="blue")
# POLN
t.20.LN = c()
for(i in 1:length(s)){
  prop.20.fun = function(x) prop.resp.POLN(m,s[i],x) - 0.2
  prop.20 = uniroot(prop.20.fun,c(0,1))$root
  t.20.LN[i] = prop.20
}
lines(s,t.20.LN,col="purple")
# POWB
t.20.WB = c()
for(i in 1:length(s)){
  prop.20.fun = function(x) prop.resp.POWB(m,s[i],x) - 0.2
  prop.20 = uniroot(prop.20.fun,c(0,1))$root
  t.20.WB[i] = prop.20
}
lines(s,t.20.WB,col="orange")
# POGG
t.20.GG = c()
a.g = c(0.3,0.69,1.355)
p.g = c(1.1,0.71,0.5)
d.g = c(3.3,0.68,0.21)
s = c(1,1.5,3)
for(i in 1:length(s)){
  a = a.g[i]
  p = p.g[i]
  d = d.g[i]
  prop.20.fun = function(x) prop.resp.POGG(m,a,p,d,x) - 0.2
  prop.20 = uniroot(prop.20.fun,c(0,1))$root
  t.20.GG[i] = prop.20
}
points(s,t.20.GG,col="red",pch=18)
legend(x=1,y=1,legend=c("NB","POLN","POWB","POGG"), lty=c(1,1,1,NA), pch=c(NA,NA,NA,18), col=c("blue","purple","orange","red"))
# dev.off()

### Account for uncertainty
R = 0.8; m = R
s = c(1,1.5,3)
# NB
t.NB = list()
for(i in 1:length(s)){
  t.NB[[i]] = prop.resp.POGA.discrete(R,s[i])
}
jpeg("Fig1c.jpeg", width = 18, height = 15, units = 'cm', res = 1200)
plot(seq(1,3,0.1),t.20.NB,type="l",xlab=expression(paste("Standard deviation (",sigma,")")),
     ylab="Realized proportion of transmission due to the 20% most infectious cases",
     col="blue")
x1 = t.NB[[1]]$x1[59:97]; x2 = t.NB[[1]]$x2[59:97]; Prop = t.NB[[1]]$Prop[59:97]
Prop = seq(Prop[1],Prop[length(Prop)],0.001)
points(rep(1,length(Prop)), Prop,col=rgb(0.2,0.6,0.8,0.2))
x1 = t.NB[[2]]$x1[76:97]; x2 = t.NB[[2]]$x2[76:97]; Prop = t.NB[[2]]$Prop[76:97]
Prop = seq(Prop[1],Prop[length(Prop)],0.001)
points(rep(1.5,length(Prop)), Prop,col=rgb(0.2,0.6,0.8,0.2))
t.NB[[3]] # nothing >= 20%
# POLN
lines(seq(1,3,0.1),t.20.LN,col="purple")
s = c(1,1.5,3)
t.LN = list()
for(i in 1:length(s)){
  t.LN[[i]] = prop.resp.POLN.discrete(R,s[i])
}
x1 = t.LN[[1]]$x1[59:97]; x2 = t.LN[[1]]$x2[59:97]; Prop = t.LN[[1]]$Prop[59:97]
Prop = seq(Prop[1],Prop[length(Prop)],0.001)
points(rep(1.02,length(Prop)), Prop,col=rgb(0.6,0.4,0.8,0.2))
x1 = t.LN[[2]]$x1[69:97]; x2 = t.LN[[2]]$x2[69:97]; Prop = t.LN[[2]]$Prop[69:97]
Prop = seq(Prop[1],Prop[length(Prop)],0.001)
points(rep(1.52,length(Prop)), Prop,col=rgb(0.6,0.4,0.8,0.2))
x1 = t.LN[[3]]$x1[78:97]; x2 = t.LN[[3]]$x2[78:97]; Prop = t.LN[[3]]$Prop[78:97]
Prop = seq(Prop[1],Prop[length(Prop)],0.001)
points(rep(3,length(Prop)), Prop,col=rgb(0.6,0.4,0.8,0.2))
# POWB
lines(seq(1,3,0.1),t.20.WB,col="orange")
s = c(1,1.5,3)
t.WB = list()
for(i in 1:length(s)){
  t.WB[[i]] = prop.resp.POWB.discrete(R,s[i])
}
x1 = t.WB[[1]]$x1[60:97]; x2 = t.WB[[1]]$x2[60:97]; Prop = t.WB[[1]]$Prop[60:97]
Prop = seq(Prop[1],Prop[length(Prop)],0.001)
points(rep(1.04,length(Prop)), Prop,col=rgb(1,0.8,0.4,0.3))
x1 = t.WB[[2]]$x1[74:97]; x2 = t.WB[[2]]$x2[74:97]; Prop = t.WB[[2]]$Prop[74:97]
Prop = seq(Prop[1],Prop[length(Prop)],0.001)
points(rep(1.54,length(Prop)), Prop,col=rgb(1,0.8,0.4,0.3))
x1 = t.WB[[3]]$x1[85:97]; x2 = t.WB[[3]]$x2[85:97]; Prop = t.WB[[3]]$Prop[85:97]
Prop = seq(Prop[1],Prop[length(Prop)],0.001)
points(rep(3.02,length(Prop)), Prop,col=rgb(1,0.8,0.4,0.3))
# POGG
a.g = c(0.3,0.69,1.355)
p.g = c(1.1,0.71,0.5)
d.g = c(3.3,0.68,0.21)
s = c(1,1.5,3)
points(s,t.20.GG,col="red",pch=15)
t.GG = list()
for(i in 1:length(s)){
  t.GG[[i]] = prop.resp.POGG.discrete(R,a.g[i],p.g[i],d.g[i])
}
x1 = t.GG[[1]]$x1[59:97]; x2 = t.GG[[1]]$x2[59:97]; Prop = t.GG[[1]]$Prop[59:97]
Prop = seq(Prop[1],Prop[length(Prop)],0.001)
points(rep(1.06,length(Prop)), Prop,col=rgb(1,0,0,0.15))
x1 = t.GG[[2]]$x1[74:97]; x2 = t.GG[[2]]$x2[74:97]; Prop = t.GG[[2]]$Prop[74:97]
Prop = seq(Prop[1],Prop[length(Prop)],0.001)
points(rep(1.56,length(Prop)), Prop,col=rgb(1,0,0,0.15))
x1 = t.GG[[3]]$x1[87:97]; x2 = t.GG[[3]]$x2[87:97]; Prop = t.GG[[3]]$Prop[87:97]
Prop = seq(Prop[1],Prop[length(Prop)],0.001)
points(rep(3.04,length(Prop)), Prop,col=rgb(1,0,0,0.15))
legend(x=2.5,y=0.75,legend=c("NB","POLN","POWB","POGG"), lty=c(1,1,1,NA), pch=c(NA,NA,NA,15), col=c("blue","purple","orange","red"))
dev.off()

########################################################
### Plot prop of cases vs prop of all transmission (cdf)

### FIGURE 1b
# parameters GG
a.g = c(0.3,0.69,1.355)
p.g = c(1.1,0.71,0.5)
d.g = c(3.3,0.68,0.21)
# sigma = 3
s = 3
R = m = 0.8
v = (s^2) - m
k = m^2 / (s^2-m)
v = m + m^2/k
a = seq(0.0001,1,0.01)
prop.NB<-c()
for(i in 1:length(a)){prop.NB[i] = prop.resp.NB(a[i],s,m)}
jpeg("Fig1b.jpeg", width = 18, height = 15, units = 'cm', res = 1200)
plot(a,prop.NB,xlim=c(0,1),type="l",ylim=c(0,1),xlab="Proportion of infectious cases",ylab="Expected proportion of transmission",col="blue")
prop.LN<-c()
for(i in 1:length(a)){prop.LN[i] = prop.resp.LN(a[i],s,m)}
lines(a,prop.LN,col="purple")
prop.WB<-c()
for(i in 1:length(a)){prop.WB[i] = prop.resp.WB(a[i],s,m)}
lines(a,prop.WB,col="orange")
prop.GG<-c()
for(i in 1:length(a)){prop.GG[i] = prop.resp.GG(a[i],m,a.g[3],p.g[3],d.g[3])}
lines(a,prop.GG,col="red")
# sigma = 1.5
s = 1.5
prop.NB<-c()
for(i in 1:length(a)){prop.NB[i] = prop.resp.NB(a[i],s,m)}
lines(a,prop.NB,col="blue",lty=2)
prop.LN<-c()
for(i in 1:length(a)){prop.LN[i] = prop.resp.LN(a[i],s,m)}
lines(a,prop.LN,col="purple",lty=2)
prop.WB<-c()
for(i in 1:length(a)){prop.WB[i] = prop.resp.WB(a[i],s,m)}
lines(a,prop.WB,col="orange",lty=2)
prop.GG<-c()
for(i in 2:length(a)-1){prop.GG[i] = prop.resp.GG(a[i],m,a.g[2],p.g[2],d.g[2])}
lines(a[1:99],prop.GG,col="red",lty=2)
# sigma = 1
s = 1
prop.NB<-c()
for(i in 1:length(a)){prop.NB[i] = prop.resp.NB(a[i],s,m)}
lines(a,prop.NB,col="blue",lty=3)
prop.LN<-c()
for(i in 1:length(a)){prop.LN[i] = prop.resp.LN(a[i],s,m)}
lines(a,prop.LN,col="purple",lty=3)
prop.WB<-c()
for(i in 1:length(a)){prop.WB[i] = prop.resp.WB(a[i],s,m)}
lines(a,prop.WB,col="orange",lty=3)
prop.GG<-c()
for(i in 1:length(a)){prop.GG[i] = prop.resp.GG(a[i],m,a.g[1],p.g[1],d.g[1])}
lines(a,prop.GG,col="red",lty=3)
legend(x=0.7, y =0.25,legend=c("Gamma","Lognoral","Weibull","GenGamma"),lty=c(1,1,1,1),col=c("blue","purple","orange","red"))
dev.off()

##################
### FIGURE 1d 

# parameters GG
a.g = c(0.3,0.69,1.355)
p.g = c(1.1,0.71,0.5)
d.g = c(3.3,0.68,0.21)
# sigma = 3
s = 3
R = m = 0.8
v = (s^2) - m
# k = m^2/(v-m)
# prop = seq(0.01,1,0.02)
prop = seq(0,1,0.01)
prop.NB<-c()
for(i in 1:length(prop)){prop.NB[i] = prop.resp.POGA(m,s,prop[i])}
# prop.NB = c(prop.NB,1)
jpeg("Fig1d.jpeg", width = 18, height = 15, units = 'cm', res = 1200)
plot(prop.NB,prop,type="l",ylim=c(0,1),xlim=c(0,0.55),xlab="Proportion of infectious cases",ylab="Realized proportion of transmission",col="blue")
prop.LN<-c()
for(i in 3:length(prop)){prop.LN[i] = prop.resp.POLN(m,s,prop[i])}
lines(prop.LN,prop,col="purple")
prop.WB<-c()
for(i in 3:length(prop)){prop.WB[i] = prop.resp.POWB(m,s,prop[i])}
lines(prop.WB,prop,col="orange")
prop.GG<-c()
for(i in 3:length(prop)){prop.GG[i] = prop.resp.POGG(m,a.g[3],p.g[3],d.g[3],prop[i])}
lines(prop.GG,prop,col="red")
# sigma = 1.5
s = 1.5
prop.NB<-c()
for(i in 1:length(prop)){prop.NB[i] = prop.resp.POGA(m,s,prop[i])}
lines(prop.NB,prop,col="blue",lty=2)
prop.LN<-c()
for(i in 1:length(prop)){prop.LN[i] = prop.resp.POLN(m,s,prop[i])}
lines(prop.LN,prop,col="purple",lty=2)
prop.WB<-c()
for(i in 1:length(prop)){prop.WB[i] = prop.resp.POWB(m,s,prop[i])}
lines(prop.WB,prop,col="orange",lty=2)
prop.GG<-c()
for(i in 1:length(prop)){prop.GG[i] = prop.resp.POGG(m,a.g[2],p.g[2],d.g[2],prop[i])}
lines(prop.GG,prop,col="red",lty=2)
# sigma = 1
s = 1
prop.NB<-c()
for(i in 1:length(prop)){prop.NB[i] = prop.resp.POGA(m,s,prop[i])}
lines(prop.NB,prop,col="blue",lty=3)
prop.LN<-c()
for(i in 1:length(prop)){prop.LN[i] = prop.resp.POLN(m,s,prop[i])}
lines(prop.LN,prop,col="purple",lty=3)
prop.WB<-c()
for(i in 1:length(prop)){prop.WB[i] = prop.resp.POWB(m,s,prop[i])}
lines(prop.WB,prop,col="orange",lty=3)
prop.GG<-c()
for(i in 1:length(prop)){prop.GG[i] = prop.resp.POGG(m,a.g[1],p.g[1],d.g[1],prop[i])}
lines(prop.GG,prop,col="red",lty=3)
legend(x=0.35, y =0.25,c("NB","POLN","POWB","POGG"),lty=c(1,1,1,1),col=c("blue","purple","orange","red"))
dev.off()


##################################
## Comparing Endo vs Lloyd-Smith 

jpeg("FigA1.jpeg", width = 24, height = 15, units = 'cm', res = 1200)
par(mfrow=c(1,2))

############################
## Figure A.1(a): varying R

# R = 2.5; k=0.4
LLS = c(0,0.4824,0.6961,0.8237,0.9025,0.9504,0.9779,0.9921,0.9981,0.9998,1)
END = c(0,0.5120,0.7331,0.8647,0.9411,0.9811,rep(1,5))
a = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
plot(a,LLS,col="blue",lty=2,type="l",main="(a)", xlab="Proportion of infectious cases",ylab="Proportion of transmission")
lines(a,END,col="red",lty=2)

R=0.8; k=0.4
v = R + R^2/k
s=sqrt(v)
a=1 # a in 0 to 1 by 0.1
propEndo = function(x) propresponsible(R,k,x) - a
prop = uniroot(propEndo,c(0.0001,1))$root; prop
prop.resp.NB(a,s,R)

LLS = c(0,0.4824,0.6961,0.8237,0.9025,0.9504,0.9779,0.9920,0.9981,0.9998,1)
END = c(0,0.5713,0.8055,0.9305,rep(1,7))
a = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
lines(a,LLS,col="blue",lty=3)
lines(a,END,col="red",lty=3)

R=3.4; k=0.4
v = R + R^2/k
s=sqrt(v)
a=0.9
propEndo = function(x) propresponsible(R,k,x) - a
prop = uniroot(propEndo,c(0.0001,1))$root; prop
prop.resp.NB(a,s,R)

LLS = c(0,0.4824,0.6961,0.8237,0.9025,0.9504,0.9779,0.9920,0.9981,0.9998,1)
END = c(0,0.5040,0.7246,0.8533,0.9289,0.9725,rep(1,5))
a = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
lines(a,LLS,col="blue",lty=5)
lines(a,END,col="red",lty=5)

R=5; k=0.4
v = R + R^2/k
s=sqrt(v)
a=0.7
propEndo = function(x) propresponsible(R,k,x) - a
prop = uniroot(propEndo,c(0.0001,1))$root; prop
prop.resp.NB(a,s,R)

LLS = c(0,0.4824,0.6961,0.8237,0.9025,0.9504,0.9779,0.9920,0.9981,0.9998,1)
END = c(0,0.4975,0.7153,0.8435,0.9211,0.9674,0.9906,rep(1,4))
a = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
lines(a,LLS,col="blue",lty=4)
lines(a,END,col="red",lty=4)

R=0.4; k=0.4
v = R + R^2/k
s=sqrt(v)
a=0.3
propEndo = function(x) propresponsible(R,k,x) - a
prop = uniroot(propEndo,c(0.0001,1))$root; prop
prop.resp.NB(a,s,R)

LLS = c(0,0.4824,0.6961,0.8237,0.9025,0.9504,0.9779,0.9920,0.9981,0.9998,1)
END = c(0,0.6447,0.8947,rep(1,8))
a = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
lines(a,LLS,col="blue",lty=6)
lines(a,END,col="red",lty=6)



legend(x=0.6,y=0.35,legend=c("Lloyd-Smith (R=2.5)","Endo et al (R=2.5)","Lloyd-Smith (R=0.8)","Endo et al (R=0.8)","Lloyd-Smith (R=3.4)",
                              "Endo et al (R=3.4)","Lloyd-Smith (R=5)","Endo et al (R=5)","Lloyd-Smith (R=0.4)","Endo et al (R=0.4)"),
       col=c("blue","red","blue","red","blue","red","blue","red","blue","red"),lty=c(2,2,3,3,5,5,4,4,6,6),cex=0.6)



### Figure A.1(b)
# k = 0.1; R = 2.5
a = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
LLS = c(0,0.8056,0.9510,0.9890,0.9980,0.9997,0.999977,0.9999986,0.9999996,1,1)
END = c(0,0.8194,0.9653,rep(1,8))
plot(a,LLS,col="blue", lty=3,type="l", xlab="Proportion of infectious cases",ylab="Proportion of transmission",main="(b)")
lines(a,END,col="red",lty=3)

R=2.5; k=0.4
v = R + R^2/k
s=sqrt(v)
a=0.1 # a in 0 to 1 by 0.1
propEndo = function(x) propresponsible(R,k,x) - a
prop = uniroot(propEndo,c(0.0001,1))$root; prop
prop.resp.NB(a,s,R)

LLS = c(0,0.4824,0.6961,0.8237,0.9025,0.9504,0.9779,0.9921,0.9981,0.9998,1)
END = c(0,0.5120,0.7331,0.8647,0.9411,0.9811,rep(1,5))
a = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
lines(a,LLS,col="blue",lty=2)
lines(a,END,col="red",lty=2)

R=2.5; k=1.2
v = R + R^2/k
s=sqrt(v)
a=0.5
propEndo = function(x) propresponsible(R,k,x) - a
prop = uniroot(propEndo,c(0.0001,1))$root; prop
prop.resp.NB(a,s,R)

LLS = c(0,0.3075,0.4932,0.6320,0.7397,0.8238,0.8888,0.9375,0.9716,0.9924,1)
END = c(0,0.3522,0.5560,0.7011,0.8111,0.8911,0.9436,0.9836,rep(1,3))
a = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
lines(a,LLS,col="blue",lty=1)
lines(a,END,col="red",lty=1)

R=2.5; k=10^10
v = R + R^2/k
s=sqrt(v)
a=1
propEndo = function(x) propresponsible(R,k,x) - a
prop = uniroot(propEndo,c(0.0001,1))$root; prop
prop.resp.NB(a,s,R)

LLS = c(0,0.0950,0.1933,0.2962,0.4051,0.5033,0.6022,0.6978,0.8038,0.8964,1)
END = c(0,0.2248,0.3883,0.5253,0.6453,0.7478,0.8278,0.9077,0.9528,0.9928,1)
a = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
lines(a,LLS,col="blue",lty=5)
lines(a,END,col="red",lty=5)

legend(x=0.6,y=0.3,legend=c("Lloyd-Smith (k=0.1)","Endo et al (k=0.1)","Lloyd-Smith (k=0.4)","Endo et al (k=0.4)","Lloyd-Smith (k=1.2)","Endo et al (k=1.2)",
                              "Lloyd-Smith (k=1e10)","Endo et al (k=1e10)"),
       col=c("blue","red","blue","red","blue","red","blue","red"),lty=c(3,3,2,2,1,1,5,5),cex=0.6)
dev.off()

###############################################
###############################################
### WHAT IF R>1  ? 
###

### FIGURE 1b --> only changes when k changes
# parameters GG?

R = 3; k = 0.4
m = R
v = m + m^2/k
s = sqrt(v)

a = seq(0.0001,1,0.01)
prop.NB<-c()
for(i in 1:length(a)){prop.NB[i] = prop.resp.NB(a[i],s,m)}
# jpeg("Fig1b.jpeg", width = 18, height = 15, units = 'cm', res = 1200)
plot(a,prop.NB,xlim=c(0,1),type="l",ylim=c(0,1),xlab="Proportion of infectious cases",ylab="Expected proportion of transmission",col="blue")
prop.LN<-c()
for(i in 1:length(a)){prop.LN[i] = prop.resp.LN(a[i],s,m)}
lines(a,prop.LN,col="purple")
prop.WB<-c()
for(i in 1:length(a)){prop.WB[i] = prop.resp.WB(a[i],s,m)}
lines(a,prop.WB,col="orange")
prop.GG<-c()
for(i in 1:length(a)){prop.GG[i] = prop.resp.GG(a[i],m,a.g[3],p.g[3],d.g[3])}
lines(a,prop.GG,col="red")

legend(x=0.7, y =0.25,legend=c("Gamma","Lognoral","Weibull","GenGamma"),lty=c(1,1,1,1),col=c("blue","purple","orange","red"))
# dev.off()

##################
### FIGURE 1d 

R = 1.2; k = 0.08
m = R
v = m + m^2/k
s = sqrt(v)
prop = seq(0,1,0.01)
prop.NB<-c()
for(i in 1:length(prop)){prop.NB[i] = prop.resp.POGA(m,s,prop[i])}
jpeg("FigA2.jpeg", width = 18, height = 15, units = 'cm', res = 1200)
plot(prop.NB,prop,type="l",ylim=c(0,1),xlab="Proportion of infectious cases", xlim=c(0,0.8),ylab="Realized proportion of transmission",col="blue")
prop.LN<-c()
for(i in 3:length(prop)){prop.LN[i] = prop.resp.POLN(m,s,prop[i])}
lines(prop.LN,prop,col="purple")
# prop.WB<-c()
# for(i in 3:length(prop)){prop.WB[i] = prop.resp.POWB(m,s,prop[i])}
# lines(prop.WB,prop,col="orange")

R = 1.2; k = 0.44
m = R
v = m + m^2/k
s = sqrt(v)
prop.NB<-c()
for(i in 1:length(prop)){prop.NB[i] = prop.resp.POGA(m,s,prop[i])}
lines(prop.NB,prop,col="blue",lty=2)
prop.LN<-c()
for(i in 1:length(prop)){prop.LN[i] = prop.resp.POLN(m,s,prop[i])}
lines(prop.LN,prop,col="purple",lty=2)
# prop.WB<-c()
# for(i in 1:length(prop)){prop.WB[i] = prop.resp.POWB(m,s,prop[i])}
# lines(prop.WB,prop,col="orange",lty=2)

R = 1.2; k = 3.2
m = R
v = m + m^2/k
s = sqrt(v)
prop.NB<-c()
for(i in 1:length(prop)){prop.NB[i] = prop.resp.POGA(m,s,prop[i])}
lines(prop.NB,prop,col="blue",lty=3)
prop.LN<-c()
for(i in 1:length(prop)){prop.LN[i] = prop.resp.POLN(m,s,prop[i])}
lines(prop.LN,prop,col="purple",lty=3)
# prop.WB<-c()
# for(i in 1:length(prop)){prop.WB[i] = prop.resp.POWB(m,s,prop[i])}
# lines(prop.WB,prop,col="orange",lty=3)

legend(x=0.6, y =0.25,c("NB","POLN"),lty=c(1,1),col=c("blue","purple"))
dev.off()



