
###################################################################################
### Functions for obtaining prop responsible for certain amount of transmission ###
###################################################################################

########################################
#### BASED ON THE DISTRIBUTION OF NU ###

## Calculate proportion of transmission due to a% of most infectious cases

# Gamma
prop.resp.NB = function(a,s,R){
  m = R
  v = (s^2)-m + 0.000001 
  shape = m^2/v; scale = v/m # parameters Gamma
  fun.gamma = function(x) 1 - pgamma(x,shape=shape,scale=scale) # cdf 1-F(x): proportion of cases with R > x
  x.fun = function(x) fun.gamma(x) - a # find x for which 20% of cases has a larger reprod number = most infectious cases
  x.a.NB = uniroot(x.fun,c(-100,100))$root # 
  fun.int.NB = function(u) u * dgamma(u,shape=shape,scale=scale) # u * pdf 
  F.trans.NB = function(x) 1/R * integrate(fun.int.NB,0,x)$value # CDF for disease transmission: expected prop of transmission due to cases with R < x
  t.a.NB = 1 - F.trans.NB(x.a.NB) # transmission due to a% of most infectious cases
  return(t.a.NB) 
}
# Lognormal
prop.resp.LN = function(a,s,R){
  m = R
  v = (s^2)-m + 0.000001 
  mlog = log(m^2/sqrt(v + m^2))
  slog = sqrt(log(1+(v/m^2)))
  fun.logn = function(x) 1 - plnorm(x,mlog,slog) - a 
  x.a.LN = uniroot(fun.logn,c(-100,100))$root
  fun.int.LN = function(u) u * dlnorm(u,mlog,slog) 
  F.trans.LN = function(x) 1/R * integrate(fun.int.LN,0,x)$value
  t.a.LN = 1 - F.trans.LN(x.a.LN) # transmission due to a% of most infectious cases
  return(t.a.LN)
}
# Weibull
prop.resp.WB = function(a,s,R){
  m = R
  v = (s^2)-m
  fun.k = function(x) v + m^2 - ((m^2*(gamma(1+2/x)))/(gamma(1+1/x)^2))
  k_w = uniroot(fun.k,c(0.1,100))$root
  l_w = m/(gamma(1+1/k_w))
  fun.weib = function(x) 1 - pweibull(x,shape=k_w, scale=l_w) - a #prop responsible
  x.a = uniroot(fun.weib,c(-100,100))$root
  fun.int = function(u) u * dweibull(u,shape=k_w, scale=l_w) 
  F.trans = function(x) 1/R * integrate(fun.int,0,x)$value
  t.a.WB = 1 - F.trans(x.a) # transmission due to a% of most infectious cases
  return(t.a.WB)
}
# Generalized Gamma
prop.resp.GG = function(a,R,a.g,p.g,d.g){
  Q.g = 1/sqrt(d.g/p.g)
  mu.g = log(a.g)+(1/p.g)*log(1/Q.g^2)
  sigma.g = Q.g/p.g
  fun.gengam = function(x) 1 - pgengamma(x,mu=mu.g,sigma=sigma.g,Q=Q.g) - a
  x.a = uniroot(fun.gengam,c(-100,100))$root
  fun.int = function(u) u * dgengamma.orig(u,shape=p.g,scale=a.g,k=d.g/p.g)
  F.trans = function(x) 1/R * integrate(fun.int,0,x)$value
  t.a.GG = 1 - F.trans(x.a) # transmission due to 20% of most infectious cases
  return(t.a.GG)
}

#####################################
#### BASED ON THE POISSON MIXTURE ###

## Calculate proportion responsible for X% (prop) of total transmissions

## Poisson-Gamma (=NB)
prop.resp.POGA = function(R,s,prop){ 
  m = R
  v = (s^2)
  shape = m^2/(v-m); rate = m/(v-m) 
  if(prop==0) return(0)
  fun = function(x) floor(x)*dgampois(floor(x),shape,rate)
  F.trans = function(x) 1/R * integrate(fun, lower=0, upper=x, rel.tol = 0.000000000001, subdivisions=10000,stop.on.error = F)$value # cdf of disease transmission: 
  fun2 = function(x) F.trans(x) - (1-prop) 
  x.a.GA = uniroot(fun2,c(0,100))$root  # x such that (1-p)% of transmission is due to cases with R < x
  pmf = function(x) dgampois(floor(x),shape,rate) 
  p.a = function(x) integrate(pmf, lower=0, upper=x, rel.tol = 0.000000000001, subdivisions = 10000)$value # F(X-1) = P(X <= x-1) = P(X < x)
  p.a.GA = 1-p.a(x.a.GA) # 1-F(X-1) = P(X>x-1) = P(X >= x)
  return(p.a.GA)
}
# discrete version (uncertainty)
prop.resp.POGA.discrete = function(R,s){
  m = R
  v = (s^2)
  shape = m^2/(v-m); rate = m/(v-m) 
  fun.int.GA = function(u){
    f=c()
    u = seq(0,u)
    for(i in 1:length(u)){f[i] = u[i] * dgampois(u[i],shape,rate)}
    return(sum(f))
  } 
  F.trans = function(x) 1/R * fun.int.GA(x)
  d<-seq(0.03,1,0.01)
  prop.NB<-data.frame("Prop"=rep(0,length(d)),"x1"=rep(0,length(d)),"x2"=rep(0,length(d)))
  k<-1
  for (i in d){
    fun1 = function(x) F.trans(x) - (1-i)
    temp.inv <- floor(uniroot(fun1,c(0,100))$root)
    if (F.trans(temp.inv)>(1-i)){ 
      x1.temp<-temp.inv-1
      x2.temp<-temp.inv
    }
    if (F.trans(temp.inv)<(1-i)){
      x1.temp<-temp.inv
      x2.temp<-temp.inv+1
    }
    if (F.trans(temp.inv)==(1-i)){
      x1.temp<-temp.inv
      x2.temp<-temp.inv
    }
    prop.NB[k,]<-c(i,(1- pnbinom(x2.temp,size=shape,mu=m)) ,(1- pnbinom(x1.temp,size=shape,mu=m)))
    k<-k+1
  }
  return(prop.NB)
}

## Poisson-LN 
prop.resp.POLN = function(R,s,prop){ 
  m = R
  v = (s^2)-m
  mlog = log(m^2/sqrt(v + m^2))
  slog = sqrt(log(1+(v/m^2)))
  if(prop==0) return(0)
  fun = function(x) floor(x)*dlogpois(floor(x),slog,mlog)
  F.trans = function(x) 1/R * integrate(fun, lower=0, upper=x, rel.tol = 0.000000001, subdivisions=10000,stop.on.error = F)$value # cdf of disease transmission: 
  fun2 = function(x) F.trans(x) - (1-prop) 
  x.a.LN = uniroot(fun2,c(0,10),extendInt = "upX",maxiter=10000)$root  # x such that (1-p)% of transmission is due to cases with R < x
  pmf = function(x) dlogpois(floor(x),slog,mlog) 
  p.a = function(x) integrate(pmf, lower=0, upper=x, rel.tol = 0.000000001, subdivisions=10000)$value # F(X-1) = P(X <= x-1) = P(X < x)
  p.a.LN = 1-p.a(x.a.LN) # 1-F(X-1) = P(X>x-1) = P(X >= x)
  return(p.a.LN)
}
# discrete version
prop.resp.POLN.discrete = function(R,s){
  m = R
  v = (s^2)-m
  mlog = log(m^2/sqrt(v + m^2))
  slog = sqrt(log(1+(v/m^2)))
  fun.int.LN = function(u){
    f=c()
    u = seq(0,u)
    for(i in 1:length(u)){f[i] = u[i] * dlogpois(u[i],slog,mlog)}
    return(sum(f))
  } 
  F.trans = function(x) 1/R * fun.int.LN(x)
  d<-seq(0.03,1,0.01)
  prop.LN<-data.frame("Prop"=rep(0,length(d)),"x1"=rep(0,length(d)),"x2"=rep(0,length(d)))
  k<-1
  for (i in d){
    fun1 = function(x) F.trans(x) - (1-i)
    temp.inv <- floor(uniroot(fun1,c(0,3),extendInt = 'upX')$root)
    if (F.trans(temp.inv)>(1-i)){ 
      x1.temp<-temp.inv-1
      x2.temp<-temp.inv
    }
    if (F.trans(temp.inv)<(1-i)){
      x1.temp<-temp.inv
      x2.temp<-temp.inv+1
    }
    if (F.trans(temp.inv)==(1-i)){
      x1.temp<-temp.inv
      x2.temp<-temp.inv
    }
    prop.LN[k,]<-c(i,(1- plogpois(x2.temp,slog,mlog)) ,(1- plogpois(x1.temp,slog,mlog)))
    k<-k+1
  }
  return(prop.LN)
}

## Poisson-WB
prop.resp.POWB = function(R,s,prop){ 
  m = R
  v = (s^2)-m
  fun.k = function(x) v + m^2 - ((m^2*(gamma(1+2/x)))/(gamma(1+1/x)^2))
  k_w = uniroot(fun.k,c(0.1,100))$root
  l_w = m/(gamma(1+1/k_w))
  if(prop==0) return(0)
  fun = function(x) floor(x)*dbullpois(floor(x),k_w,l_w)
  F.trans = function(x) 1/R * integrate(fun, lower=0, upper=x, rel.tol = 0.000000001, subdivisions=10000,stop.on.error = F)$value # cdf of disease transmission: 
  fun2 = function(x) F.trans(x) - (1-prop) 
  x.a.WB = uniroot(fun2,c(0,100),extendInt = "upX",maxiter=10000)$root  # x such that (1-p)% of transmission is due to cases with R < x
  pmf = function(x) dbullpois(floor(x),k_w,l_w) 
  p.a = function(x) integrate(pmf, lower=0, upper=x, rel.tol = 0.000000001, subdivisions=10000)$value # F(X-1) = P(X <= x-1) = P(X < x)
  p.a.WB = 1-p.a(x.a.WB) # 1-F(X-1) = P(X>x-1) = P(X >= x)
  return(p.a.WB)
}
# discrete version
prop.resp.POWB.discrete = function(R,s){
  m = R
  v = (s^2)-m
  fun.k = function(x) v + m^2 - ((m^2*(gamma(1+2/x)))/(gamma(1+1/x)^2))
  k_w = uniroot(fun.k,c(0.1,100))$root
  l_w = m/(gamma(1+1/k_w))
  fun.int.WB = function(u){
    f=c()
    u = seq(0,u)
    for(i in 1:length(u)){f[i] = u[i] * dbullpois(u[i],k_w,l_w)}
    return(sum(f))
  } 
  F.trans = function(x) 1/R * fun.int.WB(x)
  d<-seq(0.03,1,0.01)
  prop.WB<-data.frame("Prop"=rep(0,length(d)),"x1"=rep(0,length(d)),"x2"=rep(0,length(d)))
  k<-1
  for (i in d){
    fun1 = function(x) F.trans(x) - (1-i)
    temp.inv <- floor(uniroot(fun1,c(0,100))$root)
    if (F.trans(temp.inv)>(1-i)){ 
      x1.temp<-temp.inv-1
      x2.temp<-temp.inv
    }
    if (F.trans(temp.inv)<(1-i)){
      x1.temp<-temp.inv
      x2.temp<-temp.inv+1
    }
    if (F.trans(temp.inv)==(1-i)){
      x1.temp<-temp.inv
      x2.temp<-temp.inv
    }
    prop.WB[k,]<-c(i,(1- pbullpois(x2.temp,k_w,l_w)) ,(1- pbullpois(x1.temp,k_w,l_w)))
    k<-k+1
  }
  return(prop.WB)
}

## Poisson-GG
prop.resp.POGG = function(R,a,p,d,prop){ 
  if(prop==0) return(0)
  fun = function(x) floor(x)*dgengampois(floor(x),a,p,d)
  F.trans = function(x) 1/R * integrate(fun, lower=0, upper=x, rel.tol = 0.000000001, subdivisions=10000,stop.on.error = F)$value # cdf of disease transmission: 
  fun2 = function(x) F.trans(x) - (1-prop) 
  x.a.GG = uniroot(fun2,c(0,5),extendInt = "upX")$root  # x such that (1-p)% of transmission is due to cases with R < x
  pmf = function(x) dgengampois(floor(x),a,p,d)
  p.a = function(x) integrate(pmf, lower=0, upper=x, rel.tol = 0.000000001, subdivisions=10000)$value # F(X-1) = P(X <= x-1) = P(X < x)
  p.a.GG = 1-p.a(x.a.GG) # 1-F(X-1) = P(X>x-1) = P(X >= x)
  return(p.a.GG)
}
# discrete version
prop.resp.POGG.discrete = function(R,a,p,d){
  m = R
  fun.int.GG = function(u){
    f=c()
    u = seq(0,u)
    for(i in 1:length(u)){f[i] = u[i] * dgengampois(u[i],a,p,d)}
    return(sum(f))
  } 
  F.trans = function(x) 1/R * fun.int.GG(x)
  pr<-seq(0.03,1,0.01)
  prop.GG<-data.frame("Prop"=rep(0,length(pr)),"x1"=rep(0,length(pr)),"x2"=rep(0,length(pr)))
  k<-1
  for (i in pr){
    fun1 = function(x) F.trans(x) - (1-i)
    temp.inv <- floor(uniroot(fun1,c(0,3),extendInt = "upX")$root)
    if (F.trans(temp.inv)>(1-i)){ 
      x1.temp<-temp.inv-1
      x2.temp<-temp.inv
    }
    if (F.trans(temp.inv)<(1-i)){
      x1.temp<-temp.inv
      x2.temp<-temp.inv+1
    }
    if (F.trans(temp.inv)==(1-i)){
      x1.temp<-temp.inv
      x2.temp<-temp.inv
    }
    prop.GG[k,]<-c(i,(1- pgengampois(x2.temp,a,p,d)) ,(1- pgengampois(x1.temp,a,p,d)))
    k<-k+1
  }
  return(prop.GG)
}

###########################
## Based on Endo et al
## Calculate proportion responsible for X% (prop) of total transmissions
propresponsible=function(R0,k,prop){
  qm1=qnbinom(1-prop,k+1,mu=R0*(k+1)/k) # smallest x such that F(x)=P(X<=x) >= 1-p = 0.2
  remq=1-prop-pnbinom(qm1-1,k+1,mu=R0*(k+1)/k) # 1-prop-sum(dgampois(seq(0,qm1-1),shape+1,rate))
  remx=remq/dnbinom(qm1,k+1,mu=R0*(k+1)/k) # remq/dgampois(qm1,shape+1,rate)
  q=qm1+1
  1-pnbinom(q-1,k,mu=R0)-dnbinom(q,k,mu=R0)*remx # 1-sum(dgampois(seq(0,q-1),shape,rate))-dgampois(q,shape,rate)*remx
}
# proportion of transmission due to a% of most infectious cases
prop.END = function(a,R,k){
  if(a==0) return(0)
  propEndo = function(x) propresponsible(R,k,x) - a
  prop.endo = try(uniroot(propEndo,c(0.0001,1))$root)
  if(class(prop.endo)=="try-error") prop.endo = 1
  return(prop.endo)
}
