
####################################
# Density & distribution functions #
# Poisson mixtures                 #
####################################

## Poisson-Gamma (=NB)
dgampois = function(x,a,b){
  return(((b^a)/(gamma(a)*factorial(x))) * (gamma(x+a)/((1+b)^(x+a))))
}
dgampois <- Vectorize(dgampois)
pgampois = function(q,a,b){
    f=c()
    u = seq(0,q)
    for(i in 1:length(u)){f[i] = dgampois(u[i],a,b)} # a=shape,b=rate
    return(sum(f))
} 

## Poisson-Lognormal
dlogpois = function(x,sigma,mu){
  fun = function(lambda) ifelse(lambda==Inf|x==Inf,0,exp(log(1) - (lfactorial(x) + log(sigma) + log(sqrt(2*pi))) + (x-1)*log(lambda) + (-lambda-((log(lambda)-mu)^2)/(2*sigma^2))))
  integrate(fun,0,Inf)$value}
dlogpois <- Vectorize(dlogpois)
plogpois = function(q,sigma,mu){
    f=c()
    u = seq(0,q)
    for(i in 1:length(u)){f[i] = dlogpois(u[i],sigma,mu)} 
    return(sum(f))
}

## Poisson-Weibull
dbullpois = function(x,p,l){
  fun = function(lambda) ifelse(lambda==Inf|x==Inf,0,exp(log(p)-(p*log(l)) + log(1)-lfactorial(x) + (p-1+x)*log(lambda) + (-(lambda/l)^p)-lambda ))
  integrate(fun,0,Inf)$value}
dbullpois <- Vectorize(dbullpois)
pbullpois = function(q,p,l){
    f=c()
    u = seq(0,q)
    for(i in 1:length(u)){f[i] = dbullpois(u[i],p,l)} 
    return(sum(f))
}

## Poisson-GeneralizedGamma
dgengampois = function(x,a,p,d){
  fun = function(lambda) ifelse(lambda==Inf|lambda==0|x==Inf,0,exp((log(p) - (d*log(a) + lfactorial(x) + lgamma(d/p))) + (x+d-1)*log(lambda) + -(lambda/a)^p-lambda))
  integrate(fun,0,Inf)$value}
dgengampois <- Vectorize(dgengampois)
pgengampois = function(q,a,p,d){
    f=c()
    u = seq(0,q)
    for(i in 1:length(u)){f[i] = dgengampois(u[i],a,p,d)} 
    return(sum(f))
}

