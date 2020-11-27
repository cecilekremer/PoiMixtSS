
###############################################################################################
### Code used to obtain proportion of simulations where the fitted model has the lowest AIC ###
###############################################################################################


B.gam.pois = 1000 # number of simulated datasets

load("Simulation_logpois_gampoisv3.RData") # NB fit
out=matrix(NA,nrow=B.gam.pois,ncol=6)
for(i in 1:B.gam.pois){
  out[i,] = c(as.numeric(res.sim.gam.pois[[i]]))
}
res.sim.gam.pois = out
res.sim.gam.pois = data.frame(res.sim.gam.pois)
names(res.sim.gam.pois) = c("mu","size","sd mu","sd size","aic","loglik")
head(res.sim.gam.pois)

NB.AIC = res.sim.gam.pois$aic

load("Simulation_logpois_logpoisv3.RData") # POLN fit
out=matrix(NA,nrow=B.gam.pois,ncol=6)
for(i in 1:B.gam.pois){
  out[i,] = c(as.numeric(res.sim.log.pois[[i]]))
}
res.sim.log.pois = data.frame(out)
names(res.sim.log.pois) = c("mlog","slog","sd mlog","sd slog","aic","loglik")
head(res.sim.log.pois)

LN.AIC = res.sim.log.pois$aic

load("Simulation_logpois_bullpoisv3.RData") # POWB fit
out=matrix(NA,nrow=B.gam.pois,ncol=6)
for(i in 1:B.gam.pois){
  out[i,] = c(as.numeric(res.sim.bull.pois[[i]]))
}
res.sim.bull.pois = data.frame(out)
names(res.sim.bull.pois) = c("shape","scale","sd shape","sd scale","aic","loglik")
head(res.sim.bull.pois)

WB.AIC = res.sim.bull.pois$aic

load("Simulation_logpois_gengampoisv3.RData") # POGG fit
out=matrix(NA,nrow=B.gam.pois,ncol=5)
for(i in 1:B.gam.pois){
  out[i,] = c(as.numeric(res.sim.gengam.pois[[i]]))
}
res.sim.gengam.pois = data.frame(out)
names(res.sim.gengam.pois) = c("a","p","d","aic","loglik")
head(res.sim.gengam.pois)

GG.AIC = res.sim.gengam.pois$aic

sum(NB.AIC < GG.AIC & NB.AIC < WB.AIC & NB.AIC < LN.AIC)/B.gam.pois
sum(LN.AIC < NB.AIC & LN.AIC < WB.AIC & LN.AIC < GG.AIC)/B.gam.pois
sum(WB.AIC < NB.AIC & WB.AIC < LN.AIC & WB.AIC < GG.AIC)/B.gam.pois
sum(GG.AIC < NB.AIC & GG.AIC < WB.AIC & GG.AIC < LN.AIC)/B.gam.pois


