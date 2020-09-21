## Author: Tom Miller
## Purpose: toy example illustrating the role of theory in ecological prediction
## with emphasis on pitfalls associated with timescale
## Last update: 18 Sept 2020

## load libraries
library(scales)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

## Using simulated data for now (so no model uncertainty) but will likely update
## with a real data set

## deterministic skeleton of population dynamics
discrete.logistic <- function(N,r,K){
  N * (1 + r * (1 - N / K))
}

## generate a time series with simple demographic stochasticity and geometric growth (K=Inf)
time.steps <- 100
N <- rep(NA,time.steps)
## start with an arbitrarily small number of discrete individuals
N[1] <- 10
for(t in 2:time.steps){
  N[t] = rpois(1,discrete.logistic(N[t-1],r=0.05,K=Inf))
}

## save so that this exact analysis can be reproduced (do this once)
#write.csv(data.frame(t=1:time.steps,N=N),"ts_sim_r0.05KInf.csv",row.names = F)
geom_ts<-read.csv("ts_sim_r0.05KInf.csv")
plot(geom_ts$t,log(geom_ts$N),type="b",pch=16,xlab="Time",ylab="log(N)")

## Take the first 10 years of population growth, fit a model, and predict 1, 10, and 100 years ahead
obs_years<-10
plot(geom_ts$t[geom_ts$t<=obs_years],log(geom_ts$N[geom_ts$t<=obs_years])
     ,type="b",pch=16,xlab="Time",ylab="log(N)")

lm(log(N)~t,data=subset(geom_ts,t<=obs_years))

## prep data for Stan
geom_dat <- list(years=obs_years,
                 t = geom_ts$t[geom_ts$t<=obs_years],
                 logN = log(geom_ts$N[geom_ts$t<=obs_years]))
# MCMC parameters
mcmc_pars <- list(
  warmup = 1000, 
  iter = 5000, 
  thin = 3, 
  chains = 3
)
# fit model
geom_fit <- stan(
  file = 'Stan/geom_ts.stan',
  data = geom_dat,
  warmup = mcmc_pars$warmup,
  iter = mcmc_pars$iter,
  thin = mcmc_pars$thin,
  chains = mcmc_pars$chains )
geom_params_post <- rstan::extract(geom_fit,pars=c("logN0","r","sigma"))

## compare posterior of r to true value
plot(density(geom_params_post$r),xlab="r",ylab="Posterior density",main=" ",lwd=3,cex.lab=1.4)
abline(v=0.05,lty=2)

## plot the data and a sampling of the posterior estimates
plot(geom_ts$t[geom_ts$t<=(obs_years+1)],
     c((geom_ts$N[geom_ts$t<=obs_years]),NA),type="b",pch=16,ylim=c(0,60))

## estimate N 1,10, and 100 years ahead of most recent observation 
## using simple linear extrapolation with the expected sampling variance
pred_lin <- data.frame(ahead1 = rep(NA,nrow(geom_params_post$logN0)),
                       ahead10 = rep(NA,nrow(geom_params_post$logN0)),
                       ahead90 = rep(NA,nrow(geom_params_post$logN0)))
for(i in 1:nrow(geom_params_post$logN0)){
  pred_lin$ahead1[i] <- rnorm(1,
                              mean = geom_params_post$logN0[i] + geom_params_post$r[i] * (obs_years+1),
                              sd = geom_params_post$sigma[i])
  pred_lin$ahead10[i] <- rnorm(1,
                              mean = geom_params_post$logN0[i] + geom_params_post$r[i] * (obs_years+10),
                              sd = geom_params_post$sigma[i])
  pred_lin$ahead90[i] <- rnorm(1,
                              mean = geom_params_post$logN0[i] + geom_params_post$r[i] * (obs_years+90),
                              sd = geom_params_post$sigma[i])
}
points((obs_years+1.2),mean(exp(pred_lin$ahead1)),col="black",pch=2,cex=2)
points((obs_years+1.2),median(exp(pred_lin$ahead1)),col="black",pch=3,cex=2)

## what if instead we predicted future abundance as a stochastic realization
## of the geometric model, sampling over uncertainty in lambda and N0?
geom_stoch <- data.frame(ahead1 = rep(NA,nrow(geom_params_post$logN0)),
                       ahead10 = rep(NA,nrow(geom_params_post$logN0)),
                       ahead90 = rep(NA,nrow(geom_params_post$logN0)))
N_pred<-rep(NA,times = time.steps)
N_pred[obs_years] <- geom_ts$N[geom_ts$t==obs_years]
for(i in 1:nrow(geom_params_post$logN0)){
  for(t in (obs_years+1):time.steps){
    N_pred[t] = rpois(1,discrete.logistic(N_pred[t-1],r=geom_params_post$r[i],K=Inf))
  }
  geom_stoch$ahead1[i] <- N_pred[(obs_years+1)]
  geom_stoch$ahead10[i] <- N_pred[(obs_years+10)]
  geom_stoch$ahead90[i] <- N_pred[(obs_years+90)]
}
points((obs_years+.8),mean((geom_stoch$ahead1)),col="blue",pch=2,cex=2)
points((obs_years+.8),median((geom_stoch$ahead1)),col="blue",pch=3,cex=2)

points(obs_years+1,(geom_ts$N[geom_ts$t==(obs_years+1)]),pch=16,col="red",cex=3)


hist(geom_stoch$ahead1);abline(v=geom_ts$N[geom_ts$t==(obs_years+1)],col="red")

## put it all together in one plot
plot(geom_ts$t,log(geom_ts$N),pch=16,type="b",ylim=c(-20,20),
     xlab="Time",ylab="log(N)",cex.lab=1.4,
     col=c(rep("black",obs_years),rep("white",(time.steps-obs_years))))
lines(geom_ts$t,log(geom_ts$N),col="black")
for(i in 1:nrow(geom_stoch)){
  points((obs_years+1),log(geom_stoch$ahead1[i]),col=alpha("cornflowerblue",0.1))
  points((obs_years+1)+1,pred_lin$ahead1[i],col=alpha("tomato",0.1))

  points((obs_years+10),log(geom_stoch$ahead10[i]),col=alpha("cornflowerblue",0.1))
  points((obs_years+10)+1,pred_lin$ahead10[i],col=alpha("tomato",0.1))

  points((obs_years+90),log(geom_stoch$ahead90[i]),col=alpha("cornflowerblue",0.1))
  points((obs_years+90)+1,pred_lin$ahead90[i],col=alpha("tomato",0.1))
}
legend("topleft",legend=c("Observations","Linear regression","Geometric model"),
       pch=c(16,1,1),col=c("black","tomato","cornflowerblue"),bty="n")

par(mfrow=c(3,1))
plot(density((pred_lin$ahead1)),lwd=3,col="tomato",ylim=c(0,1.5),
     main="One year ahead",xlab="log(N)",cex.lab=1.5)
lines(density(log(geom_stoch$ahead1[geom_stoch$ahead1>0])),lwd=3,col="cornflowerblue")
abline(v=log(geom_ts$N[geom_ts$t==(obs_years+1)]),lwd=3)

plot(density((pred_lin$ahead10)),lwd=3,col="tomato",
     main="Ten years ahead",xlab="log(N)",cex.lab=1.5)
lines(density(log(geom_stoch$ahead10[geom_stoch$ahead10>0])),lwd=3,col="cornflowerblue")
abline(v=log(geom_ts$N[geom_ts$t==(obs_years+10)]),lwd=3)

plot(density((pred_lin$ahead90)),lwd=3,col="tomato",ylim=c(0,0.15),
     main="Ninety years ahead",xlab="log(N)",cex.lab=1.5)
lines(density(log(geom_stoch$ahead90[geom_stoch$ahead90>0])),lwd=3,col="cornflowerblue")
abline(v=log(geom_ts$N[geom_ts$t==(obs_years+90)]),lwd=3)

## extinction pr from geometric model
1-mean(geom_stoch$ahead1>0)
1-mean(geom_stoch$ahead10>0)
1-mean(geom_stoch$ahead90>0)


## the log scale masks the zeros (predicted extinctions)
plot(geom_ts$t,(geom_ts$N),pch=16,type="b",
     col=c(rep("gray",obs_years),rep("white",(time.steps-obs_years))))
lines(geom_ts$t,(geom_ts$N),
      col=c(rep("gray",obs_years),rep("gray",(time.steps-obs_years))))
for(i in 1:nrow(geom_stoch)){
  points((obs_years+1),(geom_stoch$ahead1[i]),col=alpha("cornflowerblue",0.1))
  points((obs_years+1)+1,exp(pred_lin$ahead1[i]),col=alpha("tomato",0.1))
  
  points((obs_years+10),(geom_stoch$ahead10[i]),col=alpha("cornflowerblue",0.1))
  points((obs_years+10)+1,exp(pred_lin$ahead10[i]),col=alpha("tomato",0.1))
  
  points((obs_years+90),(geom_stoch$ahead90[i]),col=alpha("cornflowerblue",0.1))
  points((obs_years+90)+1,exp(pred_lin$ahead90[i]),col=alpha("tomato",0.1))
}
